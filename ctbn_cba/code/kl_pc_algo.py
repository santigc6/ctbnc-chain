import pandas as pd
import numpy as np
#import numba
#from numba import jit, prange
import json
import logging
import time
from itertools import combinations
from scipy.stats import f as f_dist
from scipy.stats import chi2 as chi2_dist
from tqdm import tqdm
import sys
#logging.basicConfig(level=logging.DEBUG)

class Ctbn_cb:
    def prepare_trajectories(self,trajectories, variables):
        dimensions = np.array([x.shape[0] - 1 for x in trajectories], dtype=np.int)
        ret_array = np.zeros([dimensions.sum(),trajectories[0].shape[1]*2])
        cum_dim = np.zeros(len(trajectories)+1, dtype=np.int)
        cum_dim[1:] = dimensions.cumsum()
        dimensions.cumsum()
        for it in range(len(trajectories)):
            tmp = trajectories[it].to_numpy()
            dim = tmp.shape[1]
            ret_array[cum_dim[it]:cum_dim[it+1],0:dim] = tmp[:-1]
            ret_array[cum_dim[it]:cum_dim[it+1],dim] = np.diff(tmp[:,0])
            ret_array[cum_dim[it]:cum_dim[it+1], dim+1:] = np.roll(tmp[:,1:],-1, axis=0)[:-1]
        self.trajectories = ret_array
        self.variables = variables
    
    @staticmethod
#    @jit(nopython=True, parallel=False, boundscheck=False, debug=False, cache=True)
    def _compute_cim(trajectories, child_id, parents_id, T_vector, M_vector, parents_comb, M, T):
        diag_indices = np.array([x*M.shape[1]+x%M.shape[1] for x in range(M.shape[0]*M.shape[1])],dtype=np.int64)
        T_filter = np.array([child_id,*parents_id], dtype=np.int) + 1
        T[:] = np.bincount(np.sum(trajectories[:,T_filter]*T_vector/T_vector[0],axis=1).astype(np.int),\
                           trajectories[:,int(trajectories.shape[1]/2)], minlength=T_vector[-1]).reshape(-1,T.shape[1])

        trj_tmp = trajectories[trajectories[:,int(trajectories.shape[1]/2)+1+child_id].astype(np.int)>=0]
        M_filter = np.array([child_id, child_id, *parents_id], dtype=np.int) + 1
        M_filter[0] += int(trj_tmp.shape[1]/2)

        M[:] = np.bincount(np.sum(trj_tmp[:,M_filter]*M_vector/M_vector[0],axis=1).astype(np.int),\
                            minlength=M_vector[-1]).reshape(-1,M.shape[1], M.shape[2])
        M_raveled = M.ravel()
        M_raveled[diag_indices] = 0
        M_raveled[diag_indices] = np.sum(M,axis=2).ravel()

        
        #for comb_idx in prange(parents_comb.shape[0]):
        #    comb = parents_comb[comb_idx]
        #    trj_tmp1 = trajectories
        #    if parents_id.shape[0] > 0:
        #        trj_tmp1 = trj_tmp1[np.all(trj_tmp1[:,parents_id+1] == comb, axis=1)]

            #T[comb_idx] = np.bincount(trj_tmp1[:,child_id+1].astype(np.int), trj_tmp1[:,int(trj_tmp1.shape[1]/2)])
            
        #    trj_tmp1 = trj_tmp1[trj_tmp1[:,int(trj_tmp1.shape[1]/2)+1+child_id].astype(np.int)>=0]
        #    M[comb_idx] = np.bincount(T.shape[1]*trj_tmp1[:,child_id+1].astype(np.int)+\
        #                                         trj_tmp1[:,int(trj_tmp1.shape[1]/2)+1+child_id].astype(np.int),\
        #                                         minlength=np.power(T.shape[1],2)).reshape(-1,T.shape[1])
        #    np.fill_diagonal(M[comb_idx],0)
        #    np.fill_diagonal(M[comb_idx],np.sum(M[comb_idx],axis=1))

        q = (M.ravel()[diag_indices].reshape(-1,M.shape[1])+1)/(T+1)
        theta = (M + 1)/(M.ravel()[diag_indices].reshape(-1,M.shape[2],1) + 1 )
        negate_main_diag = np.ones((M.shape[1],M.shape[2]))
        np.fill_diagonal(negate_main_diag,-1)
        theta = np.multiply(theta,negate_main_diag)
        return theta * q.reshape(-1,M.shape[2],1)

    def compute_cim(self, child_id, parents_id):
        tmp = []
        child_id = int(child_id)
        parents_id = np.array(parents_id, dtype=np.int)
        parents_id.sort()
        for idx in parents_id:
            tmp.append([x for x in range(self.variables.loc[idx,"Value"])])
        if len(parents_id) > 0:
            parents_comb = np.array(np.meshgrid(*tmp)).T.reshape(-1,len(parents_id))
            if len(parents_id) > 1:
                tmp_comb = parents_comb[:,1].copy()
                parents_comb[:,1] = parents_comb[:,0].copy()
                parents_comb[:,0] = tmp_comb
        else:
            parents_comb = np.array([[]], dtype=np.int)
        M = np.zeros([max(1,parents_comb.shape[0]),\
                     self.variables.loc[child_id,"Value"],\
                     self.variables.loc[child_id,"Value"]], dtype=np.int)

        T = np.zeros([max(1,parents_comb.shape[0]),\
                     self.variables.loc[child_id, "Value"]], dtype=np.float)

        T_vector = np.array([self.variables.iloc[child_id,1].astype(np.int)])
        T_vector = np.append(T_vector, [self.variables.iloc[x,1] for x in parents_id])
        T_vector = T_vector.cumprod().astype(np.int)


        M_vector = np.array([self.variables.iloc[child_id,1], self.variables.iloc[child_id,1].astype(np.int)])
        M_vector = np.append(M_vector, [self.variables.iloc[x,1] for x in parents_id])
        M_vector = M_vector.cumprod().astype(np.int)


        CIM = self._compute_cim(self.trajectories, child_id, parents_id, T_vector, M_vector, parents_comb, M, T)
        return parents_comb, M,T,CIM
 #       self._compute_cim.parallel_diagnostics(level=4)
    
    
    def independence_test(self, to_var, from_var, sep_set, alpha):
        parents = np.array(sep_set)
        parents = np.append(parents,from_var)
        parents.sort()
        parents_no_from_mask = parents != from_var
        
        parents_comb_from, M_from, T_from, CIM_from = self.compute_cim(to_var, parents)
        
        
        parents_comb, M, T, CIM = self.compute_cim(to_var, parents[parents_no_from_mask])

        to_card = self.variables.loc[to_var, "Value"]
        diag_imatrix = np.eye(to_card,dtype = np.bool)
        

        for comb_id in range(parents_comb.shape[0]):

            total_t = np.sum(T[comb_id])
            #Bad code, inefficient
            if parents.shape[0] > 1:
                tmp_parents_comb_from_ids = np.argwhere(np.all(parents_comb_from[:,parents_no_from_mask] == parents_comb[comb_id],axis=1)).ravel()
            else:
                tmp_parents_comb_from_ids = np.array([x for x in range(parents_comb_from.shape[0])])

            Q = CIM[comb_id]
            ckl = 0.0
            for comb_from_id in tmp_parents_comb_from_ids:
                Q_from = CIM_from[comb_from_id]
                dkl = Q_from[diag_imatrix]-Q[diag_imatrix]\
                        + np.sum(Q_from[~diag_imatrix].reshape(to_card,-1) \
                                    * np.log(Q_from[~diag_imatrix].reshape(to_card,-1)\
                                            /Q[~diag_imatrix].reshape(to_card,-1)),axis=1)
                ckl += np.sum(T_from[comb_from_id]) * np.sum(T_from[comb_from_id]*dkl)
            ckl /= (total_t**2)
            if 0.5 * (1 + np.sqrt(1 - np.exp(-2*ckl))) > alpha:
                return False

                


        return True

    def cb_structure_algo(self, alpha):
        adj_matrix = np.ones((self.variables.shape[0], self.variables.shape[0]), dtype=np.bool)
        np.fill_diagonal(adj_matrix, False)
        for to_var in tqdm(range(self.variables.shape[0])):
            n = 0
            tested_variables = np.argwhere(adj_matrix[:,to_var]).ravel()
            while n < tested_variables.shape[0]:
                for from_var in tested_variables:
                    if from_var not in tested_variables:
                        continue
                    if n >= tested_variables.shape[0]:
                        break
                    sep_set_vars = tested_variables[tested_variables!=from_var]
                    for comb in combinations(sep_set_vars,n):
                        if self.independence_test(to_var, from_var, comb,alpha):
                            adj_matrix[from_var, to_var] = False
                            tested_variables = np.argwhere(adj_matrix[:,to_var]).ravel()
                            break
                n += 1
        self.structure = adj_matrix

    

def confusion_matrix(real,predicted):
    p_minus_t = predicted - real
    
    fp = np.sum(p_minus_t > 0)
    fn = np.sum(p_minus_t < 0)

    tp = np.sum(real==1) - fn
    tn = np.sum(real==0) - fp - real.shape[0]
    #ret = [{"Edge":tp, "Non-Edge":fp, "_row":"Edge"},\
    #        {"Edge":fn, "Non-Edge":tn, "_row":"Non-Edge"}]
    ret = np.array([[tp,fp],[fn,tn]], dtype=np.int)

    return ret

def adj_list_to_adj_matrix(adj_list,variables):
    adj_matrix = np.zeros((variables.shape[0], variables.shape[0]))
    for edge in adj_list:
        adj_matrix[variables[variables["Name"] == edge["From"]].index[0],\
                    variables[variables["Name"] == edge["To"]].index[0]] = 1
    return adj_matrix




if __name__=="__main__":
    with open(sys.argv[1]) as f:
        data = json.load(f)
    print(data["cardinality"])
    ret_object = {}
    for dataset_path in data["datasets_path"]:
        print(dataset_path)
        dataset_ret = {}
        with open(dataset_path) as f:
            networks = json.load(f)

        for subsample in data["subsamples"]:
            print("{} trajectories".format(subsample))
            cf_matrix = np.zeros((2,2), dtype=np.int) 
            execution_time = 0.0

            for network in networks:
                ctbn_cb = Ctbn_cb()
                traj = [pd.DataFrame(x) for x in network["samples"]]
                a = time.time()
                ctbn_cb.prepare_trajectories(traj[0:subsample],pd.DataFrame(network["variables"]))
                ctbn_cb.cb_structure_algo(alpha=0.55)
                b = time.time()
                execution_time += b-a
                cf_matrix += confusion_matrix(adj_list_to_adj_matrix(network["dyn.str"], pd.DataFrame(network["variables"])),\
                            ctbn_cb.structure)
            dataset_ret["{}_trajectories".format(subsample)] = {"cf_matrix": [\
                        {"Edge": int(cf_matrix[0,0]), "Non-Edge": int(cf_matrix[0,1]), "_row":"Edge"},
                        {"Edge": int(cf_matrix[1,0]), "Non-Edge": int(cf_matrix[1,1]), "_row":"Non-Edge"}],
                        "execution.time": execution_time/len(networks)}
        ret_object[dataset_path[5:-5]] = dataset_ret

    
    with open("metrics/kl_pc_based_{}.json".format(data["cardinality"]),"w") as f:
        json.dump(ret_object,f)







