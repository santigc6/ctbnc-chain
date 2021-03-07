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
    
    
    def independence_test(self, to_var, from_var, sep_set, alpha_exp, alpha_chi2, thumb_threshold):
        parents = np.array(sep_set)
        parents = np.append(parents,from_var)
        parents.sort()
        parents_no_from_mask = parents != from_var
        
        parents_comb_from, M_from, T_from, CIM_from = self.compute_cim(to_var, parents)
        
        if self.variables.loc[to_var, "Value"] > 2:
            df = (self.variables.loc[to_var, "Value"]  - 1 ) ** 2
            df = df * (self.variables.loc[from_var,"Value"] )
            for v in sep_set:
               df = df * (self.variables.loc[v,"Value"])

            if np.all(np.sum(np.diagonal(M_from,axis1=1,axis2=2),axis=1)/df <  thumb_threshold): 
                return False

            chi_2_quantile = chi2_dist.ppf(1-alpha_chi2,self.variables.loc[to_var, "Value"]  - 1 )

        parents_comb, M, T, CIM = self.compute_cim(to_var, parents[parents_no_from_mask])
        
        
        for comb_id in range(parents_comb.shape[0]):
            #Bad code, inefficient
            if parents.shape[0] > 1:
                tmp_parents_comb_from_ids = np.argwhere(np.all(parents_comb_from[:,parents_no_from_mask] == parents_comb[comb_id],axis=1)).ravel()
            else:
                tmp_parents_comb_from_ids = np.array([x for x in range(parents_comb_from.shape[0])])

            for comb_from_id in tmp_parents_comb_from_ids:
                diag = np.diag(CIM[comb_id])
                diag_from = np.diag(CIM_from[comb_from_id])
                r1 = np.diag(M[comb_id])
                r2 = np.diag(M_from[comb_from_id])
                stats = diag_from/diag
                for id_diag in range(diag.shape[0]):
                    if stats[id_diag] < f_dist.ppf(alpha_exp/2, r1[id_diag], r2[id_diag]) or\
                        stats[id_diag] > f_dist.ppf(1-alpha_exp/2, r1[id_diag], r2[id_diag]):
                        return False
 
                if diag.shape[0] > 2:

                    # https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/chi2samp.htm
                    K_from = np.sqrt(M[comb_id].diagonal() / M_from[comb_from_id].diagonal())
                    K = np.sqrt(M_from[comb_from_id].diagonal() / M[comb_id].diagonal())
                    
                    M_no_diag = M[comb_id][~np.eye(diag.shape[0], dtype=np.bool)].reshape(diag.shape[0],-1)
                    M_from_no_diag = M_from[comb_from_id][~np.eye(diag.shape[0], dtype=np.bool)].reshape(diag.shape[0],-1)

                    chi_stats = np.sum((np.power((M_no_diag.T * K).T - (M_from_no_diag.T * K_from).T, 2)\
                                /(M_no_diag + M_from_no_diag)), axis=1)                

                    if np.any(chi_stats > chi_2_quantile):
                        return False

        return True

    def cb_structure_algo(self, alpha_exp, alpha_chi2, thumb_threshold):
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
                        if self.independence_test(to_var, from_var, comb,alpha_exp, alpha_chi2, thumb_threshold):
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
                ctbn_cb.cb_structure_algo(alpha_exp=0.1, alpha_chi2=0.1, thumb_threshold=25)
                b = time.time()
                execution_time += b-a
                cf_matrix += confusion_matrix(adj_list_to_adj_matrix(network["dyn.str"], pd.DataFrame(network["variables"])),\
                            ctbn_cb.structure)
            dataset_ret["{}_trajectories".format(subsample)] = {"cf_matrix": [\
                        {"Edge": int(cf_matrix[0,0]), "Non-Edge": int(cf_matrix[0,1]), "_row":"Edge"},
                        {"Edge": int(cf_matrix[1,0]), "Non-Edge": int(cf_matrix[1,1]), "_row":"Non-Edge"}],
                        "execution.time": execution_time/len(networks)}
        ret_object[dataset_path[5:-5]] = dataset_ret

    
    with open("metrics/exp_and_chi2_test_pc_based_{}.json".format(data["cardinality"]),"w") as f:
        json.dump(ret_object,f)







