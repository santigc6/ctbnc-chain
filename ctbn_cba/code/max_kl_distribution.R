library(igraph)
library(jsonlite)
library(dplyr)
library(parallel)
source("code/utils_algo.R")
ck_dependency_max_kl <- function(trjs,variables, to,from,sep_set){
  cimsout_trjs = list()
  cimsout_trjs_without_to = list()
  for(i in 1:length(trjs)){
    cimsout_trjs[[i]] = trjs[[i]][c("Time",to,from,sep_set)]
    cimsout_trjs_without_to[[i]] = trjs[[i]][c("Time",to,sep_set)]
  }
  cimsout = compute_cim_for_var(cimsout_trjs,variables[variables$Name %in% c(to,from,sep_set),],to)
  cimsout_without_to = compute_cim_for_var(cimsout_trjs_without_to,variables[variables$Name %in% c(to,sep_set),],to)
  comb = list()
  for(y in sep_set)
    comb = c(comb,list(0:(variables[variables$Name==y,"Value"]-1)))
  comb = expand.grid(comb)
  Ckl = 0
  for(i in 1:nrow(comb)){
    if(i == 0)
      break()
    molt_var = 1
    molt_var2 = 0
    molt_var_q0 = 1
    index_value = 0
    index_value_q0 = 0
    index = 1
    for(s in unlist(variables[variables$Name %in% c(from,sep_set),"Name"])){
      if(s == from){
        molt_var2 = molt_var
        molt_var = molt_var * variables[variables$Name==s,"Value"]
      }
      else{
        index_value = index_value + molt_var*comb[i, index]
        index_value_q0 = index_value_q0 + molt_var_q0*comb[i, index]
        index = index + 1
        molt_var = molt_var * variables[variables$Name==s,"Value"]
        molt_var_q0 = molt_var_q0 * variables[variables$Name==s,"Value"]
      }
    }
    PQ0 = cimsout_without_to[[to]][[index_value_q0+1]]
    for(x0 in 0:(variables[variables$Name==from,"Value"]-1)){
      PQ = cimsout[[to]][[index_value+x0*molt_var2+1]]
      Ckl = max(Ckl, kl_max_distance(PQ,PQ0))
    }
  }
  return(Ckl)  
}
Ckl_distribution <- function(trjs,variables,dyn.str){
  from_ret = c()
  to_ret = c()
  sep_set_ret = c()
  n_ret = c()
  ck_ret = c()
  is_edge_ret = c()
  real_edge_ret = c()
  vars = variables$Name
  net.struct = graph_from_data_frame(dyn.str)
  n = 0
  # Variables wich are already completly explored (from the entering edges point of view)
  exausted_vars = c()
  while(length(vars) > length(exausted_vars)){
    for(to in setdiff(vars,exausted_vars)){
      from_variables = vars[vars!=to] 
      if(length(from_variables) <= n){
        exausted_vars = c(exausted_vars, to)
        next()
      }
      for(from in from_variables){
        if(length(from_variables) <= n){
          exausted_vars = c(exausted_vars, to)
          break()
        }
        #print(paste(from_variables," from:",from," n:",n))
        sep_set_comb = t(combn(setdiff(from_variables,c(from)),n))
        for(index_set in 1:nrow(sep_set_comb)){
          if(length(sep_set_comb[index_set,]) == 0)
            sep_set = c()
          else
            sep_set = unlist(sep_set_comb[index_set,])

	  from_ret = c(from_ret, from)
	  to_ret = c(to_ret,to)
	  sep_set_ret = c(sep_set_ret,toString(sep_set))
	  n_ret = c(n_ret, n)
	  ck_ret = c(ck_ret,ck_dependency_max_kl(trjs, variables, to,from,sep_set))
	  real_edge_ret = c(real_edge_ret,are.connected(net.struct,from,to))
	  asp = all_simple_paths(net.struct,from,to)
	  is_edge = FALSE
	  for(sp in asp){
	    if(length(intersect(as_ids(sp),sep_set)) == 0){
	    	is_edge = TRUE
	    	break()
	    }
	    
	  }
	  is_edge_ret = c(is_edge_ret, is_edge)
        }
      }
    }
    n = n + 1
  }

  return(data.frame(From=from_ret, To=to_ret, n=n_ret, sep_set=sep_set_ret,ck=ck_ret,is_edge=is_edge_ret,real_edge=real_edge_ret,stringsAsFactors = FALSE))
  
}


distribution_experiment <- function(dataset_path){
  ret = list()

  dataset_name = substr(dataset_path, 6, nchar(dataset_path)-6)
  print(dataset_name)
  load(dataset_path)
  for(i in seq_along(networks)){
    print(paste(dataset_name, " - ",i,"/",length(networks)))
    samples =  networks[[i]][["samples"]]
    variables = networks[[i]][["variables"]]

    dyn.str = networks[[i]][["dyn.str"]]
    ret[[i]] = list()
    for(subsample in c(10,20,50,100)){
	tmp_samples = samples[1:subsample]
    	ret[[i]][[paste(subsample,"_trajectories",sep="")]] = Ckl_distribution(tmp_samples,variables,dyn.str)
    }
  }
  tmp_ret = list()
  tmp_ret[[dataset_name]] = ret
  return(tmp_ret)
}

distribution_experiments <- function(datasets_path){	
  #distribution_experiment(datasets_path[[1]])
  results = mclapply(datasets_path, distribution_experiment, mc.cores=8)
  results = unlist(results, recursive=FALSE)
  save(results,file="data/Ckl_distribution.RData")
  rj = toJSON(results)
  write(rj, "data/Ckl_distribution.json")
}

datasets_path =c()
for(i in list(3,4,5,6)){
	datasets_path = c(datasets_path, paste("data/networks_and_trajectories_binary_data_",i,".RData",sep=""))
	datasets_path = c(datasets_path, paste("data/networks_and_trajectories_ternary_data_",i,".RData",sep=""))
}

distribution_experiments(datasets_path)






