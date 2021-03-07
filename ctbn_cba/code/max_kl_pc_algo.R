source("code/utils_algo.R")
source("code/experiment_utils.R")

ck_dependency_max_kl <- function(trjs,variables, to,from,sep_set,kx){
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
    Ckl = 0
    PQ0 = cimsout_without_to[[to]][[index_value_q0+1]]
    for(x0 in 0:(variables[variables$Name==from,"Value"]-1)){
      dkl = 0
      PQ = cimsout[[to]][[index_value+x0*molt_var2+1]]
      Ckl = kl_max_distance(PQ,PQ0)
      if(0.5*(1+sqrt(1-exp(-2*Ckl))) > kx){
        return(TRUE)
      }
    }
  }
  return(FALSE)
  
  
}
max_kl_pc_algo <- function(trjs,variables,kx_list=0.6){
  vars = variables$Name
  net_struct = expand.grid(rep(list(vars),2),stringsAsFactors=FALSE)
  colnames(net_struct) = c("From", "To")
  net_struct = net_struct[!net_struct$From == net_struct$To,]
  row.names(net_struct) <- NULL
  
  trjs_with_diff = compute_samples_with_diff(trjs)
  n = 0
  # Variables wich are already completly explored (from the entering edges point of view)
  exausted_vars = c()
  while(length(vars) > length(exausted_vars)){
    for(to in setdiff(vars,exausted_vars)){
      from_variables = net_struct[net_struct$To==to,"From"]
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
	kx = kx_list[[1]] +  kx_list[[2]]*n
        for(index_set in 1:nrow(sep_set_comb)){
          if(length(sep_set_comb[index_set,]) == 0)
            sep_set = c()
          else
            sep_set = unlist(sep_set_comb[index_set,])
          if(!ck_dependency_max_kl(trjs_with_diff, variables, to,from,sep_set,kx)){
            net_struct = net_struct[!(net_struct$From==from & net_struct$To==to),]
            from_variables = net_struct[net_struct$To==to,"From"]
            
          }
        }
      }
    }
    n = n + 1
  }
  return(net_struct)
  
}

args = commandArgs(trailingOnly=TRUE)
source(args[1])
algo_experiments(datasets_path, max_kl_pc_algo, "max_kl_pc_based",cardinality_data,subsamples,list(0.53,0.005))
