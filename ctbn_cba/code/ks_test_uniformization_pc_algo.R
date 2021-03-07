source("code/utils_algo.R")
source("code/experiment_utils.R")


ks_quantile <- function(alpha){
  return(sqrt(-0.5*log(alpha/2)))
}

compute_samples_stats <- function(samples){
  samples_with_diff = list()
  for(x in 1:length(samples)){
    s = samples[[x]]
    tmp_diff = diff(as.matrix(s))
    s[1:nrow(s)-1,"Diff"] = tmp_diff[,"Time"]
    tmp_diff = abs(tmp_diff )>0
    tmp_diff = data.frame(tmp_diff)
    #tmp_diff$Time <- NULL
    for(col_name in colnames(tmp_diff))
      if(col_name != "Time")
	    s[1:nrow(s)-1,paste(col_name,"_diff",sep="")] = tmp_diff[col_name]
    samples_with_diff[[x]] = s
  }
  #print(samples_with_diff)
  return(samples_with_diff)
}

filter_trjs_by_values <- function(samples, variables,values){
  filtered_samples = samples
	for(x in seq(samples)){

		for(var_id in seq(variables)){
            if(nrow(filtered_samples[[x]]) == 0)
              break
			filtered_samples[[x]] = filtered_samples[[x]][filtered_samples[[x]][variables[[var_id]]] == values[[var_id]],]
		}
	}
  return(filtered_samples)
}

dependency_ks_test_uniformization <- function(trjs,variables, to,from,sep_set,alpha){
  cimsout_trjs = list()
  cimsout_trjs_without_from = list()
  for(i in 1:length(trjs)){
    cimsout_trjs[[i]] = trjs[[i]][c("Time",to,from,sep_set)]
    cimsout_trjs_without_from[[i]] = trjs[[i]][c("Time",to,sep_set)]
  }
  cimsout = compute_cim_for_var(cimsout_trjs,variables[variables$Name %in% c(to,from,sep_set),],to)
  cimsout_without_from = compute_cim_for_var(cimsout_trjs_without_from,variables[variables$Name %in% c(to,sep_set),],to)
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
    
    PQ0 = cimsout_without_from[[to]][[index_value_q0+1]]
    if(i <= nrow(comb)){
      filtered_samples = filter_trjs_by_values(trjs, sep_set, comb[i,]) 
    }
    else
      filtered_samples = trjs
    r0_values = list()
    for(x in 1:variables[variables$Name == to, "Value"]){
      tmp_sum = 0
      tmp_filtered = filter_trjs_by_values(filtered_samples, c(to), c(x-1))
      for(s_index in seq(filtered_samples)){
        tmp_sum = tmp_sum + sum(as.matrix(tmp_filtered[[s_index]][paste(to,"_diff",sep="")]))
      }
      r0_values[[x]] = tmp_sum
   }

    for(x0 in 0:(variables[variables$Name==from,"Value"]-1)){
      dkl = 0
      PQ = cimsout[[to]][[index_value+x0*molt_var2+1]]
      alpha_norm = abs(min(min(PQ),min(PQ0)))
      MPQ = PQ/alpha_norm + diag(3)
      MPQ0 = PQ0/alpha_norm + diag(3)

      D_ks = abs(MPQ-MPQ0)
	    for(x in  0:(variables[variables$Name==to, "Value"]-1)){
		  tmp_filtered = filter_trjs_by_values(filtered_samples, c(to, from), c(x, x0))
		  r1 = 0
		  for(s_index in seq(tmp_filtered))
		    r1 = r1 + sum(as.matrix(tmp_filtered[[s_index]][paste(to,"_diff",sep="")]))
          quantile = ks_quantile(alpha)*sqrt((r0_values[[x+1]]+r1)/(r0_values[[x+1]]*r1))
          if(!is.na(max(D_ks[x+1,])) & !is.na(quantile))
		    if(max(D_ks[x+1,])>quantile)
		      return(TRUE)
		}
    }
  }
  return(FALSE)
  
  
}

ks_test_uniformization_pc_algo <- function(trjs,variables,alpha=0.05){
  vars = variables$Name
  net_struct = expand.grid(rep(list(vars),2),stringsAsFactors=FALSE)
  colnames(net_struct) = c("From", "To")
  net_struct = net_struct[!net_struct$From == net_struct$To,]
  row.names(net_struct) <- NULL
  trjs_with_diff  =  compute_samples_stats(trjs)
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
        for(index_set in 1:nrow(sep_set_comb)){
          if(length(sep_set_comb[index_set,]) == 0)
            sep_set = c()
          else
            sep_set = unlist(sep_set_comb[index_set,])
          if(!dependency_ks_test_uniformization(trjs_with_diff, variables, to,from,sep_set,alpha)){
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
algo_experiments(datasets_path, ks_test_uniformization_pc_algo, "ks_test_uniformization_pc_based",cardinality_data,subsamples,0.1)



