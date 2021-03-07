source("code/utils_algo.R")
source("code/experiment_utils.R")

matrix_distance_test_kl <- function(m_a, m_b,alpha){
  kl = max (kl_max_distance(m_a,m_b),kl_max_distance(m_b,m_a))
  return((0.5*(1+sqrt(1-exp(-2*kl))) > alpha))
}

max_kl_algo <- function(trjs,variables,alpha=0.01){
  vars = variables$Name
  alpha = alpha[[1]] + alpha[[2]]*(length(vars)-2)
  From = c()
  To = c()
  for(var in vars){
    cimsout = compute_cim_for_var(trjs,variables,var)
    tmp_vars = vars[!vars==var]
    for(var2 in tmp_vars){
      comb = list()
      for(y in tmp_vars[!tmp_vars==var2]){
        comb = c(comb,list(0:(variables[vars==y,"Value"]-1)))
      }
      comb = expand.grid(comb)
      break_var = FALSE
      for(i in 1:nrow(comb)){
        molt_var = 1
        molt_var2 = 0
        index_value = 0
        index = 1
        for(var3 in tmp_vars){
          if(var3 == var2){
            molt_var2 = molt_var
            molt_var = molt_var * variables[vars==var3,"Value"]
          }
          else{
            index_value = index_value + molt_var*comb[i, index]
            index = index + 1
            molt_var = molt_var * variables[vars==var3,"Value"]
          }
        }
        for(x1 in 0:(variables[vars==var2,"Value"]-2)){
          for(x2 in (x1+1):(variables[vars==var2,"Value"]-1)){
            if(matrix_distance_test_kl(cimsout[[var]][[index_value+x1*molt_var2+1]],cimsout[[var]][[index_value+x2*molt_var2+1]],alpha))
              break_var = TRUE
              break()
          }
          if(break_var)
            break()
        }
        if(break_var)
          break()
      }
      if(break_var){
        From = c(From,var2)
        To = c(To, var)
      }
    }
  }
  ret = data.frame(From,To,stringsAsFactors = FALSE)
  return(ret)
}

args = commandArgs(trailingOnly=TRUE)
source(args[1])
algo_experiments(datasets_path, max_kl_algo, "max_kl_based",cardinality_data,subsamples,list(0.56,0.01))

