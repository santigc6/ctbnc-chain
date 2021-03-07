library(ctbn)
library(dplyr)

kl_max_distance <- function(m_a,m_b){
  m_a = as.matrix(m_a)
  m_b = as.matrix(m_b)
  kl = diag(m_a)- diag(m_b)
  tmp = log(m_a/m_b)
  diag(tmp)=0
  diag(m_a) = 0
  kl = kl + rowSums(m_a*tmp)
  return(max(kl))
}

compute_cim_for_var <- function(trjs,variables,to){
  vars = variables$Name
  From = c()
  To = c()
  for(var in vars){
    if(var != to){
      From = c(From,var)
      To = c(To,to)
    }
  }
  if(length(From)>0)
    struct = data.frame(From,To,stringsAsFactors = FALSE)
  else
    struct = data.frame(From=character(), To=character())
  tmp_ctbn = NewCtbn(variables)
  SetDynStruct(tmp_ctbn,struct)
  LearnCtbnParams(tmp_ctbn,trjs,inf.type="exact")
  cimsOut  <- GetDynIntMats(tmp_ctbn)
  tmp_ctbn <- DeleteCtbn(tmp_ctbn)
  garbage <- gc()
  
  return(cimsOut)
}

compute_samples_with_diff <- function(samples){
  ret = list()
  for(x in 1:length(samples)){
    s = samples[[x]]
    s[1:nrow(s)-1,"Diff"] = diff(s$Time)
    ret[[x]] = s
  }
  return(ret)
}



