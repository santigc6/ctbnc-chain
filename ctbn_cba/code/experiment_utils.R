library(jsonlite)
library(dplyr)

algo_experiment <- function(networks, learn_function,subsamples, ...){
  ret = list()
  for(ss in subsamples){
    print(paste(ss, "trajectories"))
    cf_matrix = data.frame(c(0,0),c(0,0), row.names = c("Edge","Non-Edge"))
    colnames(cf_matrix) = c("Edge", "Non-Edge")
    execution.time=0
    
    
    for(i in seq_along(networks)){
      print(paste(i,"/",length(networks)))
      dyn.str = networks[[i]][["dyn.str"]]
      variables = networks[[i]][["variables"]]
      dyn.cims = networks[[i]][["dyn.cims"]]
      samples =  networks[[i]][["samples"]]
      samples = samples[1:ss]
      
      ptm <- proc.time()
      learned.dyn.str = learn_function(samples,variables,...)
      ptm = proc.time() - ptm
      execution.time = execution.time + ptm
      tp = 0
      if(nrow(learned.dyn.str) > 0)
        tp = nrow(intersect(dyn.str,learned.dyn.str))
      fp = nrow(learned.dyn.str)-tp
      fn = nrow(dyn.str)-tp
      tn = nrow(variables)*(nrow(variables)-1) - nrow(dyn.str) - fp
      cf_matrix["Non-Edge","Non-Edge"] = cf_matrix["Non-Edge","Non-Edge"] + tn
      cf_matrix["Edge","Edge"] = cf_matrix["Edge","Edge"] + tp
      cf_matrix["Edge","Non-Edge"] = cf_matrix["Edge","Non-Edge"] + fp
      cf_matrix["Non-Edge","Edge"] = cf_matrix["Non-Edge","Edge"] + fn
    }
    ret[[paste(ss,"_trajectories",sep="")]] = list(cf_matrix=cf_matrix, execution.time=execution.time[["elapsed"]]/length(networks))
  }
  return(ret)
}

algo_experiments <- function(datasets_path, learn_function,function_name,data_name,subsamples,...){
	results = list()
	for(dataset_path in datasets_path){
		dataset_name = substr(dataset_path, 6, nchar(dataset_path)-6)
		print(dataset_name)
		load(dataset_path)
		results[[dataset_name]] = algo_experiment(networks, learn_function,subsamples,...)
	}
	write(toJSON(results),paste("metrics/",function_name,"_",data_name,".json", sep=""))
	
}
