library(ctbn)
source("code/experiment_utils.R")

args = commandArgs(trailingOnly=TRUE)

score_based <- function(samples,variables){
    tmp_ctbn = NewCtbn(variables)
    LearnCtbnStruct(tmp_ctbn,samples)
    score.learned.dyn.str = GetDynStruct(tmp_ctbn)
    tmp_ctbn <- DeleteCtbn(tmp_ctbn)
    garbage <- gc()
    return(score.learned.dyn.str)
}

source(args[1])

algo_experiments(datasets_path, score_based,"score_based",cardinality_data,subsamples)



