# Constraint-based Learning for Continuous-Time Bayesian Network

## Dependencies

During this project we used both R and python.

### R dependencies

- [CTBN-RLE](http://rlair.cs.ucr.edu/ctbnrle/Rinterface/)
- [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)

### Python3 dependencies:

- [dvc](https://dvc.org/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://www.scipy.org/)
- [tqdm](https://pypi.org/project/tqdm/)


## DVC

We decided to use the Data Version Control (dvc) software in order to define
a pipeline for each experiment. 

"_DVC is built to make ML models shareable and reproducible._ 
_It is designed to handle large files, data sets, machine learning models,_ 
_and metrics as well as code._" [dvc.org](https://dvc.org/)


## Dataset

We have generated  datasets combining the following parameters:

- Number of nodes: 3,4,5,6,10,15,20
- Network density: 0.1,0.2,0.3
- Node cardinality: 2,3
- Time end:100
- Number of trajectories: 100,200,300

## Reproduce the experiments

During our research we developed and tested multiple algorithms but, for brevity, we decided
to present only the best two:

- ![CTPC{\chi ^2}](https://render.githubusercontent.com/render/math?math=CTPC_%7B%5Cchi%20%5E2%7D)
- ![CTPC{KS}](https://render.githubusercontent.com/render/math?math=CTPC_%7BKS%7D)

We assess the performance of these two algorithms against that of the score-based algorithm implemented
in the CTBN-RLE library.


### Algorithm ![CTPC{\chi ^2}](https://render.githubusercontent.com/render/math?math=CTPC_%7B%5Cchi%20%5E2%7D)

  This algorithm has been implemented in python and the source code can be retrieved from folder _code_. 
  [code/exp\_and\_chi2\_test\_pc\_based\_algo.py](code/exp_and_chi2_test_pc_based_algo.py).
  
  The experiments for ![CTPC{\chi ^2}](https://render.githubusercontent.com/render/math?math=CTPC_%7B%5Cchi%20%5E2%7D)
  can be reproduced using the following commands:

- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_binary_01.dvc`
- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_binary_02.dvc`
- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_binary.dvc`
- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_ternary_01.dvc`
- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_ternary_02.dvc`
- `dvc repro algo_dvc_files/exp_and_chi2_test_pc_based_algo/exp_and_chi2_test_pc_based_algo_ternary.dvc`

### Algorithm ![CTPC{KS}](https://render.githubusercontent.com/render/math?math=CTPC_%7BKS%7D)

  This algorithm has been implemented in python and the source code can be retrieved from folder _code_. 
  [code/exp\_and\_ks\_test\_pc\_based\_algo.py](code/exp_and_ks_test_pc_based_algo.py)
  
  The experiments for ![CTPC{KS}](https://render.githubusercontent.com/render/math?math=CTPC_%7BKS%7D)  
  can be reproduced using the following commands:

- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_binary_01.dvc`
- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_binary_02.dvc`
- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_binary.dvc`
- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_ternary_01.dvc`
- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_ternary_02.dvc`
- `dvc repro algo_dvc_files/exp_and_ks_test_pc_based_algo/exp_and_ks_test_pc_based_algo_ternary.dvc`

### Score based

Score-based learning has been performed by the [CTBN-RLE](http://rlair.cs.ucr.edu/ctbnrle/Rinterface/) package.

The following commands reproduces the  experiments for this algorithm.

- `dvc repro algo_dvc_files/score_based/score_based_binary_01.dvc`
- `dvc repro algo_dvc_files/score_based/score_based_binary_02.dvc`
- `dvc repro algo_dvc_files/score_based/score_based_binary.dvc`
- `dvc repro algo_dvc_files/score_based/score_based_ternary_01.dvc`
- `dvc repro algo_dvc_files/score_based/score_based_ternary_02.dvc`
- `dvc repro algo_dvc_files/score_based/score_based_ternary.dvc`


## Results

All the results can be retrieved from folder  _metrics_ in json format.

For each combination of network density, node cardinality and algorithm there is
a single file structured as follow:

```
  "networks_and_trajectories_ternary_data_01_3": {  <- cardinality, netowrk density and number of nodes.
    "100_trajectories":     <- number of trajectories
    {
      "cf_matrix": [  <- confusion matrix
        {
          "Edge": TP,  <- True Positive
          "Non-Edge": FP <- False Positive
        },
        {
          "Edge": FN, <- False Negative
          "Non-Edge": TN <- True Negative
        ],
      "esecution.time": time <- execution time in seconds
    },
    ...
  }
```

