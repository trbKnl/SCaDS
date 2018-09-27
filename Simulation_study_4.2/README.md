#### README file for the model selection simulation

##### Description of files in folder:

CVfunction.R 

- Contains the function to generated data 
- Function to generate all possible common and distinctive component weight matrices
- Crossvalidation function with the EigenVector method 

SPARSE_PCA_wRandomStart.R

- Contains old version of the algorithm
- and auxillary functions (Tucker congruence function, correctly classified function

cppTesting

- contains the C++ version of the SCaDS algorithm

modelselection.R

- Contains simulation study for model selection in the high and low dimensional case

modelselectionResults.rds

- Robject containing the results of the high dimensional condition (This file is small so it is included here).

modelselectionResults_Good.rds

- Results of the low dimensional condition (This file is small so it is included here).

produceREsultsForPaper.R

- Make table and plot for high dimensional condition

produceREsultsForPaper_good.R

- Make plot and table for the lowdimensional condition

ridgeVec

- Contains ridge parameters tuned for the used conditions in the high dimensional case

##### Version of R and packages used

- R version 3.5.1 (2018-07-02)
- gtools 3.8.1
- combinat 0.0.8
