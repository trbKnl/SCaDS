#### README for the simulation study in section 4.1: Recovery of the model parameters under the correct model 

##### Description of folder content:

SPARSE_PCA_wRandomStart.R

- Contains the **old** SCaDS algorithm (usable but not used in the paper) and some auxillary functions (Tucker congruence)

findLasso.R

- Contains the bisection method to tune the lasso parameter  

generateData.R

- Function to generate data under the model with sparse components weights and a common and distinctive block structure also for the component weights

produceBoxplots.R

- Scripts to create the boxplots in the paper (For the raw simulation result .rds objects either run the simulation yourself or send me an email).  

ridgeTuningForConditions

- Folder in which the ridge parameter is tuned for the individual conditions 


simulationStudyScript.R

- Script that executes the simulation study

sparseSCAcppPackage

- Folder to create a package of the SCaDS algorithm in cpp so it can be used in do parallel
- This is the version of the algorithm used for the simulation study and other analysis


##### R version and package versions

R version 3.5.1 (2018-07-02)

combinat 0.0.8
doParallel 1.0.11
Rcpp 0.12.18
elasticnet 1.1



