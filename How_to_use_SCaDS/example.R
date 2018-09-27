#################################################
# Example of data analysis with SCaDS
#################################################

# Load in needed packages
library(Rcpp)
library(MASS)

source("./CVfunction.R")
sourceCpp("./sparseSCA.cpp")


# Generate our data set X for our example, with 100 objects and 30 variables
n <- 100
x <- 30
X <- mvrnorm(n = 100, mu = rep(0, x), Sigma = diag(x)) 

# To use the SCaDS method we need to pick values for the lasso
# and the ridge parameter. A method of choosing these values is by crossvalidation 
# we implemented cross valdiation with the EigenVector method (see the paper for a reference)


# This cross validation function can be used as follows by specifying the parameters:
# inputDat: the input X 
# ridge: value for the ridge penalty
# lasso: value for the lasso penalty
# nrFolds: number of folds, a sensible default is 10
# nfactors: the number of components
# fixW: The imposed constraints weights (0 indicates a constraint, a non zero value indicates a free weight)
# nScale: biased or unbased estimator of the variance 1 is unbiased
# maxItr: maximum iterations of SCaDS
# scaleDat: can be used if the data is scaled


# Let's do a crossvaludation for a set of parameters
nfactors <- 3
fixW <- matrix(1, nrow = x, ncol = nfactors)

EigenVectorCV2(inputDat = X , ridge = 0.1, lasso = 0.1, nrFolds = 10, nfactors = nfactors, 
                           fixW = fixW, nScale = 0, maxItr = 10000, scaleDat = FALSE)

# The function returns:
# MSE: The mean squared prediction error 
# MSEkthFold: The mean squared prediction error per fold 
# stdError: The standard error of the mean squared prediction error
# cvError: The mean squared prediction error for each individual data point


# To pick the "best" ridge value given the other parameters
# a range of ridge values need to be examined and the best one needs to be picked, This can be done with CV this can be done as follows

ridgeRange <- seq(0, 2, by = 0.01)
MSERidge <- rep(NA, length(ridgeRange)) 
stdErrMSERidge <- rep(NA, length(ridgeRange)) 

# try out all the ridge parameters in ridgeRange
for(i in 1:length(ridgeRange)){
    CVout <- EigenVectorCV2(inputDat = X , ridge = ridgeRange[i], lasso = 0.1, nrFolds = 10,
                            nfactors = nfactors, fixW = fixW, nScale = 0,
                            maxItr = 10000, scaleDat = FALSE)
    MSERidge[i] <- CVout$MSE
    stdErrMSERidge[i] <- CVout$stdError
}

# The same process can be used to pick the lasso value, the number of components, and the common and distinctive structure


# To get all common and distinctive structures for the components weights
# you can use the following function: 
allCommonSpecific(vars, components)

# vars you need to specify the number of variables in each block in a vector:
# c(10, 5, 5, 10), the first block has 10 variables, the second block 5, the third block 5, the fourth block 10 variabes
# c(15, 15) the first block has 15 variables the second block has 15 variables
# Components is the number of components you want to estimate

# Lets say we have 2 blocks where the first 15 variables belong to the first block, and remaning 15 variables belong to the  second block. To get all possible component weight structures do: 

allCommonSpecific(vars = c(15, 15), components = 3)
allCommonSpecific(vars = c(10, 5, 5, 10), components = 3)

# The function returns all possible structures in a list

# To pick the "best" structure according to CV we can repeat the above proces
allCommonDisStructures <- allCommonSpecific(vars = c(15, 15), components = 3)
MSEComDisStruc <- rep(NA, length(allCommonDisStructures)) 
stdErrMSEComDisStruc <- rep(NA, length(allCommonDisStructures)) 

# Try out all the comon and distinctive structures 
for(i in 1:length(allCommonDisStructures)){
    CVout <- EigenVectorCV2(inputDat = X , ridge = 0.1, lasso = 0.1, nrFolds = 10,
                            nfactors = nfactors, fixW = allCommonDisStructures[[i]], nScale = 0,
                            maxItr = 10000, scaleDat = FALSE)
    MSEComDisStruc[i] <- CVout$MSE
    stdErrMSEComDisStruc[i] <- CVout$stdError
}


# If all the parameters have been chosen final analysis can be done with the sparseSCAcpp function:
# The function has the following parameters:

# X: input data
# Q: the number of factors
# RIDGE: the ridge parameter
# LASSO: the lasso parameter vector with a lasso parameter for each component
# fixW: the chosen constraints zero constraints on W
# maxItrOuterloop: number of iterations
# nStarts: the amount of starts (DOES NOT WORK YET), the algorithm uses a warm start based on the first Q columns of V from the SVD of X
# print: print loss function value per iteration 
# tol: the algorithm terminates when the difference between loss functionsvalues of the iterations is less than tol 


# Let's do an analysis
# Pick a structure at random
fixW <- allCommonDisStructures[[5]]
nfactors <- 3

sparseSCAcpp(X = X, Q = nfactors, RIDGE = 0.1, LASSO = rep(0.1, nfactors), 
                       fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
                       print = TRUE, tol = 10^-8)
 
# This function returns:
# W: The component weight matrix
# P: the component loading matrix
# loss: The loss function value upon termination


