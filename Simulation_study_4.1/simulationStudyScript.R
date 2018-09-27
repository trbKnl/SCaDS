############
# Begin simulation study

rm(list=ls())

require(devtools)
require(doParallel)
require(Rcpp)
require(elasticnet)
source('./SPARSE_PCA_wRandomStart.R')
source('./findLasso.R')
source('./generateData.R')
sourceCpp('./sparseSCA.cpp')


##############################################################################################
#functions for the simulation study

#function that distributes the simulation conditions
simulationStudy <- function(pos, repCondition, conditions, seedVec, scaleDat){
  
  set.seed(seedVec[pos])
  
  cat("i = ", pos, fill = TRUE)
  conditionRes <- vector('list', repCondition)

  noise <- conditions[pos, 1]
  sparsity <- conditions[pos, 2]
  ridge <- conditions[pos, 3]

  if( sparsity == "high"){
      sparsity  <- c(0.6, 0.6, 0.6)
  } else {
      sparsity  <- c(0.02, 0.52, .52)
  }
  
  for(j in 1:repCondition){
    conditionRes[[j]] <- tryCatch({
      simRes(noise,  sparsity, ridge, scaleDat)
    },
    error = function(cond){
      message("Here's the original error message:")
      message(cond)
      return(NA)
    }
    )
  }
  
  return(conditionRes)
} 

#given a set of conditions, this function generate data and estimates the parameters
simRes <- function(noise, sparsity, ridge, scaleDat = TRUE){
    
    cases <- 100
    variables <- 500

    dat <- generateCommonSpecific(x = variables, nx = cases, nfactors = 3, p = noise,
                                   coefFixed = TRUE, sparsity = sparsity)
    dat$X <- scale(dat$X, scale = FALSE)

    if(scaleDat){
        dat$X <- scale(dat$X, center = TRUE)
    }

      lasso <- findLasso(dat, maxItr = 1000, ridge)

      resSparseSCA <- sparseSCAcpp(dat$X, Q = 3, RIDGE = ridge, LASSO = rep(lasso$lasso, 3),
                                   fixW = dat$fixW, maxItrOuterloop = 100000, nStarts = 1, 
                                   print = TRUE, tol = 10^-8)
                       
      amountZeroes <- round(sparsity * variables)
      resElasticnet <- spca(dat$X, K  = 3, para = variables - amountZeroes, type=c("predictor"),
                  sparse="varnum", use.corr=FALSE, lambda = ridge, max.iter = 100000,
                  trace = TRUE, eps.conv = 10^-3)

      correctClassElasticnet <-   correctlyClassified(dat$W, resElasticnet$loadings)
      correctClassSparseSCA <-   correctlyClassified(dat$W, resSparseSCA$W)
      tuckerConElasticnet <- tuckerCongruence(dat$W, resElasticnet$loadings)
      tuckerConSparseSCA <- tuckerCongruence(dat$W, resSparseSCA$W)


      return(list(dat = dat, 
                  correctClassElasticnet = correctClassElasticnet,
                  correctClassSparseSCA = correctClassSparseSCA,
                  tuckerConElasticnet = tuckerConElasticnet,
                  tuckerConSparseSCA = tuckerConSparseSCA,
                  resSparseSCA = resSparseSCA,
                  resElasticnet = resElasticnet))
}


out <- simRes(noise = 0.05, sparsity = c(.6,.6,.6), ridge = .2)
out$correctClassElasticnet
out$correctClassSparseSCA
out$tuckerConElasticnet
out$tuckerConSparseSCA


############################################################################
# create conditions 
error <- c(0.05, 0.25, 0.5) #corresponding to an error ratio of 5, 25, and 50%
sparsity <- c("high", "low") #high sparsity = c(.6, .6, .6), low sparsity = c(.2, .52, .52)

conditions <- expand.grid(error, sparsity)
colnames(conditions) <- c("error",  "sparsity")
conditions

#load in the ridge parameters for the simulation study
ridgeVec <- readRDS('./ridgeTuningForConditions/ridgeVec')
ridgeVec


#optimize script for 12 cores
#conditions <- rbind(conditions, conditions)
#ridgeVec <- c(ridgeVec, ridgeVec)
conditions  <- cbind(conditions, ridgeVec)
conditions

#the number of reps per condition
repCondition <- 20 

#seedList
seedVec  <- 1:nrow(conditions)

################### start simulation ############################
cores <- detectCores() - 1
cores
cl <- makeCluster(cores)
registerDoParallel(cl)


simres <- foreach(i = 1:nrow(conditions), 
                  .packages = c('sparseSCAcppFunction', 'elasticnet', 'combinat')) %dopar% {

    source('./SPARSE_PCA_wRandomStart.R', local = TRUE)
    source('./findLasso.R', local = TRUE)
    source('./generateData.R', local = TRUE)

    simulationStudy(pos = i, repCondition = repCondition,
                    conditions = conditions, seedVec = seedVec, scaleDat = FALSE) 
}

stopCluster(cl)

simres
saveRDS(simres, './simres4.rds')

######################################################################
