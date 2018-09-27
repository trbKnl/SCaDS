###############################################
# Model selection simulation study
#################################################
rm(list=ls())



##############################################################################################
# Start Model selection simulation study

#functions for the simulation study


#This function takes the conditions as input 
#and is used in the parallel package to distribute over the cores
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
#is used in the function "simulationStudy"
simRes <- function(noise, sparsity, ridge, scaleDat){
    
    cases <- 100
    variables <- 500
    allCombn <- allCommonSpecific(c(variables/2, variables/2), 3)

    dat <- generateCommonSpecific(x = variables, nx = cases, nfactors = 3, p = noise,
                                   coefFixed = TRUE, sparsity = sparsity)
    dat$X <- scale(dat$X, scale = FALSE)

    if(scaleDat){
        dat$X <- scale(dat$X, center = TRUE)
    }

    res <- rep(NA, length(allCombn)) 
    stdErrorRes <- rep(NA, length(allCombn))

    for(i in 1:length(allCombn)){
        error <- EigenVectorCV2(dat$X, ridge = ridge , lasso = 0, nrFolds = 10,
                                nfactors = 3, fixW = allCombn[[i]],
                                maxItr = 10000, scaleDat = scaleDat)
        print(ridge)
        res[i] <- error$MSE
        stdErrorRes[i] <- error$stdError
    }

    bestModel <- which.min(res)
    withinOneStdError <- res < res[bestModel] + stdErrorRes[bestModel] 
    orderOfModels <- order(res)
    distanceInStdErr <- (res[5] - res[bestModel]) / stdErrorRes[bestModel] 

    return(list(res = res,
                stdErrorRes = stdErrorRes,
                withinOneStdError = withinOneStdError,
                orderOfModels = orderOfModels,
                distanceInStdErr = distanceInStdErr))
}

#initialize conditions for the simulation study

#create the conditions
error <- c(0.05, 0.25, 0.5) #corresponding to an error ratio of 5, 25, and 50%
sparsity <- c("high", "low") #high sparsity = c(.6, .6, .6), low sparsity = c(.2, .52, .52)

conditions <- expand.grid(error, sparsity)
colnames(conditions) <- c("error",  "sparsity")
conditions

#load in the ridge parameters for the simulation study
ridgeVec <- readRDS('../simulation_study_surfLisa/ridgeTuningForConditions/ridgeVec') 
conditions <- cbind(conditions, ridgeVec)
conditions



# Start simulation study
require(doParallel)

#the number of reps per condition
repCondition <- 20 

#seedList
seedVec  <- 1:nrow(conditions)

cores <- detectCores() - 1
cores
#write output of the CV function to a log to monitor the progress of the simulation
cl <- makeCluster(cores, outfile = "./log")
registerDoParallel(cl)

#scaleDat data is TRUE
simres <- foreach(i = 1:nrow(conditions), 
                  .packages = c('sparseSCAcppFunction', 'gtools', 'combinat')) %dopar% {

    source('./SPARSE_PCA_wRandomStart.R', local = TRUE)
    source('./CVfunction.R', local = TRUE)
    simulationStudy(pos = i, repCondition = repCondition,
                    conditions = conditions, seedVec = seedVec, scaleDat = FALSE) 
}


saveRDS(simres, './modelselectionResults.rds')

stopCluster(cl)


######################################################################
######################################################################
# add simulation study to see whether model selection performs well 
# in the low dimensional setting (same functions different parameters)

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
#drastically reduce the number of variables
simRes <- function(noise, sparsity, ridge, scaleDat){
    
    #number of variables and cases changed
    cases <- 195
    variables <- 20
    allCombn <- allCommonSpecific(c(variables/2, variables/2), 3)

    dat <- generateCommonSpecific(x = variables, nx = cases, nfactors = 3, p = noise,
                                   coefFixed = TRUE, sparsity = sparsity)
    dat$X <- scale(dat$X, scale = FALSE)

    if(scaleDat){
        dat$X <- scale(dat$X, center = TRUE)
    }

    res <- rep(NA, length(allCombn)) 
    stdErrorRes <- rep(NA, length(allCombn))

    for(i in 1:length(allCombn)){
        error <- EigenVectorCV2(dat$X, ridge = 0, lasso = 0, nrFolds = 10,
                                nfactors = 3, fixW = allCombn[[i]],
                                maxItr = 10000, scaleDat = scaleDat)
        res[i] <- error$MSE
        stdErrorRes[i] <- error$stdError
    }

    bestModel <- which.min(res)
    withinOneStdError <- res < res[bestModel] + stdErrorRes[bestModel] 
    orderOfModels <- order(res)
    distanceInStdErr <- (res[5] - res[bestModel]) / stdErrorRes[bestModel] 

    return(list(res = res,
                stdErrorRes = stdErrorRes,
                withinOneStdError = withinOneStdError,
                orderOfModels = orderOfModels,
                distanceInStdErr = distanceInStdErr))
}



#create the conditions
error <- c(0.05, 0.20, 0.50) #corresponding to an error ratio of 5, 25, and 50%
sparsity <- c("high", "low") #high sparsity = c(.6, .6, .6), low sparsity = c(.2, .52, .52)

conditions <- expand.grid(error, sparsity)
colnames(conditions) <- c("error",  "sparsity")
conditions

#load in the ridge parameters for the simulation study
#ridge is set to zero
ridgeVec <- rep(0, 6)
conditions <- cbind(conditions, ridgeVec)
conditions


# Start simulation study
require(doParallel)

#the number of reps per condition
repCondition <- 20 


#seedList
seedVec  <- 1:6
cores <- 3
cores

cl <- makeCluster(cores, outfile = "./log")
registerDoParallel(cl)

#scaleDat data is FALSE
simResGood <- foreach(i = 1:nrow(conditions), 
                  .packages = c('sparseSCAcppFunction', 'gtools', 'combinat')) %dopar% {

    source('./SPARSE_PCA_wRandomStart.R', local = TRUE)
    source('./CVfunction.R', local = TRUE)
    simulationStudy(pos = i, repCondition = repCondition,
                    conditions = conditions, seedVec = seedVec, scaleDat = FALSE) 
}

saveRDS(simResGood, './modelselectionResults_Good.rds')

stopCluster(cl)



