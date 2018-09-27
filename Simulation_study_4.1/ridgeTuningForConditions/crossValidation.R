###########################################################################################
# crossvalidation function 
###########################################################################################


# Eigenvector crossvalidation method
EigenVectorCV2 <- function(inputDat, ridge, lasso, nrFolds, nfactors,
                           fixW, nScale = 0, maxItr, scaleDat = FALSE){

    folds <- rep_len(1:nrFolds, nrow(inputDat))
    #folds <- sample(folds)
    cvError  <- matrix(NA, nrow(inputDat), ncol(inputDat))
    MSEkthFold <- rep(NA, nrFolds)

    for(a in 1:nrFolds){

        # actual split of the data
        fold <- which(folds == a)
        trainDat <- inputDat[-fold, ]
        testDat <- inputDat[fold, , drop = FALSE]

        if(scaleDat){
            #Scale training and the testing data
            means <- apply( trainDat, 2, mean )		
            trainDat <- t( t(trainDat) - means)
            sdTrain <- apply(trainDat, 2, function(x) sqrt(sum( x^2, na.rm = TRUE ) / (length(na.omit(x)) - nScale )))
            trainDat <- t( t(trainDat)/sdTrain )
            
            testDat <-  t( t(testDat) - means)
            testDat <- t( t(testDat)/sdTrain ) 
            
            rm(means, sdTrain)		
        } 

        #model estimation
        res <- sparseSCAcpp(trainDat, Q = 3, RIDGE = ridge, LASSO = rep(lasso , 3),
                            fixW = fixW, maxItrOuterloop = maxItr, nStarts = 1,
                            print = FALSE, tol = 10^-8)


        #Eigenvector crossvalidation Bro Kiers
        pred <- matrix(NA, nrow(testDat), ncol(testDat))
        print(a)
        for(i in 1:ncol(inputDat)){
            TMinusVariableJ <- testDat[, -i] %*% res$W[-i, ]
            pred[, i] <- TMinusVariableJ %*% res$P[i, ] 
            }
        cvError[fold, ] <- (testDat - pred)^2
        MSEkthFold[a]  <-  mean(cvError[fold, ]) 
        }

    return(list(cvError = cvError, MSE = mean(MSEkthFold), MSEkthFold = MSEkthFold, 
                stdError = sd(MSEkthFold) / sqrt(nrFolds)))
}


