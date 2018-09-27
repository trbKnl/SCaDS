################################################
# file contains:
# more versatile function to simulate data with
# crossvalidation function
################################################
library(gtools)

#adapt the data generating to be more suitable
#so you can give the sparcity in percentage, and you can set the amount of error
generateCommonSpecific <- function(x, nx, nfactors, p, 
                                   coefFixed = TRUE, sparsity){
    X <- matrix(rnorm(x*nx, 0, 3), nx, x)
    X <- scaleData(X)
    XTX <- t(X) %*% X
    V <- svd(X)$v
    W <- V[, 1:nfactors]
    W[1:floor(nrow(W)/2), 2] <- 0
    W[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0

    if(coefFixed){
       fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
       fixW[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0
       fixW[1:(floor(nrow(W)/2)), 2] <- 0
    } else {
        fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
    }

    #create sparsity in the columns of W based on a percentage
    for(i in 1:ncol(W)){
        spar <- quantile(abs(W[, i]), probs = sparsity[i])
        W[abs(W[, i]) < spar, i] <- 0
    }

    XandP <- betterXandP(X, W)
    Xtrue <- XandP$X
    P <- XandP$P

    E <- matrix(rnorm(x*nx, 0, 1), nx, x)
    g = sqrt((var(as.vector(Xtrue))*p) / (var(as.vector(E)) * (1 - p))) 
    X <- Xtrue + E*g

    SStrue <- var(as.vector(Xtrue))
    SSX <- var(as.vector(X))

    return(list(X = X, W = W, P = P, errorRatio = 1 - SStrue/SSX,
    percentageZeroes = apply(W, 2, countZero), fixW = fixW))
}


betterXandP <- function(X, W){
    for(i in 1:1){
        X <- scaleData(X)
        XTX <- t(X) %*% X
        SVD <- svd(XTX %*% W)   
        P <- SVD$u %*% t(SVD$v)  
        X <- X %*% W %*% t(P)
    }
    return(list(X = X, P = P))
}

#function to get all common and specific block structures
allCommonSpecific <- function(vars, components){

    blocks <- length(vars)

    W  <- matrix(NA, sum(vars), components)
    cd <- as.matrix(expand.grid(rep(list(0:1), blocks))[-1, ])
    commonSpecific <- combinations(n = nrow(cd), r = components,
                                  v = c(1:nrow(cd)), repeats.allowed = TRUE)
    allpossibleWmatrices <- rep(list(NA), nrow(commonSpecific))

    for(i in 1:nrow(commonSpecific)){
        for(j in 1:ncol(commonSpecific)){
            W[, j] <-  rep(cd[commonSpecific[i,j], ], times = vars)
            allpossibleWmatrices[[i]] <- W
    }
    }
    return(allpossibleWmatrices)
}



###########################################################################
# Eigenvector crossvalidation method
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
        res <- sparseSCAcpp(trainDat, Q = nfactors, RIDGE = ridge, LASSO = rep(lasso , nfactors),
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





