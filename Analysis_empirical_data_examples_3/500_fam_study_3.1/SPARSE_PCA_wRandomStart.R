##################################################################
# Coordinate descent for sparse PCA
##################################################################

# Function as described in the techincal report

sparsePCAMultistart <- function(data, nfactors, RIDGE, LASSO, 
                      percentage, fixW, maxItrOuterloop, percentageMode = FALSE){
  
  #data preprocessing
  data <- scale(data, scale = FALSE)
  n <- nrow(data)          #determine n
  q <- nfactors            #determine the amount of factors q
  k <- dim(data)[2]        #determine the amount of items k
  XTX <- t(data) %*% data  #X'X 
  XTXdiag <- diag(XTX)

  lossFunctionValue <- c(Inf, rep(Inf, maxItrOuterloop))
  
  #initialize, W
  #for a random W matrix choose:
  #W  <- matrix(rnorm(k*q, 100), k, q)
  W  <- svd(data)$v
  W <- W[, 1:q]    
  W[fixW == 0] <- 0
  WINIT <- W
  
  N <- n

  if(percentageMode == FALSE){
      if(length(LASSO) == 1){
          LASSO <- rep(LASSO[1], q)
      }
  } else {
      LASSO <-  determineLasso(it = 200, XTX, W, LASSO, k, q, N,  percentage, data, RIDGE, WINIT, XTXdiag)
  }

    resultList <- vector('list', 1) 
    minLossVector <- rep(NA, 1)

    for(z in 1:1){

        lossFunctionValue <- c(Inf, rep(Inf, maxItrOuterloop))
        W  <- svd(data)$v
        W <- W[, 1:q]     
        W <- W + matrix(rnorm(dim(W)[1]*dim(W)[2], 0, 0.0001 * (z - 1)), dim(W)[1], dim(W)[2])
        W[fixW == 0] <- 0

        for(j in 2:(maxItrOuterloop + 1)){ 
            #update P
            P <- updateP(XTX, W)

            #update Lossfunction, break if smaller than epsilon
            lossFunctionValue[j] <- updateLossfunction(data, N, W, P, RIDGE, LASSO)
            #print(lossFunctionValue[j])
            if(lossFunctionValue[j - 1] - lossFunctionValue[j]  < 10^-10){ break;}

            # coordinate descent           
            W1 <- coordinateDescent(data, P, W, q, k, n = N, RIDGE, LASSO, WINIT, XTXdiag)
            W <- W1$W
            if(j %% 1000 == 0){print(j)} 
        }
        loss <- lossFunctionValue[lossFunctionValue != Inf]
        minLoss <- loss[length(loss)]

        minLossVector[z] <- minLoss
        resultList[[z]] <- list(loss = loss, minLoss = minLoss, W = W, P = P)
    }
    return(resultList[[which(minLossVector == min(minLossVector))]])
}


determineLasso <- function(it = 200, XTX, W, LASSO, k, q, N,  percentage, data, RIDGE, WINIT, XTXdiag){
    for(i in 1:it){
        P <- updateP(XTX, W)
        LASSO <- percentageThresholding(W, k, q, percentage)
        #print(LASSO)
        W1 <- coordinateDescent(data, P, W, q, k, n = N, RIDGE, LASSO, WINIT, XTXdiag)
        W <- W1$W
    }
    return(LASSO)
}


updateP <- function(XTX, W){
  svd <- svd(XTX %*% W)   
  P <- svd$u %*% t(svd$v)  
  return(P)
}

updateLossfunction <- function(data, N, W, P, RIDGE, LASSO){
  out <- (1/(N*2))*sum((data - data %*% W %*% t(P))^2) + 
    (1/2)*RIDGE*sum(W^2) + LASSO %*% apply(abs(W), 2, sum)
  return(out)
}

percentageThresholding <- function(W, k, q, percentage){
  W <- abs(W)
  LASSO <- rep(0, q)
  for (i in 1:q){
    index <- ceiling(k * (1-percentage[i]))
    LASSO[i] <- W[ order(W[, i], decreasing = TRUE)[index], i ]
    #print(index)
    #print(LASSO)
  }
  return(LASSO)
}

coordinateDescent <- function(X, P, W, q, J, n, RIDGE, LASSO, WINIT, XTXdiag){

  for(d in 1:q){
    sumPRES <- X %*% ( P[, d] - W[, d] )
    
    for (j in 1:J){
      if(WINIT[j, d] != 0){
        wold <- W[j, d]
        CP <- (sum(sumPRES * X[, j]) + (wold * XTXdiag[j])) / n
        wols <- sign(CP) * (abs(CP) - LASSO[d])
        w <- (wols) / (RIDGE + (XTXdiag[j] /n ))
        if(abs(CP) < LASSO[d]){
          w <- 0
          sumPRES <- sumPRES + (X[, j] * wold)
        } else {
          sumPRES <- sumPRES + (X[, j] * (wold - w))
        }
        W[j, d] <- w
      }
    }
  }
  return(list(W = W))
}

################################################################################
# auxillary function to center and scale the data
# set value to 1 to scale with the unbiased version

scaleData <- function(X, value = 0){
  
  X <- scale(X, scale = FALSE)   
  attr(X, "scaled:center") <- NULL 
  sdX <-  apply(X, 2, function(x) sqrt( sum( x^2 ) / (length(x) - value )   ))  #compute the sd for each column  
  sdX <- matrix(sdX, nrow(X), ncol(X), byrow = T)                     #put all the sd's in a matrix
  
  X <- X * (1 / sdX)      #divide each entry in X by its sd
  return(X)
}

##################################################################

# scrips for the simulation study 

generateData <- function(x, nx, nfactors){
  X <- matrix(rnorm(x*nx, 0, 3), nx, x)
  X <- scaleData(X)
  XTX <- t(X) %*% X
  V <- svd(X)$v
  W <- V[, 1:nfactors]
  W[ abs(W) < 0.1 ] <- 0

  E <- matrix(rnorm(x*nx, 0, 0.1), nx, x)
  
  svd <- svd(XTX %*% W)   
  P <- svd$u %*% t(svd$v)  
  
  X <- X %*% W %*% t(P) + E
  return(list(X = X, W = W, P = P))
}

tuckerCongruence <- function(A, B){                                                                 
    combinationList <- combinat::permn(1:(ncol(A)))                                                 
    res <- -Inf                                                                                     
    fixA <- A                                                                                       
    changeInSign <- c(1, -1)
    placeHolder <- rep(NA, 2)

        for(combination in combinationList){                                                        
            A <- fixA[, combination]                                                                
            tucCon <- rep(NA, ncol(A))                                                              
               for( i in 1:ncol(A) ){                                                                  
                    for(j in c(1, 2)){
                    placeHolder[j] <- changeInSign[j]*A[, i] %*% B[, i]  / sqrt(sum(A[, i]^2) * sum(B[, i]^2)) 
                }
                tucCon[i] <- max(placeHolder)
            }                                                                                        
            candidate <- mean(tucCon)                                                               
            if(candidate > res){                                                                    
                res <- candidate                                                                    
            }                                                                                       
        }                                                                                           
    return(res)                                                                                     
}         



countZero <- function(A){
  return(sum( A == 0 ) / length(A))
} 

generateCommonSpecific <- function(x, nx, nfactors, errorFactor, 
                                   coefFixed = TRUE, sparsity){
  X <- matrix(rnorm(x*nx, 0, 3), nx, x)
  X <- scaleData(X)
  XTX <- t(X) %*% X
  V <- svd(X)$v
  W <- V[, 1:nfactors]
  W[ abs(W) < sparsity] <- 0
  W[1:floor(nrow(W)/2), 2] <- 0
  W[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0
  
  if(coefFixed){
    fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
    fixW[(floor(nrow(W)/2)+1):nrow(W), 3] <- 0
    fixW[1:(floor(nrow(W)/2)), 2] <- 0
  } else {
    fixW <- matrix(rnorm(dim(X)[2]*nfactors, 100), dim(X)[2], nfactors) 
  }
  
  E <- matrix(rnorm(x*nx, 0, 1), nx, x)
  
  svd <- svd(XTX %*% W)   
  P <- svd$u %*% t(svd$v)  
  Xtrue <- X %*% W %*% t(P)
  Xtrue <- Xtrue * errorFactor
  X <- Xtrue + E
  
  SStrue <- var(as.vector(Xtrue))
  SSX <- var(as.vector(X))
  
  return(list(X = X, W = W, P = P, errorRatio = 1 - SStrue/SSX,
              percentageZeroes = apply(W, 2, countZero), fixW = fixW))
}


correctlyClassified <- function(A, B){
  counter <- 0
  for( i in 1:nrow(A) ){
    for( j in 1:ncol(A) ){
      if(A[i, j] == B[i, j]){
        counter <- counter + 1
      }
      if(A[i, j] != 0 && B[i, j] != 0){
        counter <- counter + 1
      }
    }
  }
  return(counter / (nrow(A)*ncol(A)))
}



simRes <- function(error, fixednotFixed, sparsity){
  dat <- generateCommonSpecific(60, 60, 3, error, coefFixed = fixednotFixed, sparsity = sparsity)
  res <- sparsePCAMultistart(dat$X, nfactors = 3, RIDGE = 0.1, LASSO = c(0, 0, 0), 
                             percentage = dat$percentageZeroes, dat$fixW,
                             maxItrOuterloop = 1000000, percentageMode = TRUE)
  correctClass <- correctlyClassified(dat$W, res$W)
  tuckerCon <- tuckerCongruence(res$W, dat$W)
  return(list(dat = dat, res = res, correctClass =  correctClass, tuckerCon = tuckerCon))
}




