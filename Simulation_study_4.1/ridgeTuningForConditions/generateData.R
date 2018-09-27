##################################################################################
#  Ffunction to generate common specific data with
##################################################################################
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


