#################################################################################################
# Function to tune lasso, so that the chosen lasso leads to
# The number of zeroes in the generated data
#################################################################################################

findLasso <- function(dat, maxItr, ridge){

    percentageZeroesInData <- mean(dat$percentageZeroes)
    percentageInW <- 0
    i  <- 0
    lassou <- 1 
    lassol <- 0 

    converged <- FALSE

    while( abs(percentageZeroesInData - percentageInW) > 0.001 && i <= maxItr ){

        lasso <- (lassol + lassou) / 2
        fit <- sparseSCAcpp(dat$X, Q = 3, RIDGE = ridge, LASSO = rep(lasso , 3), fixW = dat$fixW, 
                     maxItrOuterloop = 100000, nStarts = 1, print = FALSE, tol = 10^-8)
        percentageInW <- countZero(fit$W)
        if( percentageZeroesInData > percentageInW){
            lassol  <- lasso

        } else {
            lassou  <- lasso
        }
       print(lasso) 
       i <- i + 1 
    }

    if( i < maxItr ){
        converged <- TRUE
    } 

    return(list(lasso = lasso, converged = converged))
}


