################################################
# cross validation of ridge for the conditions
################################################

rm(list=ls())
################################################################################################
require(Rcpp)
source('./SPARSE_PCA_wRandomStart.R')
source('./crossValidation.R')
source('./generateData.R')
sourceCpp('./sparseSCA.cpp')



################################################################################################

#create conditions
pvec <- c(0.05, 0.25, 0.5)
sparvec <- c(TRUE, FALSE)
conditions <- expand.grid(pvec, sparvec)
conditions


ridgeTuning <- list() 

set.seed(123)

ridgeVec <- seq(1, 0, by = -0.01)



for(j in 1:6){
        if(conditions[j, 2]){
            sparsity <- c(.6, .6, .6)
        } else {

            sparsity <- c(.02, .52, .52)
        }

    dat <- generateCommonSpecific(500, 100, 3, p = conditions[j, 1],
                              coefFixed = TRUE, sparsity = sparsity)
    dat$X <- scale(dat$X, scale = FALSE) 

    res <- matrix(NA, length(ridgeVec), 2)

    for(i in 1:length(ridgeVec)){

        error <- EigenVectorCV2(dat$X, ridge = ridgeVec[i], lasso = 0,
                                nrFolds = 10, nfactors = 3, fixW = dat$fixW,
                                maxItr = 100000, scaleDat = FALSE)
        res[i, 1] <- error$MSE
        res[i, 2] <- error$stdError

    }
    ridgeTuning[[j]] <- res
    print(j)

}

plot(res[, 1])

ridgeTuning


res <- ridgeTuning[[1]]


x <- 1:length(res[, 1])
plot(x, res[, 1],
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)
arrows(x, res[,1]-res[,2], x, res[,1]+res[,2], length=0.05, angle=90, code=3)


ridge <- rep(NA, 6)
for(i in 1:length(ridgeTuning)){

    res <- ridgeTuning[[i]]

    withinOneStdErr <- res[, 1] < (res[which.min(res[, 1]), 1] + res[which.min(res[, 1]), 2] )
    ridge[i]  <- ridgeVec[withinOneStdErr][1] 
}


saveRDS(ridge, './ridgeVec')

##################################################################################################


