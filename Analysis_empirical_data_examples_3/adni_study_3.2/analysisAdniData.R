rm(list=ls())

source("./SPARSE_PCA_wRandomStart.R")
source("./CVfunction.R")
library(sparseSCAcppFunction)

#################################################################################################
# pre prepping of the data

load("./ADNI/DataUsedforPaper/merge data/ADNI_final.RData")
ls()

dim(genes_data)
dim(neuropsy_data)

class(genes_data)
class(neuropsy_data)

colnames(genes_data)
colnames(neuropsy_data)


#scale the data so the SS in each block is the same
#Divide by sqrt of Q
genes_data <- scaleData(genes_data) * ( 1 / sqrt(ncol(genes_data)))
neuropsy_data <- scaleData(neuropsy_data) * ( 1 / sqrt(ncol(neuropsy_data)))


#bind the data again together
dat <- cbind(genes_data, neuropsy_data)
tail(dat)



varnames <- c(paste("1:", colnames(genes_data)), paste("2:", colnames(neuropsy_data)))
varnames

 
###################################################################################################
#MSEf <- rep(NA, 10)


#choose factors to estimate
MSEf <- rep(NA, 10)
MSEStdErrf <- rep(NA, 10)


for(i in 1:10){
    res <- EigenVectorCV2(inputDat = dat, ridge = 0, lasso = 0,
                   fixW = matrix(1, ncol(dat), i) ,nrFolds = 10, nfactors = i, maxItr = 10000 ) 
    MSEf[i] <- res$MSE
    MSEStdErrf[i] <- res$stdError
}

#a 4 factor solution


# component plot for the paper
MSEcomponents <- as.data.frame(cbind(MSEf, MSEStdErrf))
names(MSEcomponents) <- c("MSE", "StandardError")

library(ggplot2)

#setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(MSEcomponents, aes(x=as.factor(1:10), y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_hline(yintercept = MSEcomponents$MSE[which.min(MSEcomponents$MSE)] + MSEcomponents$StandardError[which.min(MSEcomponents$MSE)], linetype = "dashed")+ 
    geom_point()+ 
    xlab("Factors")+
    ylab("MPRESS")+
    theme_bw() 
    #ggsave('factorPlot_adni.pdf', dpi=300)

#system(paste("pdfcrop", 'factorPlot_adni.pdf', 'factorPlot_adni.pdf'))


####################################################################################################
# Tuning of the ridge parameter
seqRidge <- seq(2, 0, by = -0.01)
MSE <- rep(NA, length(seqRidge))
MSEStdErr <- rep(NA, length(seqRidge))

for(i in 1:length(seqRidge)){
    res <- EigenVectorCV2(inputDat = dat, ridge = seqRidge[i], lasso = 0,
                   fixW = matrix(1, ncol(dat), 4) ,nrFolds = 10, nfactors = 4, maxItr = 10000 ) 
    MSE[i] <- res$MSE
    MSEStdErr[i] <- res$stdError
}

plot(MSE)
plot(MSEStdErr)
bestRidge <- MSE[which.min(MSE)]

withinOneStdErr <- MSE < MSE[which.min(MSE)] + MSEStdErr[which.min(MSE)]
highestRidgeWithinOneStdError <- seqRidge[withinOneStdErr][1]
highestRidgeWithinOneStdError

#check the plot for the ridge parameter
x <- 1:length(seqRidge)
plot(x, MSE,
    ylim=range(c(MSE - MSEStdErr, MSE + MSEStdErr)),
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)

arrows(x, MSE - MSEStdErr, x, MSE + MSEStdErr, length=0.05, angle=90, code=3)



#make ridge plot for paper
frameForPlotRidge <- as.data.frame(cbind(seqRidge, MSE, MSEStdErr))
frameForPlotRidge

colnames(frameForPlotRidge) <- c("ridge", "MSE", "StandardError")
frameForPlotRidge
#setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(frameForPlotRidge, aes(x=ridge, y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_point()+ 
    geom_hline(yintercept= frameForPlotRidge$MSE[which.min(frameForPlotRidge$MSE)]
               + frameForPlotRidge$StandardError[which.min(frameForPlotRidge$MSE)] ,
               linetype="dashed", color = "black") +
    geom_point(data = frameForPlotRidge[ frameForPlotRidge$ridge == highestRidgeWithinOneStdError, ]
               ,size = 3, shape = 3) +
    geom_vline(xintercept = frameForPlotRidge$ridge[which( frameForPlotRidge$ridge == highestRidgeWithinOneStdError)])+
    xlab("Ridge") + 
    ylab("MPRESS") +
    theme_bw() +
    #ggsave('ridgePlot_adni.pdf', dpi=300)

#system(paste("pdfcrop", 'ridgePlot_adni.pdf', 'ridgePlot_adni.pdf'))

###################################################################################################
# Tuning of the common distinctive structure
allCombn <- allCommonSpecific(c(388, 12), 4)
MSEcombn <- rep(NA, length(allCombn))
MSEStdErrCombn <- rep(NA, length(allCombn))

for(i in 1:length(allCombn)){
    res <- EigenVectorCV2(inputDat = dat, ridge = highestRidgeWithinOneStdError, lasso = 0,
                   fixW = allCombn[[i]], nrFolds = 10, nfactors = 4, maxItr = 10000 ) 
    MSEcombn[i] <- res$MSE
    MSEStdErrCombn[i] <- res$stdError
}

bestModel <- which.min(MSEcombn)
#the best model is model 4, which is a model with 4 distinctive components two for each block
allCombn[which.min(MSEcombn)]


withinOneStdErr <- MSEcombn < MSEcombn[which.min(MSEcombn)] + MSEStdErrCombn[which.min(MSEcombn)]
allCombn[withinOneStdErr]




x <- 1:length(allCombn)
plot(x, MSEcombn,
    ylim=range(c(MSEcombn - MSEStdErrCombn, MSEcombn + MSEStdErrCombn)),
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)

arrows(x, MSEcombn - MSEStdErrCombn, x, MSEcombn + MSEStdErrCombn, length=0.05, angle=90, code=3)


# make a plot for the paper
# first create data frame for ggplot2
frameForPlot <- as.data.frame(cbind(1:15, MSEcombn, MSEStdErrCombn, unlist(lapply(allCombn, countZero)) ))
colnames(frameForPlot) <- c("model", "MSE", "StandardError", "percentageOfZeroes")
frameForPlot$percentageOfZeroes <- as.factor(frameForPlot$percentageOfZeroes) 
frameForPlot$model <- as.factor(frameForPlot$model)

#create nice X labels for plot
coding <- as.character(unlist(lapply(allCombn,  colSums )))
coding[coding == 388]  <- "D1"
coding[coding == 12] <- "D2"
coding[coding == 400] <- "C"

commonDistinctive <- rep(NA, length(coding) / 4)

for(i in 1:(length(coding) /4)){
    index <- which(rep(1:15, each = 4) ==  i)
    commonDistinctive[i] <- paste(coding[index], collapse = " ")
}

frameForPlot$commonDistinctive <- commonDistinctive



library(ggplot2)

setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(frameForPlot, aes(x=commonDistinctive, y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    ylab("MPRESS") +
    xlab("Model") +
    geom_point() + 
    geom_point(data = frameForPlot[which.min(frameForPlot$MSE), ], shape=8, size = 4) +
    geom_hline(yintercept = frameForPlot$MSE[which.min(frameForPlot$MSE)] +
               frameForPlot$StandardError[which.min(frameForPlot$MSE)], linetype="dashed") +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5))  
    ggsave('modelsADNI.pdf', dpi=300)

system(paste("pdfcrop", 'modelsADNI.pdf', 'modelsADNI.pdf'))


#model 5 seems appropriate
lapply(allCombn, countZero)
order(MSEcombn)

MSEcombn

###################################################################################################
# Tuning of the lasso parameter
lassoSeq <- seq(0.1, 0, by = -0.0001)
MSElasso <- rep(NA, length(lassoSeq))
MSEStdErrLasso <- rep(NA, length(lassoSeq))

for(i in 1:length(lassoSeq)){
    res <- EigenVectorCV2(inputDat = dat, ridge = highestRidgeWithinOneStdError,
                          lasso = lassoSeq[i], fixW = allCombn[[bestModel]], nrFolds = 10,
                          nfactors = 4, maxItr = 10000 ) 
    MSElasso[i] <- res$MSE
    MSEStdErrLasso[i] <- res$stdError
}


x <- 1:length(lassoSeq)
plot(x, MSElasso,
    ylim=range(c(MSElasso - MSEStdErrLasso, MSElasso + MSEStdErrLasso)),
    pch=19, xlab="Measurements", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars"
)

arrows(x, MSElasso - MSEStdErrLasso, x, MSElasso + MSEStdErrLasso, length = 0.05, angle = 90, code=3)

alarm()

withinOneStdErr <- MSElasso < MSElasso[which.min(MSElasso)] + MSEStdErrLasso[which.min(MSElasso)]
highestLassoWithinOneStdError <- lassoSeq[withinOneStdErr][1]
highestLassoWithinOneStdError



#make lasso plot for paper
frameForPlotLasso <- as.data.frame(cbind(lassoSeq, MSElasso, MSEStdErrLasso))


#make a subset of the dots because there are too many dots
frameForPlotLasso <- frameForPlotLasso[rep(c(TRUE, FALSE), 1001 )[1:nrow(frameForPlotLasso) ], ]
 

library(ggplot2)
colnames(frameForPlotLasso) <- c("lasso", "MSE", "StandardError")
frameForPlotLasso
setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(frameForPlotLasso, aes(x=lasso, y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_point()+ 
    geom_hline(yintercept= frameForPlotLasso$MSE[which.min(frameForPlotLasso$MSE)]
               + frameForPlotLasso$StandardError[which.min(frameForPlotLasso$MSE)] ,
               linetype="dashed", color = "black") +
    geom_point(data = frameForPlotLasso[ frameForPlotLasso$lasso == highestLassoWithinOneStdError, ]
               ,size = 3, shape = 3) +
    geom_vline(xintercept = frameForPlotLasso$lasso[which( frameForPlotLasso$lasso == highestLassoWithinOneStdError)])+
    xlab("Lasso") + 
    ylab("MPRESS") +
    theme_bw() 
    ggsave('lassoPlot_adni.pdf', dpi=300)

system(paste("pdfcrop", 'ridgePlot_adni.pdf', 'ridgePlot_adni.pdf'))


###################################################################################################
# analysis of the data with the chosen parameters

res <- sparseSCAcpp(dat, Q = 4, RIDGE = highestRidgeWithinOneStdError,
                    LASSO = rep(highestLassoWithinOneStdError , 4),
                    fixW = allCombn[[bestModel]], maxItrOuterloop = 10000, nStarts = 1,
                    print = TRUE, tol = 10^-12)

rownames(res$W) <- varnames
res$W <- scale(res$W, center = FALSE)  
res$W
 


#make heatplot table for paper
#make data frame suitable for ggplot
library(reshape2)
frameForPlot2 <- melt(abs(res$W))

group <- rep(c(rep(1, 100), rep(2, 100), rep(3, 100), rep(4, 100)), 4)
frameForPlot2 <- cbind(frameForPlot2, group) 

colnames(frameForPlot2) <- c("Variable", "Component", "value", "group")

setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(frameForPlot2, aes(Component, Variable)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(low = "white", high = "black") +
    facet_wrap(~ group, ncol = 4, scale = "free") +
    ylab("Absolute component weight") +
    theme(strip.background = element_blank(),
      strip.text.x = element_blank()) +
    ggsave('heatplotADNI.pdf', dpi=300)

system(paste("pdfcrop", 'heatplotADNI.pdf', 'heatplotADNI.pdf'))



################################################################################

