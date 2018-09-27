############################################################################################
# produce the boxplots and the table in the low dimensional condition
############################################################################################
rm(list=ls())

modelRes <- readRDS('./modelselectionResults_Good.rds')


###############
source('./CVfunction.R')
source('./SPARSE_PCA_wRandomStart.R')

models <- allCommonSpecific(c(250, 250), 3)
zeroes <- lapply(models, countZero)

percentageOfZeroes <- unlist(zeroes) 
percentageOfZeroes

conditions <- length(modelRes)
reps <- length(modelRes[[1]])


#######################################################################
# create the boxplots
dat1 <- cbind(modelRes[[1]][[1]]$res, modelRes[[1]][[1]]$stdErrorRes)
dat2 <- cbind(modelRes[[2]][[2]]$res, modelRes[[2]][[2]]$stdErrorRes)
dat3 <- cbind(modelRes[[3]][[3]]$res, modelRes[[3]][[3]]$stdErrorRes)


cond <- c(rep("high sparsity, 5% error", 10),
rep("high sparsity, 25% error", 10),
rep("high sparsity, 50% error", 10))


frameForPlot <- as.data.frame(cbind(rep(1:10, 3), rbind(dat1, dat2, dat3), cond))
colnames(frameForPlot) <- c("model", "MSE", "StandardError", "condition")

frameForPlot$MSE <- as.numeric(levels(frameForPlot$MSE)[frameForPlot$MSE]) 
frameForPlot$StandardError <- as.numeric(levels(frameForPlot$StandardError)[frameForPlot$StandardError]) 
frameForPlot$model <- factor(frameForPlot$model, levels = as.character(1:10))
frameForPlot$condition <- factor(frameForPlot$condition, c("high sparsity, 5% error", "high sparsity, 25% error", "high sparsity, 50% error"))

frameForPlot$percentageOfZeroes <- rep(round(percentageOfZeroes, 2), 3)
levelsPercentageZeroes <- as.character(round(unique(percentageOfZeroes), 2)) 
frameForPlot$percentageOfZeroes <- as.factor(frameForPlot$percentageOfZeroes)

frameForPlot

#create nice xlabels for plot 
allCombn <- allCommonSpecific(c(388,12), 3)

coding <- as.character(unlist(lapply(allCombn,  colSums )))
coding[coding == 388]  <- "D1"
coding[coding == 12] <- "D2"
coding[coding == 400] <- "C"
coding

commonDistinctive <- rep(NA, length(coding) / 3)

for(i in 1:(length(coding) /3)){
    index <- which(rep(1:15, each = 3) ==  i)
    commonDistinctive[i] <- paste(coding[index], collapse = " ")
}
frameForPlot$commonDistinctive <- commonDistinctive


library(ggplot2)


setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

# Standard error of the mean
ggplot(frameForPlot, aes(x=commonDistinctive, y=MSE, shape = percentageOfZeroes)) + 
    #geom_errorbar(aes(ymin=MSE-Standard error, ymax=MSE+Standard error), width=.1) +
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_point()+ 
    ylab("MPRESS")+
    xlab("Model") +
    facet_wrap(. ~ condition, scales="free_y")+
    guides(shape = guide_legend(title="% of zeroes"))+
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5))  

    ggsave('allComDisPlot_Good.pdf', dpi=300)


system(paste("pdfcrop", 'allComDisPlot_Good.pdf', 'allComDisPlot_Good.pdf'))

###########################################################################################
# Make table

library(stringr)

#count the amount of distinctive components
distinctiveComponentsInStructure <- str_count(commonDistinctive, "D")


bestModel <- rep(NA, conditions)
withinOneStdError <- rep(NA, conditions)
didTheBestModelComeFirst <- rep(NA, conditions)
rule <- rep(NA, conditions) 
for(i in 1:conditions){
    counter <- 0
    counter2 <- 0
    counter3 <- 0
    counter4 <- 0
    for(j in 1:reps){
        if(modelRes[[i]][[j]]$orderOfModels[1] == 5){ 
            counter3 <- counter3 + 1
        }
        if(5 %in% modelRes[[i]][[j]]$orderOfModels[1:4]){
            counter <- counter + 1
        }
        if(modelRes[[i]][[j]]$withinOneStdError[5] == TRUE){
            counter2 <- counter2 + 1
        }
        if(modelRes[[i]][[j]]$withinOneStdError[5] == TRUE &&
            max(distinctiveComponentsInStructure[modelRes[[i]][[j]]$withinOneStdError]) == 2){
            counter4 <- counter4 + 1
        }

    }
    rule[i]  <- counter4 /reps
    didTheBestModelComeFirst[i]  <- counter3 / reps
    bestModel[i] <- counter / reps
    withinOneStdError[i] <- counter2 / reps
}


rule
didTheBestModelComeFirst
bestModel
withinOneStdError


tableDf <- data.frame(sparsity = c(rep("high", 3), rep("low", 3)),
                      noise = rep(c("5%", "25%", "50%"), 2),
                      bestModel = didTheBestModelComeFirst,
                      withinOneStdError = rule
                      )

#table for paper
tableDf


require(xtable)
xtable(tableDf)



