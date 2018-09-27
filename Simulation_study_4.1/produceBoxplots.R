#############################################################################
# results simulation study 
###############################

simRes1 <- readRDS('./simres4.rds')
# see how the simulation results have been structured


simRes1[[1]][[1]]$tuckerConSparseSCA


#check what is object
#a list for the 12 conditions
class(simRes1)
length(simRes1)

#check what is in the object
class(simRes1[[1]])
length(simRes1[[1]])

#check what the iteration results look like
class(simRes1[[1]][[1]])
length(simRes1[[1]][[1]])
names(simRes1[[1]][[1]])

#check what simRes1[[1]][[1]]$dat lookslike
names(simRes1[[1]][[1]]$dat)

#check what simRes1[[1]][[1]]$res lookslike
names(simRes1[[1]][[1]]$resSparseSCA)
names(simRes1[[1]][[1]]$resElasticnet)

##########################################################################################333

tuckerConSparseSCAMat <- matrix(NA, 20, 6)
tuckerConElasticnetMat <- matrix(NA, 20, 6)

for(i in 1:length(simRes1)){
    for(j in 1:length(simRes1[[1]])){ 

        tuckerConSparseSCAMat[j, i] <- simRes1[[i]][[j]]$tuckerConSparseSCA
        tuckerConElasticnetMat[j, i] <- simRes1[[i]][[j]]$tuckerConElasticnet
    }
}

tuckerCongruenceValues <- c(as.vector(tuckerConSparseSCAMat[, 1:3]),
as.vector(tuckerConElasticnetMat[, 1:3]),
as.vector(tuckerConSparseSCAMat[, 4:6]),
as.vector(tuckerConElasticnetMat[, 4:6]))

length(tuckerCongruenceValues)

#boxplot for tucker congruence W

sparsityVector <- c(rep("High levels of sparsity", 120), rep("Low levels of sparsity", 120))
fixednotFixed <- rep(c(rep("SCaDS", 60), rep("SPCA", 60)), 2)
errorRatio <- rep(c(rep("5%", 20), rep("25%", 20), rep("50%", 20)), 4)
errorRatio <- factor(errorRatio,  levels = c("5%", "25%", "50%"), ordered = TRUE)

boxplotDf2 <- data.frame(tuckerCongruenceValues, sparsityVector, fixednotFixed, errorRatio = errorRatio)

library(ggplot2)


# for tuckercongruence W
# generate plot for the paper

setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

plot <- ggplot(boxplotDf2, aes(y = tuckerCongruenceValues, x = fixednotFixed))
plot + geom_boxplot(aes(fill = errorRatio)) +
    facet_grid(. ~ sparsityVector) +
    geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
    guides(fill=guide_legend(title="Error ratio")) +
    xlab("Estimation method") + 
    ylab("Tucker congruence") +
    coord_fixed(ratio = 26/5) +
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
    theme_bw() 

    ggsave('plotTuckerConGgplot_new.pdf', dpi=300)

system(paste("pdfcrop", 'plotTuckerConGgplot_new.pdf', 'plotTuckerConGgplot_new.pdf'))


################################################################################
# Correctly classified boxplot


names(simRes1[[1]][[1]])

classifiedSparseSCAMat <- matrix(NA, 20, 6)
classifiedElasticnetMat <- matrix(NA, 20, 6)

for(i in 1:length(simRes1)){
    for(j in 1:length(simRes1[[1]])){ 

        classifiedSparseSCAMat[j, i] <- simRes1[[i]][[j]]$correctClassSparseSCA
        classifiedElasticnetMat[j, i] <- simRes1[[i]][[j]]$correctClassElasticnet
    }
}

correctClassified <- c(as.vector(classifiedSparseSCAMat[, 1:3]),
as.vector(classifiedElasticnetMat[, 1:3]),
as.vector(classifiedSparseSCAMat[, 4:6]),
as.vector(classifiedElasticnetMat[, 4:6]))

correctClassified

#boxplot for tucker congruence W

sparsityVector <- c(rep("High levels of sparsity", 120), rep("Low levels of sparsity", 120))
fixednotFixed <- rep(c(rep("SCaDS", 60), rep("SPCA", 60)), 2)
errorRatio <- rep(c(rep("5%", 20), rep("25%", 20), rep("50%", 20)), 4)
errorRatio <- factor(errorRatio,  levels = c("5%", "25%", "50%"), ordered = TRUE)

boxplotDf2 <- data.frame(correctClassified, sparsityVector, fixednotFixed, errorRatio = errorRatio)
boxplotDf2

plot <- ggplot(boxplotDf2, aes(y = correctClassified, x = fixednotFixed))
plot + geom_boxplot(aes(fill = errorRatio)) +
    facet_grid(. ~ sparsityVector) +
    guides(fill=guide_legend(title="Error ratio")) +
    xlab("Estimation method") + 
    ylab("Correctly classified") +
    coord_fixed(ratio = 26/5) +
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
    theme_bw() +
    ggsave('correctClassified.pdf', dpi=300)


system(paste("pdfcrop", 'correctClassified.pdf', 'correctClassified.pdf'))




