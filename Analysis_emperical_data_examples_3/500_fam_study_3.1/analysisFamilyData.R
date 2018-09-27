####################
# - analysis of the family data

#load required functions
rm(list=ls())

source("./SPARSE_PCA_wRandomStart.R")
source("./CVfunction.R")
library(sparseSCAcppFunction)

#########################
# Family data analysis


#to obtain the data see README
load('./family_data.RData')

str(family_data)

child <- family_data$child
mom <- family_data$mom
dad <- family_data$dad

newNamesChild <- paste("C:", names(child)) 
names(child) <- newNamesChild

namesDadMomChild <- c(names(dad), names(mom), names(child))


ncol(child)
ncol(dad)
ncol(mom)

rownames(child) <- 1:nrow(child)
rownames(mom) <- 1:nrow(mom)
rownames(dad) <- 1:nrow(dad)

dadMomChild <- cbind(dad, mom, child)
dadMomChild <- scale(as.matrix(dadMomChild))

#scree plot
eigenV <- eigen(cov(scale(as.matrix(dadMomChild))))$values 
plot(eigenV / sum(eigenV))

MSEcomponents <- matrix(NA, 10, 2) 
for(i in 1:nrow(MSEcomponents)){
    out  <- EigenVectorCV2(inputDat = dadMomChild, ridge = 0, fixW = matrix(1, ncol(dadMomChild), i), 
                           lasso = 0, nrFolds = 10, nfactors = i, maxItr = 10000)
    MSEcomponents[i, 1] <- out$MSE 
    MSEcomponents[i, 2] <- out$stdError
}



MSEcomponents <- as.data.frame(MSEcomponents)
names(MSEcomponents) <- c("MSE", "StandardError")

library(ggplot2)

setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(MSEcomponents, aes(x=as.factor(1:10), y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_hline(yintercept = MSEcomponents$MSE[which.min(MSEcomponents$MSE)] + MSEcomponents$StandardError[which.min(MSEcomponents$MSE)], linetype = "dashed")+ 
    geom_point()+ 
    xlab("Components")+
    ylab("MPRESS")+
    theme_bw() +
    ggsave('factorPlot_family.pdf', dpi=300)

system(paste("pdfcrop", 'factorPlot_family.pdf', 'factorPlot_family.pdf'))

fac <- 6


###############################################################################
#get all common and distinctive structures
allWs3comp <- allCommonSpecific(vars = c(8, 8, 7), components = fac)

length(allWs3comp)

#10 fold on all common distinctive structures
res3comp <- matrix(NA, length(allWs3comp), 2)
for(i in 1:length(allWs3comp)){
    out  <- EigenVectorCV2(inputDat = dadMomChild, ridge = 0, fixW = allWs3comp[[i]], 
                           lasso = 0, nrFolds = 10, nfactors = fac, maxItr = 10000)
    res3comp[i, 1] <- out$MSE 
    res3comp[i, 2] <- out$stdError
}

bestModel3comp <- which.min(res3comp[, 1])
allWs3comp[[bestModel3comp]]


##############################################################################
#tune lasso best model
lasso3 <- seq(2, 0, by = -0.01)
lasso3

lasso3comp <- matrix(NA, length(lasso3), 2)

for(i in 1:length(lasso3)){
    out  <- EigenVectorCV2(inputDat = dadMomChild, ridge = 0,
                           lasso = lasso3[i], nrFolds = 10,
                           nfactors = fac, fixW = allWs3comp[[bestModel3comp]],
                           maxItr = 10000)
    lasso3comp[i, 1] <- out$MSE
    lasso3comp[i, 2] <- out$stdError
}

plot(lasso3comp[, 1])
minLas <- which.min(lasso3comp[, 1])
bestLasso <- lasso3[minLas]
bestLasso

analysisLasso <- lasso3[lasso3comp[, 1] <  lasso3comp[minLas, 1] + lasso3comp[minLas, 2]][1]

bestLasso <- analysisLasso

#make a nice plot for the lasso tuning
library(ggplot2)

frameForPlot <- as.data.frame(cbind(lasso3, lasso3comp))
frameForPlot

colnames(frameForPlot) <- c("lasso", "MSE", "StandardError")
frameForPlot



setwd('~/Desktop/PhD/sparsePCAPaper/abstract/images')

ggplot(frameForPlot, aes(x=lasso, y=MSE)) + 
    geom_errorbar(aes(ymin = MSE - StandardError, ymax = MSE + StandardError)) + 
    geom_point()+ 
    guides(shape = guide_legend(title="% of zeroes"))+
    geom_hline(yintercept= frameForPlot$MSE[which.min(frameForPlot$MSE)]
               + frameForPlot$StandardError[which.min(frameForPlot$MSE)] ,
               linetype="dashed", color = "black") +
    geom_point(data = frameForPlot[ frameForPlot$lasso == analysisLasso, ], size = 3, shape = 3) +
    geom_vline(xintercept = frameForPlot$lasso[which(frameForPlot$lasso == analysisLasso)])+
    xlab("Lasso") + 
    ylab("MPRESS") +
    theme_bw() +
    ggsave('lassoPlot_familyData.pdf', dpi=300)

system(paste("pdfcrop", 'lassoPlot_familyData.pdf', 'lassoPlot_familyData.pdf'))

######################################################################################

res <- sparsePCAMultistart(dadMomChild, nfactors = 3, RIDGE = 0, LASSO = bestLasso ,
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

rownames(res$W) <- namesDadMomChild
res$W

###################################################################################


#run the models for in Table 1
res1 <- sparsePCAMultistart(dadMomChild, nfactors = fac, RIDGE = 0, LASSO = 0,
                            percentage, fixW = matrix(1, ncol(dadMomChild), 6),
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

res2 <- sparsePCAMultistart(dadMomChild, nfactors = fac, RIDGE = 0, LASSO = 0,
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

res3 <- sparsePCAMultistart(dadMomChild, nfactors = fac, RIDGE = 0, LASSO = bestLasso,
                            percentage, fixW = allWs3comp[[bestModel3comp]],
                            maxItrOuterloop = 100000,
                            percentageMode = FALSE)

rownames(res3$W) <- namesDadMomChild
res3$W

allWs3comp[[bestModel3comp]]

#function to get the variance accounted for
VAFw <- function(inputDat, analysisResults){
    #get sum of squares of the data
    SS <- sum(scale(inputDat, scale = FALSE)^2)
    SSmodel <- sum((inputDat - inputDat %*% analysisResults$W %*% t(analysisResults$P))^2)
    return(1 - SSmodel/SS)
}


#VAFS for the three models
round(VAFw(dadMomChild, res1)*100, 1)
round(VAFw(dadMomChild, res2)*100, 1)
round(VAFw(dadMomChild, res3)*100, 1)

library(xtable)

#make the SCaDS table
table1 <- res3$W
table1 <- round(table1, 3)
rownames(table1) <- namesDadMomChild
colnames(table1) <- paste("T" , 1:6, sep = "")

#make the varimax table
table2 <- matrix(as.vector(varimax(res1$W)$loadings), 23, 6)
#table2 <- varimax(res1$W)$loadings
table2 <- round(table2, 3)
rownames(table2) <- namesDadMomChild
colnames(table2) <- paste("T" , 1:6, sep = "")

#make the DISCO table
table3 <- Wrot
table3 <- round(table3, 3)
rownames(table3) <- namesDadMomChild
colnames(table3) <- paste("T" , 1:6, sep = "")
table3


besteStructuur <-  allWs3comp[[bestModel3comp]]
rownames(besteStructuur) <- namesDadMomChild 
colnames(besteStructuur) <- 1:6


#save.image(file = './resPaper.RData')
load('./resPaper.RData')

##############################################################################
# script to get the rotateed W matrix

library(devtools)
install_github("ZhengguoGu/RegularizedSCA")
library(RegularizedSCA)

rotationmatrix <- pstr(res1$W,besteStructuur,1-besteStructuur,500,0.0001)
Wrot <- res1$W%*%rotationmatrix$Bmatrix
Prot <- res1$W%*%rotationmatrix$Bmatrix
#number of zeroes in non-zero constrained part
nrzero <- colSums(res3$W==0)#-colSums(besteStructuur==0)
#set smallest elements in rotated matrix to zero
for (q in 1:fac) {
  i <- order(abs(Wrot[,q]),decreasing = FALSE)
  Wrot[i[1:nrzero[q]],q] <- 0
}
cbind(Wrot,res3$W)#zelfde comm.dist.struct maar andere selectie

#calculate %vaf
Xhat <- dadMomChild%*%Wrot%*%t(Prot)
DEVsq <- sum(colSums((Xhat-dadMomChild)^2))
ssqX <- sum(colSums(dadMomChild^2))
vaf <- 1-DEVsq/ssqX
round(vaf*100,1)

#res4<-list(loss = res3$loss, minLoss = res3$minLoss, W = Wrot, P = Prot)
#round(VAFw(dadMomChild, res4)*100, 1)

#########################
# function pstr
pstr <- function (P, Target, W, maxiter, convergence) 
{
  n <- dim(P)[1]
  m <- dim(P)[2]
  L <- array()
  Bmat <- list()
  REFL <- reflexmat(m)
  for (i in 1:dim(REFL)[1]) {
    k <- which(REFL[i, ] == -1)
    Binit <- diag(m)
    Binit[, k] <- -1 * Binit[, k]
    B1 <- t(P) %*% P
    alpha <- max(eigen(B1)$values)
    iter <- 1
    stop <- 0
    Bcurrent <- Binit
    Lossc <- pstrLoss(Binit, P, Target, W)
    while (stop == 0) {
      Pw <- W * Target + P %*% Bcurrent - W * (P %*% Bcurrent)
      A <- -2 * t(Pw) %*% P
      Fmat <- A + 2 * t(Bcurrent) %*% t(B1) - 2 * alpha * 
        t(Bcurrent)
      F_svd <- svd(-Fmat)
      B <- F_svd$v %*% t(F_svd$u)
      if (iter == maxiter) {
        stop <- 1
      }
      Loss <- pstrLoss(B, P, Target, W)
      Diff <- Lossc - Loss
      if (abs(Diff) < convergence) {
        stop <- 1
      }
      iter <- iter + 1
      Lossc <- Loss
      Bcurrent <- B
    }
    L[i] <- Lossc
    Bmat[[i]] <- Bcurrent
  }
  k <- which(L == min(L))
  Loss <- L[k[1]]
  B <- Bmat[[k[1]]]
  results <- list()
  results$Bmatrix <- B
  results$Loss <- Loss
  return(results)
}
#function reflexmat
reflexmat <- function (m) 
{
  mat <- rep(1, m)
  for (i in 1:(m - 1)) {
    B <- utils::combn(1:m, i)
    for (j in 1:dim(B)[2]) {
      v <- rep(1, m)
      v[t(B[, j])] <- -1
      mat <- rbind(mat, v)
    }
  }
  return(mat)
}

#function prstrloss
pstrLoss <- function (B, Tmat, Target, W) 
{
  DEV <- Tmat %*% B - Target
  wDEV <- W * DEV
  Loss <- sum(wDEV^2)
  return(Loss)
}


