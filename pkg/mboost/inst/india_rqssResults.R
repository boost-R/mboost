#-----------------------------------------------------------
# Explore results for empirical risk with all algorithms 
# including rqss
#-----------------------------------------------------------


#-----------------------------------------------------------
# Path, Packages, source functions
#-----------------------------------------------------------

library("quantreg")
library("lattice")
library("mboost")

load("india.Rdata")
load("plotresults.Rdata")

source("india_rqss_lambdaOptFunc.R")

#-----------------------------------------------------------
# Quantiles and Indices
#-----------------------------------------------------------
inds <- paste("w", 1:50, sep="")
quantiles <- c(0.05, 0.1, 0.5)

pred1 <- pred2 <- pred3 <- matrix(0, ncol=5, nrow=length(inds))
colnames(pred1) <- colnames(pred2) <- colnames(pred3) <- c("rqss", "add","vcm","stump","tree")

lambda1rqss <- lambda2rqss <- lambda3rqss <- matrix(NA, nrow=length(inds), ncol=6)

predAdditive1 <- predAdditive2 <- predAdditive3 <- matrix(NA, ncol=50, nrow=nrow(india)/3)
predVCM1 <- predVCM2 <-predVCM3 <- matrix(NA, ncol=50, nrow=nrow(india)/3)
nTest <- rep(0, length(inds))

# write cvRisks into a vector
for(i in 1:length(inds)){
  print(i)
  
  indiaTrain <- india[india[,inds[i]]==1,]
  
  indiaTestOrg <- india[india[,inds[i]]==3,] 
  
  # Add variable for prediction with mboost
  indiaTestOrg <- cbind(indiaTestOrg, ff = 0)
  names(indiaTestOrg)[ncol(indiaTestOrg)] <- "FALSE"
  
  # Indices of convex hull of the test data
  indiaTestInd <- convexHullInd(dataOrg=indiaTestOrg, dataRef=indiaTrain)
  indiaTest <- convexHull(dataOrg=indiaTestOrg, dataRef=indiaTrain)
  nTest[i] <- nrow(indiaTest)
  
  #----------------
  # Prediction rqss
  load(paste("rqssObjects/model",i,"_0.05.Rdata",sep=""))
  pred <- predict(model, newdata=indiaTest)
  pred1[i,1] <- loss(indiaTest$stunting, pred, tau=0.05)
  lambda1rqss[i,] <- model$lambdas
  
  load(paste("rqssObjects/model",i,"_0.1.Rdata",sep=""))
  pred <- predict(model, newdata=indiaTest)
  pred2[i,1] <- loss(indiaTest$stunting, pred, tau=0.1)
  lambda2rqss[i,] <- model$lambdas

  load(paste("rqssObjects/model",i,"_0.5.Rdata",sep=""))
  pred <- predict(model, newdata=indiaTest)
  pred3[i,1] <- loss(indiaTest$stunting, pred, tau=0.5)
  lambda3rqss[i,] <- model$lambdas
  
  #---------------
  # Prediction additive
  load(paste("additive/model",i,"_0.05.Rdata",sep=""))
  predAdditive1[,i] <- predict(model, newdata=indiaTestOrg)
  pred1[i,2] <- loss(indiaTest$stunting, predAdditive1[indiaTestInd==0,i], tau=0.05)
  
  load(paste("additive/model",i,"_0.1.Rdata",sep=""))
  predAdditive2[,i] <- predict(model, newdata=indiaTestOrg)
  pred2[i,2] <- loss(indiaTest$stunting, predAdditive2[indiaTestInd==0,i], tau=0.1)
   
  load(paste("additive/model",i,"_0.5.Rdata",sep=""))
  predAdditive3[,i] <- predict(model, newdata=indiaTestOrg)
  pred3[i,2] <- loss(indiaTest$stunting, predAdditive3[indiaTestInd==0,i], tau=0.5)
  
  #---------------
  # Prediction VCM
  load(paste("VCM/model",i,"_0.05.Rdata",sep=""))
  predVCM1[,i] <- predict(model, newdata=indiaTestOrg)
  pred1[i,3] <- loss(indiaTest$stunting, predVCM1[indiaTestInd==0,i], tau=0.05)
  
  load(paste("VCM/model",i,"_0.1.Rdata",sep=""))
  predVCM2[,i] <- predict(model, newdata=indiaTestOrg)
  pred2[i,3] <- loss(indiaTest$stunting, predVCM2[indiaTestInd==0,i], tau=0.1)
  
  load(paste("VCM/model",i,"_0.5.Rdata",sep=""))
  predVCM3[,i] <- predict(model, newdata=indiaTestOrg)
  pred3[i,3] <- loss(indiaTest$stunting, predVCM3[indiaTestInd==0,i], tau=0.5)
  
  #---------------
  # Prediction Stumps
  load(paste("stumps/model",i,"_0.05.Rdata",sep=""))
  
  # -> The prediction on the test data was directly calculated after
  # calculating the model. Therefore not any more necessary here
  pred1[i,4] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.05)

  load(paste("stumps/model",i,"_0.1.Rdata",sep=""))
  pred2[i,4] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.1)

  load(paste("stumps/model",i,"_0.5.Rdata",sep=""))
  pred3[i,4] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.5)


  load(paste("trees/model",i,"_0.05.Rdata",sep=""))
  pred1[i,5] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.05)

  load(paste("trees/model",i,"_0.1.Rdata",sep=""))
  pred2[i,5] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.1)

  load(paste("trees/model",i,"_0.5.Rdata",sep=""))
  pred3[i,5] <- loss(indiaTest$stunting, pred[indiaTestInd==0], tau=0.5)


}


save(pred1, pred2, pred3, file="CVRisks.Rdata")
save(predAdditive1, predAdditive2, predAdditive3, file="PredAdditive.Rdata")
save(predVCM1, predVCM2, predVCM3, file="PredVCM.Rdata")

#load("CVRisks.Rdata")

# Change the indices according to performance
pred1 <- pred1[,c(2,1,5,4,3)]
pred2 <- pred2[,c(2,1,5,4,3)]
pred3 <- pred3[,c(2,1,5,4,3)]



#-----------------------------------------------------------
# Boxplots for Empirical Risks 
#-----------------------------------------------------------

#pdf("IndiaRisksAll.pdf", width=8, height=6.5)

cv <- data.frame(error = c(as.vector(pred1),
                           as.vector(pred2),
                           as.vector(pred3)),
                 algo = rep(factor(colnames(pred1), levels = colnames(pred1),
                            labels = c("additive", "rqss", "trees", "stumps", "VCM")), 
                            rep(nrow(pred1), ncol(pred1))),
                 quant = rep(factor(c("5% Quantile", "10% Quantile", "Median"),
                             levels = c("5% Quantile", "10% Quantile", "Median"),
                             labels = c("5% Quantile", "10% Quantile", "Median")),
                             rep(length(pred1), 3)))

mypanel <- function(x, y, ...) {
    panel.bwplot(x = x, y = y, ...)
    z <- as.data.frame(matrix(y, ncol = 5))
    apply(z, 1, function(x)
        llines(x = 1:5, y = x, col = rgb(0.1, 0.1, 0.1, 0.1)))
}
print(bwplot(error ~ algo | quant, data = cv, ylab = "CV risk", 
    scales = list(y = list(relation = "sliced"), x = list(rot = 25)), 
    layout = c(3, 1), panel = mypanel, fill="grey", pch = "|",
    par.settings = list(plot.symbol = list(col=1,pch=20, cex=0.7),
                        box.rectangle = list(col=1),
                        box.umbrella = list(lty=1, col=1)),
    strip= strip.custom(bg="white")))
#dev.off()


#-----------------------------------------------------------
# Boxplots for Optimized Lambda Parameters rqss
#-----------------------------------------------------------

lims <- range(c(lambda1rqss, lambda2rqss, lambda3rqss))

#pdf("IndiaLambdas.pdf", width=8, height=6)
par(mfrow=c(1,3), las=3, mar=c(8,4,2,1))
boxplot(lambda1rqss[,1:6],axes=FALSE, ylab="lambda", xlab="",  main="0.05", ylim=range(lims))
box()
axis(2)
axis(1, at=1:6, labels=c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner"))

boxplot(lambda2rqss[,1:6],axes=FALSE, ylab="lambda", xlab="",  main="0.1", ylim=range(lims))
box()
axis(2)
axis(1, at=1:6, labels=c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner"))

boxplot(lambda3rqss[,1:6],axes=FALSE, ylab="lambda", xlab="",  main="0.5", ylim=range(lims))
box()
axis(2)
axis(1, at=1:6, labels=c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner"))

#dev.off()


# Number of observations in the test data
#pdf("IndiaNTest.pdf")
boxplot(nTest)
#dev.off()

print(summary(nTest))


#-----------------------------------------------------------
# Boxplots for empirical risks without rqss
#-----------------------------------------------------------

load("PredAdditive.Rdata")
load("PredVCM.Rdata")

inds <- paste("w", 1:50, sep="")
quantiles <- c(0.05, 0.1, 0.5)

pred1 <- pred2 <- pred3 <- matrix(0, ncol=4, nrow=length(inds))
colnames(pred1) <- colnames(pred2) <- colnames(pred3) <- c("add","vcm","stump","tree")

for(i in 1:length(inds)){
  print(i)
  
  indiaTrain <- india[india[,inds[i]]==1,]
  
  indiaTest <- india[india[,inds[i]]==3,] 
  
  indiaTest <- cbind(indiaTest, ff = 0)
  names(indiaTest)[ncol(indiaTest)] <- "FALSE"
  
  #---------------
  # Prediction additive
  pred1[i,1] <- loss(indiaTest$stunting, predAdditive1[,i], tau=0.05)
  pred2[i,1] <- loss(indiaTest$stunting, predAdditive2[,i], tau=0.1)   
  pred3[i,1] <- loss(indiaTest$stunting, predAdditive3[,i], tau=0.5)
  
  #---------------
  # Prediction VCM
  pred1[i,2] <- loss(indiaTest$stunting, predVCM1[,i], tau=0.05)
  pred2[i,2] <- loss(indiaTest$stunting, predVCM2[,i], tau=0.1)
  pred3[i,2] <- loss(indiaTest$stunting, predVCM3[,i], tau=0.5)
  
  #---------------
  # Prediction Stumps
  load(paste("stumps/model",i,"_0.05.Rdata",sep=""))
  pred1[i,3] <- loss(indiaTest$stunting, pred, tau=0.05)

  load(paste("stumps/model",i,"_0.1.Rdata",sep=""))
  pred2[i,3] <- loss(indiaTest$stunting, pred, tau=0.1)

  load(paste("stumps/model",i,"_0.5.Rdata",sep=""))
  pred3[i,3] <- loss(indiaTest$stunting, pred, tau=0.5)

  load(paste("trees/model",i,"_0.05.Rdata",sep=""))
  pred1[i,4] <- loss(indiaTest$stunting, pred, tau=0.05)

  load(paste("trees/model",i,"_0.1.Rdata",sep=""))
  pred2[i,4] <- loss(indiaTest$stunting, pred, tau=0.1)

  load(paste("trees/model",i,"_0.5.Rdata",sep=""))
  pred3[i,4] <- loss(indiaTest$stunting, pred, tau=0.5)


}

save(pred1, pred2, pred3, file="CVRisksTotal.Rdata")

pred1 <- pred1[,c(1,4,3,2)]
pred2 <- pred2[,c(1,4,3,2)]
pred3 <- pred3[,c(1,4,3,2)]

#pdf("IndiaRisksWithoutRqss.pdf", width=8, height=6.5)

cv <- data.frame(error = c(as.vector(pred1),
                           as.vector(pred2),
                           as.vector(pred3)),
                 algo = rep(factor(colnames(pred1), levels = colnames(pred1),
                            labels = c("additive", "trees", "stumps", "VCM")), 
                            rep(nrow(pred1), ncol(pred1))),
                 quant = rep(factor(c("5% Quantile", "10% Quantile", "Median"),
                             levels = c("5% Quantile", "10% Quantile", "Median"),
                             labels = c("5% Quantile", "10% Quantile", "Median")),
                             rep(length(pred1), 3)))

mypanel <- function(x, y, ...) {
    panel.bwplot(x = x, y = y, ...)
    z <- as.data.frame(matrix(y, ncol = 4))
    apply(z, 1, function(x)
        llines(x = 1:4, y = x, col = rgb(0.1, 0.1, 0.1, 0.1)))
}
print(bwplot(error ~ algo | quant, data = cv, ylab = "CV risk", 
    scales = list(y = list(relation = "sliced"), x = list(rot = 25)), 
    layout = c(3, 1), panel = mypanel, fill="grey", pch = "|",
    par.settings = list(plot.symbol = list(col=1,pch=20, cex=0.7),
                        box.rectangle = list(col=1),
                        box.umbrella = list(lty=1, col=1)),
    strip= strip.custom(bg="white")))
#dev.off()



