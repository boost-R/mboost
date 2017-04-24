
library("mboost")
load("india.Rdata")

inds <- paste("w", 1:50, sep="")
quantiles <- c(0.05, 0.1, 0.5)
pred1 <- pred2 <- pred3 <- matrix(0, ncol=4, nrow=length(inds))
colnames(pred1) <- colnames(pred2) <- colnames(pred3) <- c("add","vcm","stump","tree")
add1.selected <- add2.selected <- add3.selected <- list()
vcm1.selected <- vcm2.selected <- vcm3.selected <- list()
add1.stabsel <- add2.stabsel <- add3.stabsel <- matrix(NA, ncol=length(inds), nrow=210000)
add1.its <- add2.its <- add3.its <- rep(0, length(inds))
vcm1.its <- vcm2.its <- vcm3.its <- rep(0, length(inds))
vcm1.stabsel <- vcm2.stabsel <- vcm3.stabsel <- matrix(NA, ncol=length(inds), nrow=210000)

loss <- function(y, f, tau)
{
    sum(tau*(y - f)*((y - f) >= 0) - (1-tau)*(y - f)*((y - f) < 0))
}

# read the predictions

for(i in 1:length(inds))
{
    print(i)
    indiapred <- india[india[,inds[i]]==3,]
    indiapred <- cbind(indiapred, ff = 0)
    names(indiapred)[71] <- "FALSE"
    
    load(paste("additive/model",i,"_0.05.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred1[i,1] <- loss(indiapred$stunting, pred, tau=0.05)
    add1.selected[[i]] <- factor(selected(model), levels=1:28, labels=names(model$baselearner))
    add1.stabsel[1:mstop(model),i] <- selected(model) 
    add1.its[i] <- mstop(model)
    
    load(paste("additive/model",i,"_0.1.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred2[i,1] <- loss(indiapred$stunting, pred, tau=0.1)
    add2.selected[[i]] <- factor(selected(model), levels=1:28, labels=names(model$baselearner))
    add2.stabsel[1:mstop(model),i] <- selected(model) 
    add2.its[i] <- length(selected(model))
    
    load(paste("additive/model",i,"_0.5.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred3[i,1] <- loss(indiapred$stunting, pred, tau=0.5)
    add3.selected[[i]] <- factor(selected(model), levels=1:28, labels=names(model$baselearner))
    add3.stabsel[1:mstop(model),i] <- selected(model) 
    add3.its[i] <- length(selected(model))
    
    
    load(paste("VCM/model",i,"_0.05.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred1[i,2] <- loss(indiapred$stunting, pred, tau=0.05)
    vcm1.selected[[i]] <- factor(selected(model), levels=1:54, labels=names(model$baselearner))
    vcm1.stabsel[1:mstop(model),i] <- selected(model) 
    vcm1.its[i] <- length(selected(model))
    
    load(paste("VCM/model",i,"_0.1.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred2[i,2] <- loss(indiapred$stunting, pred, tau=0.1)
    vcm2.selected[[i]] <- factor(selected(model), levels=1:54, labels=names(model$baselearner))
    vcm2.stabsel[1:mstop(model),i] <- selected(model) 
    vcm2.its[i] <- length(selected(model))
    
    load(paste("VCM/model",i,"_0.5.Rdata",sep=""))
    pred <- predict(model, newdata=indiapred)
    pred3[i,2] <- loss(indiapred$stunting, pred, tau=0.5)
    vcm3.selected[[i]] <- factor(selected(model), levels=1:54, labels=names(model$baselearner))
    vcm3.stabsel[1:mstop(model),i] <- selected(model) 
    vcm3.its[i] <- length(selected(model))
    
    
    load(paste("stumps/model",i,"_0.05.Rdata",sep=""))
    pred1[i,3] <- loss(indiapred$stunting, pred, tau=0.05)
    
    load(paste("stumps/model",i,"_0.1.Rdata",sep=""))
    pred2[i,3] <- loss(indiapred$stunting, pred, tau=0.1)
    
    load(paste("stumps/model",i,"_0.5.Rdata",sep=""))
    pred3[i,3] <- loss(indiapred$stunting, pred, tau=0.5)
    
    
    load(paste("trees/model",i,"_0.05.Rdata",sep=""))
    pred1[i,4] <- loss(indiapred$stunting, pred, tau=0.05)
    
    load(paste("trees/model",i,"_0.1.Rdata",sep=""))
    pred2[i,4] <- loss(indiapred$stunting, pred, tau=0.1)
    
    load(paste("trees/model",i,"_0.5.Rdata",sep=""))
    pred3[i,4] <- loss(indiapred$stunting, pred, tau=0.5)
}



# plotting of effects in the additive model

preddata <- preddatafm <- list()
plotdata.add1 <- plotdata.add2 <- plotdata.add3 <- list()
plotdata.vcm1 <- plotdata.vcm2 <- plotdata.vcm3 <- list()

varnames <- c("cage", "breastfeeding", "mbmi", "mage", "medu", "edupartner",
              "csex", "ctwin", "cbirthorder", "munemployed", "mreligion",
              "mresidence", "deadchildren", "wealth", "electricity", "radio",
              "television", "refrigerator", "bicycle", "motorcycle", "car")

refval <- list()
for(i in 1:length(varnames))
{
    print(i)
    if(is.numeric(india[,varnames[i]]))
        refval[[i]] <- mean(india[,varnames[i]])
    if(is.factor(india[,varnames[i]]))
        refval[[i]] <- unique(india[,varnames[i]])[which.max(table(india[,varnames[i]]))]
}

for(i in 1:length(varnames))
{
    print(i)
    vals <- unique(india[,varnames[i]])
    preddata[[i]] <- india[1:length(vals),]
    for(j in 1:length(varnames))
    {
        if(j==i)
            preddata[[i]][,varnames[j]] <- vals
        else
            preddata[[i]][,varnames[j]] <- rep(refval[[j]], length(vals))
    }
    preddata[[i]] <- cbind(preddata[[i]], ff = 0)
    names(preddata[[i]])[71] <- "FALSE"
    
    preddatafm[[i]] <- preddata[[i]]
    preddatafm[[i]][,"csex"] <- rep(unique(india[,"csex"])[which.min(table(india[,"csex"]))], length(vals))
    
    
    plotdata.add1[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.add1[[i]]) <- c(varnames[i],paste("f", 1:length(inds), sep=""))
    
    plotdata.add2[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.add1[[i]]) <- c(varnames[i],paste("f", 1:length(inds), sep=""))
    
    plotdata.add3[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.add1[[i]]) <- c(varnames[i],paste("f", 1:length(inds), sep=""))
    
    
    plotdata.vcm1[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm1[[i]]) <- c(varnames[i],paste("fm", 1:length(inds), sep=""))
    plotdata.vcm1[[length(varnames)+i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm1[[length(varnames)+i]]) <- c(varnames[i],paste("ff", 1:length(inds), sep=""))
    
    plotdata.vcm2[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm2[[i]]) <- c(varnames[i],paste("fm", 1:length(inds), sep=""))
    plotdata.vcm2[[length(varnames)+i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm2[[length(varnames)+i]]) <- c(varnames[i],paste("ff", 1:length(inds), sep=""))
    
    plotdata.vcm3[[i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm3[[i]]) <- c(varnames[i],paste("fm", 1:length(inds), sep=""))
    plotdata.vcm3[[length(varnames)+i]] <- data.frame(vals, matrix(0, nrow=length(vals), ncol=length(inds)))
    names(plotdata.vcm3[[length(varnames)+i]]) <- c(varnames[i],paste("ff", 1:length(inds), sep=""))
}

for(i in 1:length(inds))
{
    print(i)
    load(paste("additive/model",i,"_0.05.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.add1[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
    }
    
    load(paste("additive/model",i,"_0.1.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.add2[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
    }
    
    load(paste("additive/model",i,"_0.5.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.add3[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
    }
    
    
    load(paste("VCM/model",i,"_0.05.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.vcm1[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
        plotdata.vcm1[[length(varnames)+j]][,(i+1)] <- predict(model, newdata=preddatafm[[j]])
    }
    
    load(paste("VCM/model",i,"_0.1.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.vcm2[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
        plotdata.vcm2[[length(varnames)+j]][,(i+1)] <- predict(model, newdata=preddatafm[[j]])
    }
    
    load(paste("VCM/model",i,"_0.5.Rdata",sep=""))
    for(j in 1:length(varnames))
    {
        plotdata.vcm3[[j]][,(i+1)] <- predict(model, newdata=preddata[[j]])
        plotdata.vcm3[[length(varnames)+j]][,(i+1)] <- predict(model, newdata=preddatafm[[j]])
    }
}

save("pred1", "pred2", "pred3", "plotdata.add1", "plotdata.add2",
     "plotdata.add3", "plotdata.vcm1", "plotdata.vcm2", "plotdata.vcm3",
     "add1.selected", "add2.selected", "add3.selected", 
     "vcm1.selected", "vcm2.selected", "vcm3.selected", 
     "add1.its", "add2.its", "add3.its", "vcm1.its", "vcm2.its", "vcm3.its", 
     "add1.stabsel", "add2.stabsel", "add3.stabsel", 
     "vcm1.stabsel", "vcm2.stabsel", "vcm3.stabsel", 
     "varnames", file="plotresults.Rdata")