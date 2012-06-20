
load("plotresults.Rdata")

mystabsel <- function(mat, q=10, FWER=0.05, bl)
  {
  p <- length(bl)
  cutoff <- min(0.9, (q^2/(FWER * p) + 1)/2)
  topq <- rep(0, ncol(mat))
  for(i in 1:ncol(mat))
    {
    print(i)
    if(length(unique(mat[!is.na(mat[,i]),i])) < q)
      topq[i] <- sum(!is.na(mat[,i]))
    else 
      {
      topq[i] <- q
      while(length(unique(mat[1:(topq[i]),i])) < q)
        {
        topq[i] <- topq[i]+1
        }
      }
    }
  topbl <- rep(0, p)
  for(i in 1:p)
    {
    for(j in 1:ncol(mat))
      {
      if(i %in% mat[1:(topq[j]),j])
        topbl[i] <- topbl[i] + 1
      }
    }  
  topbl <- topbl / ncol(mat)
  topblind <- 1*(topbl >= cutoff)
  return(list(topbl=topbl, topblind=topblind))
  }

# additive models
q <- 10
FWER <- 0.05
bl <- levels(add1.selected[[1]])

add1.topbl <- mystabsel(add1.stabsel, q=q, FWER=FWER, bl=bl)
add2.topbl <- mystabsel(add2.stabsel, q=q, FWER=FWER, bl=bl)
add3.topbl <- mystabsel(add3.stabsel, q=q, FWER=FWER, bl=bl)

bl[add1.topbl$topblind==1]
bl[add2.topbl$topblind==1]
bl[add3.topbl$topblind==1]

# VCMs

q <- 10
FWER <- 0.05
bl <- levels(vcm1.selected[[1]])

vcm1.topbl <- mystabsel(vcm1.stabsel, q=q, FWER=FWER, bl=bl)
vcm2.topbl <- mystabsel(vcm2.stabsel, q=q, FWER=FWER, bl=bl)
vcm3.topbl <- mystabsel(vcm3.stabsel, q=q, FWER=FWER, bl=bl)

bl[vcm1.topbl$topblind==1]
bl[vcm2.topbl$topblind==1]
bl[vcm3.topbl$topblind==1]

# rescale results to match definition of the Z-Score

plotdata.add1[[1]][,-1] <- plotdata.add1[[1]][,-1]/100
plotdata.add3[[1]][,-1] <- plotdata.add3[[1]][,-1]/100

plotdata.add1[[2]][,-1] <- plotdata.add1[[2]][,-1]/100
plotdata.add3[[2]][,-1] <- plotdata.add3[[2]][,-1]/100

plotdata.add1[[4]][,-1] <- plotdata.add1[[4]][,-1]/100
plotdata.add3[[4]][,-1] <- plotdata.add3[[4]][,-1]/100

plotdata.add1[[7]][,-1] <- plotdata.add1[[7]][,-1]/100
plotdata.add3[[7]][,-1] <- plotdata.add3[[7]][,-1]/100

plotdata.add1[[8]][,-1] <- plotdata.add1[[8]][,-1]/100
plotdata.add3[[8]][,-1] <- plotdata.add3[[8]][,-1]/100

plotdata.add1[[14]][,-1] <- plotdata.add1[[14]][,-1]/100
plotdata.add3[[14]][,-1] <- plotdata.add3[[14]][,-1]/100

plotdata.vcm1[[1]][,-1] <- plotdata.vcm1[[1]][,-1]/100
plotdata.vcm3[[1]][,-1] <- plotdata.vcm3[[1]][,-1]/100

plotdata.vcm1[[length(varnames)+1]][,-1] <- plotdata.vcm1[[length(varnames)+1]][,-1]/100
plotdata.vcm3[[length(varnames)+1]][,-1] <- plotdata.vcm3[[length(varnames)+1]][,-1]/100

# parallel coordinate plots of the prediction risks

pdf("india_predictions.pdf", height = 4)

library("lattice")
cv <- data.frame(error = c(as.vector(pred1),
                           as.vector(pred2),
                           as.vector(pred3)),
                 algo = rep(factor(colnames(pred1), levels = colnames(pred1),
                            labels = c("additive", "VCM", "stumps", "trees")), 
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

dev.off()

# nonlinear effects from the aditive model shown in the main part of the paper

pdf("india_additive_nonlinear.pdf", width=8, height=10.5)
mat <- matrix(0, 6, 4)
mat[1,] <- c(0, 0.5, 2/3, 1)
mat[2,] <- c(0.5, 1, 2/3, 1)
mat[3,] <- c(0, 0.5, 1/3, 2/3)
mat[4,] <- c(0.5, 1, 1/3, 2/3)
mat[5,] <- c(0, 0.5, 0, 1/3)
mat[6,] <- c(0.5, 1, 0, 1/3)
split.screen(mat)

i <- 1
screen(1)
par(mai=c(1.0,0.5,0.5,0.1))

col <- rgb(0.1, 0.1, 0.1, 0.1)
    pr <- range(plotdata.add1[[i]][,-1])
    ind <- order(plotdata.add1[[i]][,1])
    plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="5% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add1[[i]]))
      lines(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], col=col)
 
screen(2)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add3[[i]][,-1])
    ind <- order(plotdata.add3[[i]][,1])
    plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="Median",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add3[[i]]))
      lines(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], col=col)

i <- 2
screen(3)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add1[[i]][,-1])
    ind <- order(plotdata.add1[[i]][,1])
    plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="duration of breastfeeding in months", ylab="", main="5% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add1[[i]]))
      lines(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], col=col)
 
screen(4)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add3[[i]][,-1])
    ind <- order(plotdata.add3[[i]][,1])
    plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="duration of breastfeeding in months", ylab="", main="Median",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add3[[i]]))
      lines(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], col=col)

i <- 4
screen(5)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add1[[i]][,-1])
    ind <- order(plotdata.add1[[i]][,1])
    plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the mother in years", ylab="", main="5% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add1[[i]]))
      lines(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], col=col)
 
screen(6)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add3[[i]][,-1])
    ind <- order(plotdata.add3[[i]][,1])
    plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the mother in years", ylab="", main="Median",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add3[[i]]))
      lines(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], col=col)

close.screen(all = TRUE)
dev.off()

# parametric effects from the additive model shown in the main part of the paper

pdf("india_additive_categorical.pdf", width=8, height=10.5)
mat <- matrix(0, 6, 4)
mat[1,] <- c(0, 0.5, 2/3, 1)
mat[2,] <- c(0.5, 1, 2/3, 1)
mat[3,] <- c(0, 0.5, 1/3, 2/3)
mat[4,] <- c(0.5, 1, 1/3, 2/3)
mat[5,] <- c(0, 0.5, 0, 1/3)
mat[6,] <- c(0.5, 1, 0, 1/3)
split.screen(mat)

i <- 7
screen(1)
par(mai=c(0.5,0.5,0.5,0.1))

pr <- range(plotdata.add1[[i]][,-1])
ind <- order(plotdata.add1[[i]][,1])

tmp <- t(plotdata.add1[[i]][,-1])
colnames(tmp) <- levels(plotdata.add1[[i]][,1])
boxplot(tmp, main="5% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

screen(2)
par(mai=c(0.5,0.5,0.5,0.1))

tmp <- t(plotdata.add3[[i]][,-1])
colnames(tmp) <- levels(plotdata.add3[[i]][,1])
boxplot(tmp, main="Median", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

i <- 8
screen(3)
par(mai=c(0.5,0.5,0.5,0.1))

tmp <- t(plotdata.add1[[i]][,-1])
colnames(tmp) <- levels(plotdata.add1[[i]][,1])
boxplot(tmp, main="5% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

screen(4)
par(mai=c(0.5,0.5,0.5,0.1))

tmp <- t(plotdata.add3[[i]][,-1])
colnames(tmp) <- levels(plotdata.add3[[i]][,1])
boxplot(tmp, main="Median", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

i <- 14
screen(5)
par(mai=c(0.5,0.5,0.5,0.1))

tmp <- t(plotdata.add1[[i]][,-1])
colnames(tmp) <- as.character(plotdata.add1[[i]][,1])
tmp <- tmp[, levels(plotdata.add1[[i]][,1])]
boxplot(tmp, main="5% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

screen(6)
par(mai=c(0.5,0.5,0.5,0.1))

tmp <- t(plotdata.add3[[i]][,-1])
colnames(tmp) <- as.character(plotdata.add3[[i]][,1])
tmp <- tmp[, levels(plotdata.add3[[i]][,1])]
boxplot(tmp, main="Median", cex.lab=1., cex.axis=0.9, col = "grey", 
        pars = list(outcol = "black", outpch = 19, outcex = 0.7))
out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

close.screen(all = TRUE)
dev.off()

# nonlinear effects from the VCM shown in the main paper

pdf("india_vcm_nonlinear.pdf", width=8, height=7)
mat <- matrix(0, 4, 4)
mat[1,] <- c(0, 0.5, 0.5, 1)
mat[2,] <- c(0.5, 1, 0.5, 1)
mat[3,] <- c(0, 0.5, 0, 0.5)
mat[4,] <- c(0.5, 1, 0, 0.5)
split.screen(mat)

i <- 1
screen(1)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(cbind(plotdata.vcm1[[i]][,-1],
                  plotdata.vcm1[[length(varnames)+i]][,-1]))
    ind <- order(plotdata.vcm1[[i]][,1])
    plot(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="5% Quantile (boys)",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.vcm1[[i]]))
      lines(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,j], col=col)

 
screen(2)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(cbind(plotdata.vcm3[[i]][,-1],
                  plotdata.vcm3[[length(varnames)+i]][,-1]))
    ind <- order(plotdata.vcm3[[i]][,1])
    plot(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="Median (boys)",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.vcm3[[i]]))
      lines(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,j], col=col)

screen(3)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(cbind(plotdata.vcm3[[i]][,-1],
                  plotdata.vcm3[[length(varnames)+i]][,-1]))
    pr <- range(plotdata.vcm1[[length(varnames)+i]][,-1])
    ind <- order(plotdata.vcm1[[length(varnames)+i]][,1])
    plot(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="5% Quantile (girls)",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.vcm1[[length(varnames)+i]]))
      lines(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,j], col=col)
 
screen(4)
par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(cbind(plotdata.vcm3[[i]][,-1],
                   plotdata.vcm3[[length(varnames)+i]][,-1]))
    pr <- range(plotdata.vcm3[[length(varnames)+i]][,-1])
    ind <- order(plotdata.vcm3[[length(varnames)+i]][,1])
    plot(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="age of the child in months", ylab="", main="Median (girls)",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.vcm3[[length(varnames)+i]]))
      lines(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,j], col=col)

close.screen(all = TRUE)
dev.off()

# Figures for the additive model from the eSupplement

load("plotresults.Rdata")

for(i in 1:21)
  {
  plotdata.add1[[i]][,-1] <- plotdata.add1[[i]][,-1]/100
  plotdata.add2[[i]][,-1] <- plotdata.add2[[i]][,-1]/100
  plotdata.add3[[i]][,-1] <- plotdata.add3[[i]][,-1]/100
  }
  
varlab <- c("age of the child in months",
"duration of breastfeeding in months",
"body mass index of the mother",
"age of the mother in years",
"education years of the mother",
"education years of the mother's partner",
"sex of the child",
"twin birth",
"position in the birth order",
"unemployment of the mother",
"religion of the mother",
"area of residence",
"no. of dead children",
"wealth index",
"electrictiy",
"radio",
"television",
"refrigerator",
"bicycle",
"motorcycle",
"car")

plotfnc <- function(offset=1)
  {
  mat <- matrix(0, 6, 4)
  mat[1,] <- c(0, 0.5, 2/3, 1)
  mat[2,] <- c(0.5, 1, 2/3, 1)
  mat[3,] <- c(0, 0.5, 1/3, 2/3)
  mat[4,] <- c(0.5, 1, 1/3, 2/3)
  mat[5,] <- c(0, 0.5, 0, 1/3)
  mat[6,] <- c(0.5, 1, 0, 1/3)
  split.screen(mat)

  if(is.numeric(plotdata.add1[[offset]][,1]))
    {
    col <- rgb(0.1, 0.1, 0.1, 0.1)

    screen(1)
    par(mai=c(1.0,0.5,0.5,0.1))
  
    pr <- range(plotdata.add1[[offset]][,-1])
    ind <- order(plotdata.add1[[offset]][,1])

    plot(plotdata.add1[[offset]][ind,1], plotdata.add1[[offset]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset], ylab="", main="5% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add1[[offset]]))
      lines(plotdata.add1[[offset]][ind,1], plotdata.add1[[offset]][ind,j], col=col)
  
    screen(3)
    par(mai=c(1.0,0.5,0.5,0.1))
     
    pr <- range(plotdata.add2[[offset]][,-1])
    ind <- order(plotdata.add2[[offset]][,1])
  
    plot(plotdata.add2[[offset]][ind,1], plotdata.add2[[offset]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset], ylab="", main="10% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add2[[offset]]))
      lines(plotdata.add2[[offset]][ind,1], plotdata.add2[[offset]][ind,j], col=col)
  
    screen(5)
    par(mai=c(1.0,0.5,0.5,0.1))
    
    pr <- range(plotdata.add3[[offset]][,-1])
    ind <- order(plotdata.add3[[offset]][,1])
  
    plot(plotdata.add3[[offset]][ind,1], plotdata.add3[[offset]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset], ylab="", main="Median",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add3[[offset]]))
      lines(plotdata.add3[[offset]][ind,1], plotdata.add3[[offset]][ind,j], col=col)
    }
  else
    {
    col <- rgb(0.1, 0.1, 0.1, 0.1)

    screen(1)
    par(mai=c(1.0,0.5,0.5,0.1))
  
    pr <- range(plotdata.add1[[offset]][,-1])
    ind <- order(plotdata.add1[[offset]][,1])

    tmp <- t(plotdata.add1[[offset]][,-1])
    colnames(tmp) <- as.character(plotdata.add1[[offset]][,1])
    tmp <- tmp[, levels(plotdata.add1[[offset]][,1])]
    boxplot(tmp, main="5% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))
  
    screen(3)
    par(mai=c(1.0,0.5,0.5,0.1))
     
    pr <- range(plotdata.add2[[offset]][,-1])
    ind <- order(plotdata.add2[[offset]][,1])
  
    tmp <- t(plotdata.add2[[offset]][,-1])
    colnames(tmp) <- as.character(plotdata.add2[[offset]][,1])
    tmp <- tmp[, levels(plotdata.add2[[offset]][,1])]
    boxplot(tmp, main="10% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))
  
    screen(5)
    par(mai=c(1.0,0.5,0.5,0.1))
    
    pr <- range(plotdata.add3[[offset]][,-1])
    ind <- order(plotdata.add3[[offset]][,1])
  
    tmp <- t(plotdata.add3[[offset]][,-1])
    colnames(tmp) <- as.character(plotdata.add3[[offset]][,1])
    tmp <- tmp[, levels(plotdata.add3[[offset]][,1])]
    boxplot(tmp, main="Median", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))
    }

  if(is.numeric(plotdata.add1[[offset+1]][,1]))
    {
    col <- rgb(0.1, 0.1, 0.1, 0.1)

    screen(2)
    par(mai=c(1.0,0.5,0.5,0.1))
  
    pr <- range(plotdata.add1[[offset+1]][,-1])
    ind <- order(plotdata.add1[[offset+1]][,1])

    plot(plotdata.add1[[offset+1]][ind,1], plotdata.add1[[offset+1]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset+1], ylab="", main="5% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add1[[offset+1]]))
      lines(plotdata.add1[[offset+1]][ind,1], plotdata.add1[[offset+1]][ind,j], col=col)
  
    screen(4)
    par(mai=c(1.0,0.5,0.5,0.1))
     
    pr <- range(plotdata.add2[[offset+1]][,-1])
    ind <- order(plotdata.add2[[offset+1]][,1])
  
    plot(plotdata.add2[[offset+1]][ind,1], plotdata.add2[[offset+1]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset+1], ylab="", main="10% Quantile",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add2[[offset+1]]))
      lines(plotdata.add2[[offset+1]][ind,1], plotdata.add2[[offset+1]][ind,j], col=col)
  
    screen(6)
    par(mai=c(1.0,0.5,0.5,0.1))
    
    pr <- range(plotdata.add3[[offset+1]][,-1])
    ind <- order(plotdata.add3[[offset+1]][,1])
  
    plot(plotdata.add3[[offset+1]][ind,1], plotdata.add3[[offset+1]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab=varlab[offset+1], ylab="", main="Median",
         col=col, cex.lab=1.)
    for(j in 2:ncol(plotdata.add3[[offset+1]]))
      lines(plotdata.add3[[offset+1]][ind,1], plotdata.add3[[offset+1]][ind,j], col=col)
    }
  else
    {
    col <- rgb(0.1, 0.1, 0.1, 0.1)

    screen(2)
    par(mai=c(1.0,0.5,0.5,0.1))

    pr <- range(plotdata.add1[[offset+1]][,-1])
    ind <- order(plotdata.add1[[offset+1]][,1])

    tmp <- t(plotdata.add1[[offset+1]][,-1])
    colnames(tmp) <- as.character(plotdata.add1[[offset+1]][,1])
    tmp <- tmp[, levels(plotdata.add1[[offset+1]][,1])]
    boxplot(tmp, main="5% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset+1])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))

    screen(4)
    par(mai=c(1.0,0.5,0.5,0.1))
    
    pr <- range(plotdata.add2[[offset+1]][,-1])
    ind <- order(plotdata.add2[[offset+1]][,1])
  
    tmp <- t(plotdata.add2[[offset+1]][,-1])
    colnames(tmp) <- as.character(plotdata.add2[[offset+1]][,1])
    tmp <- tmp[, levels(plotdata.add2[[offset+1]][,1])]
    boxplot(tmp, main="10% Quantile", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset+1])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))
  
    screen(6)
    par(mai=c(1.0,0.5,0.5,0.1))
    
    pr <- range(plotdata.add3[[offset+1]][,-1])
    ind <- order(plotdata.add3[[offset+1]][,1])
  
    tmp <- t(plotdata.add3[[offset+1]][,-1])
    colnames(tmp) <- as.character(plotdata.add3[[offset+1]][,1])
    tmp <- tmp[, levels(plotdata.add3[[offset+1]][,1])]
    boxplot(tmp, main="Median", cex.lab=1., cex.axis=0.9, col = "grey", 
            pars = list(outcol = "black", outpch = 19, outcex = 0.7), xlab=varlab[offset+1])
    out <- apply(tmp, 1, function(x) lines(1:ncol(tmp), x, col = col))
    }
    
  close.screen(all = TRUE)
  }

pdf("india_additive_appendix1.pdf", width=8, height=10.5)
plotfnc(1)
dev.off()

pdf("india_additive_appendix2.pdf", width=8, height=10.5)
plotfnc(3)
dev.off()

pdf("india_additive_appendix3.pdf", width=8, height=10.5)
plotfnc(5)
dev.off()

pdf("india_additive_appendix4.pdf", width=8, height=10.5)
plotfnc(7)
dev.off()

pdf("india_additive_appendix5.pdf", width=8, height=10.5)
plotfnc(9)
dev.off()

pdf("india_additive_appendix6.pdf", width=8, height=10.5)
plotfnc(11)
dev.off()

pdf("india_additive_appendix7.pdf", width=8, height=10.5)
plotfnc(13)
dev.off()

pdf("india_additive_appendix8.pdf", width=8, height=10.5)
plotfnc(15)
dev.off()

pdf("india_additive_appendix9.pdf", width=8, height=10.5)
plotfnc(17)
dev.off()

pdf("india_additive_appendix10.pdf", width=8, height=10.5)
plotfnc(19)
dev.off()
  
pdf("india_additive_appendix11.pdf", width=8, height=10.5)
plotfnc(21)
dev.off()
