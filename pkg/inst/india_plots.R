
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


pdf("india_plots.pdf")

# parallel coordinate plots of the prediction risks

pr <- range(pred1)
ind <- 1:4
plot(ind, pred1[1,ind], ylim=pr, axes=FALSE, type="l", xlab="", ylab="CV risk",
     col=gray(0.7), cex.lab=1.4, main="0.05")
box()
axis(2)
axis(1, at=ind, labels=c("additive", "VCM", "stumps", "trees")[ind])
for(i in 2:nrow(pred1))
  lines(ind, pred1[i,ind], col=gray(0.7))


pr <- range(pred2)
ind <- 1:4
plot(ind, pred2[1,ind], ylim=pr, axes=FALSE, type="l", xlab="", ylab="CV risk",
     col=gray(0.7), cex.lab=1.4, main="0.1")
box()
axis(2)
axis(1, at=ind, labels=c("additive", "VCM", "stumps", "trees")[ind])
for(i in 2:nrow(pred2))
  lines(ind, pred2[i,ind], col=gray(0.7))


pr <- range(pred3)
ind <- 1:4
plot(ind, pred3[1,ind], ylim=pr, axes=FALSE, type="l", xlab="", ylab="CV risk",
     col=gray(0.7), cex.lab=1.4, main="0.5")
box()
axis(2)
axis(1, at=ind, labels=c("additive", "VCM", "stumps", "trees")[ind])
for(i in 2:nrow(pred3))
  lines(ind, pred3[i,ind], col=gray(0.7))

# parallel coordinate plots of estimated effects

# additive model, tau=0.05:

par(mfrow=c(2,2))
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.add1[[i]][,1]))
    {
    pr <- range(plotdata.add1[[i]][,-1])
    ind <- order(plotdata.add1[[i]][,1])
    plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.05",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add1[[i]]))
      lines(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], col=gray(0.7))
    }
  else
    {
    pr <- range(plotdata.add1[[i]][,-1])
    ind <- order(plotdata.add1[[i]][,1])
    plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.05",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add1[[i]]))
      plot(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.add1[[i]]))
      lines(plotdata.add1[[i]][ind,1], plotdata.add1[[i]][ind,j], lwd=1.5)
    }
  }

par(mfrow=c(2,2))
# additive model, tau=0.1:
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.add2[[i]][,1]))
    {
    pr <- range(plotdata.add2[[i]][,-1])
    ind <- order(plotdata.add2[[i]][,1])
    plot(plotdata.add2[[i]][ind,1], plotdata.add2[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.1",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add2[[i]]))
      lines(plotdata.add2[[i]][ind,1], plotdata.add2[[i]][ind,j], col=gray(0.7))
    }
  else
    {
    pr <- range(plotdata.add2[[i]][,-1])
    ind <- order(plotdata.add2[[i]][,1])
    plot(plotdata.add2[[i]][ind,1], plotdata.add2[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.1",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add2[[i]]))
      plot(plotdata.add2[[i]][ind,1], plotdata.add2[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.add2[[i]]))
      lines(plotdata.add2[[i]][ind,1], plotdata.add2[[i]][ind,j], lwd=1.5)
    }
  }


par(mfrow=c(2,2))
# additive model, tau=0.5:
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.add3[[i]][,1]))
    {
    pr <- range(plotdata.add3[[i]][,-1])
    ind <- order(plotdata.add3[[i]][,1])
    plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.5",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add3[[i]]))
      lines(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], col=gray(0.7))
    }
  else
    {
    pr <- range(plotdata.add3[[i]][,-1])
    ind <- order(plotdata.add3[[i]][,1])
    plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="additive, 0.5",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.add3[[i]]))
      plot(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.add3[[i]]))
      lines(plotdata.add3[[i]][ind,1], plotdata.add3[[i]][ind,j], lwd=1.5)
    }
  }

par(mfrow=c(2,2))
# VCM, tau=0.05:
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.vcm1[[i]][,1]))
    {
    ind <- order(plotdata.vcm1[[i]][,1])
    pr <- range(plotdata.vcm1[[i]][,-1])
    plot(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.05, male",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm1[[i]]))
      lines(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,j], col=gray(0.7))
 
    pr <- range(plotdata.vcm1[[length(varnames)+i]][,-1])
    plot(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.05, female",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm1[[length(varnames)+i]]))
      lines(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,j], col=gray(0.7))
    }
  else
    {
    ind <- order(plotdata.vcm1[[i]][,1])
    pr <- range(plotdata.vcm1[[i]][,-1])
    plot(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.05, male",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm1[[i]]))
      plot(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.vcm1[[i]]))
      lines(plotdata.vcm1[[i]][ind,1], plotdata.vcm1[[i]][ind,j], lwd=1.5)
      
    pr <- range(plotdata.vcm1[[length(varnames)+i]][,-1])
    plot(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.05, female",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm1[[length(varnames)+i]]))
      plot(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,j], add=TRUE, border=gray(0.7))
    for(j in 2:ncol(plotdata.vcm1[[length(varnames)+i]]))
      lines(plotdata.vcm1[[length(varnames)+i]][ind,1], plotdata.vcm1[[length(varnames)+i]][ind,j], lwd=1.5)
    }   
  }

par(mfrow=c(2,2))
# VCM, tau=0.1:
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.vcm2[[i]][,1]))
    {
    ind <- order(plotdata.vcm2[[i]][,1])
    pr <- range(plotdata.vcm2[[i]][,-1])
    plot(plotdata.vcm2[[i]][ind,1], plotdata.vcm2[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.1, male",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm2[[i]]))
      lines(plotdata.vcm2[[i]][ind,1], plotdata.vcm2[[i]][ind,j], col=gray(0.7))
 
    pr <- range(plotdata.vcm2[[length(varnames)+i]][,-1])
    plot(plotdata.vcm2[[length(varnames)+i]][ind,1], plotdata.vcm2[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.1, female",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm2[[length(varnames)+i]]))
      lines(plotdata.vcm2[[length(varnames)+i]][ind,1], plotdata.vcm2[[length(varnames)+i]][ind,j], col=gray(0.7))
    }
  else
    {
    ind <- order(plotdata.vcm2[[i]][,1])
    pr <- range(plotdata.vcm2[[i]][,-1])
    plot(plotdata.vcm2[[i]][ind,1], plotdata.vcm2[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.1, male",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm2[[i]]))
      plot(plotdata.vcm2[[i]][ind,1], plotdata.vcm2[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.vcm2[[i]]))
      lines(plotdata.vcm2[[i]][ind,1], plotdata.vcm2[[i]][ind,j], lwd=1.5)
      
    pr <- range(plotdata.vcm2[[length(varnames)+i]][,-1])
    plot(plotdata.vcm2[[length(varnames)+i]][ind,1], plotdata.vcm2[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.1, female",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm2[[length(varnames)+i]]))
      plot(plotdata.vcm2[[length(varnames)+i]][ind,1], plotdata.vcm2[[length(varnames)+i]][ind,j], add=TRUE, border=gray(0.7))
    for(j in 2:ncol(plotdata.vcm2[[length(varnames)+i]]))
      lines(plotdata.vcm2[[length(varnames)+i]][ind,1], plotdata.vcm2[[length(varnames)+i]][ind,j], lwd=1.5)
    }   
  }

par(mfrow=c(2,2))
# VCM, tau=0.5:
for(i in 1:length(varnames))
  {
  if(is.numeric(plotdata.vcm3[[i]][,1]))
    {
    ind <- order(plotdata.vcm3[[i]][,1])
    pr <- range(plotdata.vcm3[[i]][,-1])
    plot(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.5, male",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm3[[i]]))
      lines(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,j], col=gray(0.7))
 
    pr <- range(plotdata.vcm3[[length(varnames)+i]][,-1])
    plot(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.5, female",
         col=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm3[[length(varnames)+i]]))
      lines(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,j], col=gray(0.7))
    }
  else
    {
    ind <- order(plotdata.vcm3[[i]][,1])
    pr <- range(plotdata.vcm3[[i]][,-1])
    plot(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.5, male",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm3[[i]]))
      plot(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,j], border=gray(0.7), add=TRUE)
    for(j in 2:ncol(plotdata.vcm3[[i]]))
      lines(plotdata.vcm3[[i]][ind,1], plotdata.vcm3[[i]][ind,j], lwd=1.5)
      
    pr <- range(plotdata.vcm3[[length(varnames)+i]][,-1])
    plot(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,2], ylim=pr,
         axes=TRUE, type="l", xlab="", ylab=varnames[i], main="VCM, 0.5, female",
         border=gray(0.7), cex.lab=1.4)
    for(j in 2:ncol(plotdata.vcm3[[length(varnames)+i]]))
      plot(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,j], add=TRUE, border=gray(0.7))
    for(j in 2:ncol(plotdata.vcm3[[length(varnames)+i]]))
      lines(plotdata.vcm3[[length(varnames)+i]][ind,1], plotdata.vcm3[[length(varnames)+i]][ind,j], lwd=1.5)
    }   
  }

dev.off()