extractfit <- function(inb, fm, its)
{

blsstr <- labels(terms(fm))
attach(india)
bls <- list()
for(i in 1:length(blsstr))
  bls[[i]] <- eval(parse(text=blsstr[i]))
detach(india)

res <- list()
res$blsstr <- blsstr
res$its <- its
res$coefs <- list()
res$fits <- list()
res$ensemble <- list()
res$meanfits <- list()
res$meancoefs <- list()

for(i in 1:length(its))
 {
 coefs <- list()
 for(k in 1:length(bls))
   coefs[[k]] <- rep(0, dim(bls[[k]])[[2]])

 for(k in 1:its[i])
  {
  coefs[[inb$ensemble[k]]] <- coefs[[inb$ensemble[k]]] + inb$coef[[k]]
  }

 fits <- list()
 for(k in 1:length(bls))
   fits[[k]] <- bls[[k]] %*% coefs[[k]]

 meanfits <- list()
 for(k in 1:length(bls))
   {
   meanfits[[k]] <- fits[[k]] + inb$offset
   for(l in (1:length(bls))[-k])
     {
     if(!is.factor(india[,get("xname", environment(attr(bls[[l]], "dpp")))]))
       {
       desX <- get("newX", environment(attr(bls[[l]], "dpp")))
       meanx <- mean(india[,get("xname", environment(attr(bls[[l]], "dpp")))])
       meanfits[[k]] <- meanfits[[k]] + (desX(meanx)[1,,drop=FALSE] %*% coefs[[l]])[1,1]
       }
     else
       meanfits[[k]] <- meanfits[[k]] + coefs[[l]][1]
     }
   }
   
 meancoefs <- list()
 for(k in 1:length(bls))
   {
   meancoefs[[k]] <- coefs[[k]] + inb$offset
   for(l in (1:length(bls))[-k])
     {
     if(!is.factor(india[,get("xname", environment(attr(bls[[l]], "dpp")))]))
       {
       desX <- get("newX", environment(attr(bls[[l]], "dpp")))
       meanx <- mean(india[,get("xname", environment(attr(bls[[l]], "dpp")))])
       meancoefs[[k]] <- meancoefs[[k]] + (desX(meanx)[1,,drop=FALSE] %*% coefs[[l]])[1,1]
       }
     else
       meancoefs[[k]] <- meancoefs[[k]] + coefs[[l]][1]
     }
   }
   
 res$ensemble[[i]] <- table(factor(inb$ensemble[1:its[i]], levels=1:length(blsstr), labels=blsstr))
 res$coefs[[i]] <- coefs
 res$fits[[i]] <- fits
 res$meanfits[[i]] <- meanfits
 res$meancoefs[[i]] <- meancoefs
 }

return(res)
}

plotnonpar <- function(name="cage", index=1, label="age of the child", fit)
{
nits <- length(fit$its)
greyseq <- seq(0.9, 0, length=nits)

ind <- (1:length(india[,name]))[!duplicated(india[,name])]
o <- order(india[ind,name])
ind <- ind[o]
india[ind,name]

plot(india[ind,name], fit$fits[[nits]][[index]][ind], type="l", xlab=label, ylab="")
for(i in 1:nits)
  lines(india[ind,name], fit$fits[[i]][[index]][ind], col=grey(greyseq[i]))

return(invisible())
}

plotnonpar2 <- function(name="cage", index=1, label="age of the child", fits)
{
nits <- length(fits[[1]]$its)
ind <- (1:length(india[,name]))[!duplicated(india[,name])]
o <- order(india[ind,name])
ind <- ind[o]

xplot <- india[ind,name]
fitplot <- matrix(0, length(fits), length(xplot))
for(i in 1:length(fits))
  fitplot[i,] <-  fits[[i]]$fits[[nits]][[index]][ind]
rfit <- range(fitplot)

plot(xplot, fitplot[1,], type="l", xlab=label, ylab="", ylim=c(rfit[1],rfit[2]))
for(i in 2:length(fits))
  lines(xplot, fitplot[i,], lty=i)

return(invisible())
}

plotnonpar3 <- function(name="cage", index=1, label="age of the child", fits)
{
nits <- length(fits[[1]]$its)
ind <- (1:length(india[,name]))[!duplicated(india[,name])]
o <- order(india[ind,name])
ind <- ind[o]

xplot <- india[ind,name]
fitplot <- matrix(0, length(fits), length(xplot))
for(i in 1:length(fits))
  fitplot[i,] <-  fits[[i]]$meanfits[[nits]][[index]][ind]
rfit <- range(fitplot)

plot(xplot, fitplot[1,], type="l", xlab=label, ylab="", ylim=c(rfit[1],rfit[2]))
for(i in 2:length(fits))
  lines(xplot, fitplot[i,], lty=i)

return(invisible())
}


plotpar <- function(name = "csex", index = 3, label="sex of the child", fit)
{
labels <- levels(india[,name])[-1]

nits <- length(fit$its)
mat <- matrix(0,nits,length(fit$coefs[[nits]][[index]])-1)
for(i in 1:nits)
  mat[i,] <- fit$coefs[[i]][[index]][-1]

plot(fit$its, mat[,1], ylim=range(mat), type="l", xlab=label, ylab="")

if(dim(mat)[2]>1)
  {
  for(i in 2:(dim(mat)[2]))
    {
    lines(fit$its, mat[,i])
    }
  }
axis(4, at=mat[nits,], labels)

return(invisible())
}

plotpar2 <- function(name = "csex", index = 3, label="sex of the child", fits)
{
labels <- factor(levels(india[,name])[-1])

nits <- length(fits[[1]]$its)
mat <- matrix(0,length(fits),length(fits[[1]]$coefs[[nits]][[index]])-1)
for(i in 1:length(fits))
  mat[i,] <- fits[[i]]$coefs[[nits]][[index]][-1]
                        
plot(labels, mat[1,], ylim=range(mat), ylab=label)
for(i in 2:length(fits))
  plot(labels, mat[i,], lty=i, add=TRUE)
return(invisible())
}

plotpar3 <- function(name = "csex", index = 3, label="sex of the child", fits)
{
labels <- factor(levels(india[,name]))

nits <- length(fits[[1]]$its)
mat <- matrix(0,length(fits),length(fits[[1]]$meancoefs[[nits]][[index]]))
for(i in 1:length(fits))
  mat[i,] <- fits[[i]]$meancoefs[[nits]][[index]]
                        
plot(labels, mat[1,], ylim=range(mat), ylab=label)
for(i in 2:length(fits))
  plot(labels, mat[i,], lty=i, add=TRUE)
return(invisible())
}

