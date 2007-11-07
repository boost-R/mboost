
library("mboost")

#for smoothing splines

fastlp <- function(object) {

    mstop <- nrow(object$ensemble)
    x <- object$data$input
    RET <- matrix(0, nrow = NROW(x[[1]]), ncol = mstop)
    nu <- object$control$nu
    f <- mboost:::fitted.baselist

    tmp <- object$offset
    for (m in 1:mstop) {
        RET[,m] <- tmp + nu * f(object$ensembless[[m]])
        tmp <- RET[,m]
    }
    return(RET)
}


AIC_sim_ss <- function(dgp) {

    m <- 200
    ctrl <- boost_control(mstop = m)

    res.df <- matrix(0, ncol = m, nrow = 100)
    res.pred <- array(0, c(100, n, m))
    res.obs <- matrix(0, ncol = n, nrow = 100)
    res.var <- matrix(0, ncol = m, nrow = 100)
    res.stop <- rep(0, m)

    for (i in 1:100) {
        print(i)
        learn <- dgp(n) 
        mod <- gamboost(y ~ . - 1, data = learn, control = ctrl)
        ht <- hatvalues(mod)
        df <- attr(ht, "trace")
        res.stop[i] <- mstop(AIC(mod))
        res.df[i,] <- df
        lp <- fastlp(mod)
        res.pred[i,,] <- lp
        for (k in 1:m)
            res.var[i,k] <- length(unique(mod[k]$xselect))
        res.obs[i,] <- learn$y
    }

    cc.cov <- matrix(0, ncol = m, nrow = n)
    for (i in 1:n) {
        for (k in 1:m)
            cc.cov[i,k] <- cov(res.obs[,i],res.pred[,i,k])
    }

    df.true <- colSums(cc.cov)

    save(df.true, res.df,
         file = paste("AIC_sim_ss_", deparse(substitute(dgp)), ".Rda", sep = ""))
}


set.seed(22)
n <- 50
p <- 20
x <- matrix(runif(n*p),ncol=p)
dgp1 <- function(n) {  
    y <- 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2 + 10*x[,4]+ 5*x[,5] + rnorm(n,0,1)
    data.frame(y = y, x)
}
AIC_sim_ss(dgp1)

set.seed(22)
x <- matrix(runif(n*p),ncol=p)
dgp2 <- function(n) {  
    y <- 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2 + 10*x[,4]+ 5*x[,5] + rnorm(n,0,sqrt(10))
    data.frame(y = y, x)
}
AIC_sim_ss(dgp2)
