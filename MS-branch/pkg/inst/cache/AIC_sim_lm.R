
library("mboost")

AIC_sim_lm <- function(dgp) {

    n <- 100
    m <- 251
    ctrl <- boost_control(mstop = m)

    res.df <- matrix(0, ncol = m, nrow = 100)
    res.pred <- array(0, c(100, n, m))
    res.obs <- matrix(0, ncol = n, nrow = 100)
    res.coef <- matrix(0, m, nrow = 100)
    res.stop <- rep(0, 100)

    for (i in 1:100) {
        print(i)
        learn <- dgp(n)
        mod <- glmboost(y ~ . - 1, data = learn, control = ctrl)
        ht <- hatvalues(mod)
        df <- attr(ht, "trace")
        res.stop[i] <- mstop(AIC(mod))
        res.df[i,] <- df + 1

        cf <- abs(mboost:::coefpath.glmboost(mod))
        res.coef[i,1:(m-1)] <- apply(cf[-1,], 1, function(x) sum(x > 0))

        lpm <- mboost:::gm.glmboost(mod)
        lp <- t(apply(lpm, 1, cumsum))
        res.pred[i,,] <- lp
        res.obs[i,] <- learn$y
    }

    cc.cov <- matrix(0, ncol = m, nrow = n)
    for (i in 1:n) {
        for (k in 1:m)
            cc.cov[i,k] <- cov(res.obs[,i],res.pred[,i,k])
    }

    df.true <- colSums(cc.cov)

    save(df.true, res.df, res.coef, 
         file = paste("AIC_sim_lm_", deparse(substitute(dgp)), ".Rda", sep = ""))
}



set.seed(22)
p <- 10
Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
x <- mvrnorm(n, rep(0, p), Sigma)
dgp1 <- function(n) {
    y <- sqrt(34.5)*x[,5]+ rnorm(n)
    data.frame(y = y, x)
}
AIC_sim_lm(dgp1)

set.seed(22)
p <- 200
Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
x <- mvrnorm(n, rep(0, p), Sigma)
dgp2 <- function(n) {
    y <- sqrt(34.5)*x[,5]+ rnorm(n)
    data.frame(y = y, x)
}
AIC_sim_lm(dgp2)

set.seed(22)
p <- 200
Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
x <- mvrnorm(n, rep(0, p), Sigma)
dgp3 <- function(n) {
    y <- 0.5*(2*apply(x[,1:10],1,sum) + apply(x[,11:20],1,sum)) + rnorm(n)
    data.frame(y = y, x)
}
AIC_sim_lm(dgp3)
