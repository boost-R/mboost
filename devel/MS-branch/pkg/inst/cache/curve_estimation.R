
library("mboost")

dgp <- function(n = 100) {
    x <- sort(runif(n) - 0.5)
    ytrue <- 0.8 * x + sin(6 * x)
    data.frame(x = x, y = ytrue + rnorm(n, sd = sqrt(2)), ytrue = ytrue)
}

nsim <- 100
dfree <- seq(from = 2, to = 40, by = 2)
mstops <- seq(from = 5, to = 1005, by = 10)

mseSS <- matrix(0, nrow = nsim, ncol = length(dfree))
mseB <- matrix(0, nrow = nsim, ncol = length(mstops))
se <- GaussReg()@risk

for (i in 1:nsim) {
    print(i)
    learn <- dgp()
    mseSS[i,] <- sapply(dfree, function(d) {
        se(learn$ytrue, predict(smooth.spline(x = learn$x, y = learn$y, df = d), 
                           x = learn$x)$y)/nrow(learn)
        })
    mod <- gamboost(y ~ x - 1, data = learn, dfbase = 2.5, 
                    control = boost_control(mstop = max(mstops)))
    
    gm <- mboost:::gm.gamboost(mod)
    fm <- t(apply(gm, 1, cumsum))[,mstops]
    mseB[i,] <- apply(fm, 2, function(f) se(learn$ytrue, f)) / nrow(learn)
}

save(dfree, mstops, mseSS, mseB, file = "curve_estimation.rda")
