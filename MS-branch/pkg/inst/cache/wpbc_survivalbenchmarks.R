library("mboost")

### Wisconsin prognostic breast cancer data, without missing values
data("wpbc")
wpbc <- wpbc[complete.cases(wpbc),]
iw <- IPCweights(Surv(wpbc$time, wpbc$status == "R"))
wpbc3 <- wpbc[,colnames(wpbc) != "status"]

wpbc3 <- wpbc3[iw > 0,]
iw <- iw[iw > 0]

### setup benchmark experiments with B = 100 bootstrap samples
n <- nrow(wpbc3)
set.seed(290875)
bs <- rmultinom(10, n, rep(1, n)/n)
obj <- boost_dpp(log(time) ~ ., data = wpbc3)

boob <- c()
grid <- seq(from = 5, to = 1000, by = 10)

### weighted squared error 
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f)
        GaussReg()@risk(obj$yfit, f, as.numeric(b == 0) * iw) / sum((b == 0) * iw)
    object <- glmboost_fit(obj, control = boost_control(mstop = max(grid)), weights = b * iw)

    gm <- mboost:::gm.glmboost(object)
    fm <- t(apply(gm, 1, cumsum))[,grid]
    boob <- rbind(boob, apply(fm, 2, function(f) foo(f)))
}

save(boob, grid, file = "wpbc_survivalbenchmarks.rda")
