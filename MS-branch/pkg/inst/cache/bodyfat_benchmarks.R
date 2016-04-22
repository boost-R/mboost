
library("mboost")

data("bodyfat", package = "mboost")
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4 - 1

### center independent variables
indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
bodyfat[indep] <- lapply(bodyfat[indep], function(x) x - mean(x))

### setup benchmark experiments with B = 100 bootstrap samples
n <- nrow(bodyfat)
set.seed(290875)
bs <- rmultinom(100, n, rep(1, n)/n)
obj <- boost_dpp(bffm, data = bodyfat)

boob <- c()
grid <- seq(from = 2, to = 100, by = 2)

for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f)
        GaussReg()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- glmboost_fit(obj, family = GaussReg(),
        control = boost_control(mstop = max(grid)), weights = b)
    gm <- mboost:::gm.glmboost(object)
    fm <- t(apply(gm, 1, cumsum))[,grid]
    boob <- rbind(boob, apply(fm, 2, function(f) foo(f)))
}


bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4

boobfull <- c()

for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f)
        GaussReg()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- lm(bffm, data = bodyfat, weights = b)
    boobfull <- c(boob, foo(predict(object, newdata = bodyfat)))
}


bffm <- DEXfat ~ hipcirc + kneebreadth + anthro3a

boobrest <- c()

for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f)
        GaussReg()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- lm(bffm, data = bodyfat, weights = b)
    boobrest <- c(boob, foo(predict(object, newdata = bodyfat)))
}

save(boob, boobfull, boobrest, grid, file = "bodyfat_benchmarks.rda")
