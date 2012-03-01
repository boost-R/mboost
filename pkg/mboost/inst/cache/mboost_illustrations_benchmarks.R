
library("mboost")

### bodyfat benchmarks

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

mod1 <- glmboost(bffm, data = bodyfat)

grid <- seq(from = 2, to = 100, by = 2)
boob <- cvrisk(mod1, folds = bs, grid = grid)

mod2 <- glmboost(bffm, data = bodyfat, control = boost_control(nu = 0.2))
boobfull <- cvrisk(mod1, folds = bs, grid = 5000)

bffm <- DEXfat ~ hipcirc + kneebreadth + anthro3a
mod3 <- glmboost(bffm, data = bodyfat, control = boost_control(nu = 0.2)) 
boobrest <- cvrisk(mod3, folds = bs, grid = 5000)

save(boob, boobfull, boobrest, grid, file = "bodyfat_benchmarks.rda")

### curve estimation benchmarks

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
    mod <- mboost(y ~ bbs(x, df = 2.5), data = learn, 
                  control = boost_control(mstop = max(mstops)))
    
    fm <- predict(mod, agg = "cumsum")[,mstops]
    mseB[i,] <- apply(fm, 2, function(f) se(learn$ytrue, f)) / nrow(learn)
}

save(dfree, mstops, mseSS, mseB, file = "curve_estimation.rda")

### WPBC glm benchmarks

### Wisconsin prognostic breast cancer data, without missing values
data("wpbc")
wpbc <- wpbc[complete.cases(wpbc),]

indep <- names(wpbc)[!(names(wpbc) %in% c("time", "status"))]
wpbc[indep] <- lapply(wpbc[indep], function(x) x - mean(x))

### setup benchmark experiments with B = 100 bootstrap samples
n <- nrow(wpbc)
set.seed(290875)
bs <- rmultinom(100, n, rep(1, n)/n)

mod <- glmboost(status ~ ., data = wpbc[,-2], family = Binomial())

boob <- c()
grid <- seq(from = 5, to = 500, by = 5)

boob <- cvrisk(mod, folds = bs, grid = grid)

goob <- c()

flink <- function(x)
    pmin(abs(x), 18) * sign(x)

### negative binom. log-lik / n for glm(...)
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    y <- c(-1, 1)[as.integer(wpbc$status)]
    foo <- function(f) 
        Binomial()@risk(y, f, as.numeric(b == 0)) / sum(b == 0)
    object <- glm(status ~  ., data = wpbc[,-2], subset = b > 0, 
                  family = binomial(), weights = b)
    p <- predict(object, newdata = wpbc[,-2], type = "link")
    goob <- c(goob, foo(flink(p)))
}

soob <- c()

flink <- function(x)
    pmin(abs(x), 18) * sign(x)

### negative binom. log-lik / n for step(glm(...))
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f) 
        Binomial()@risk(y, f, as.numeric(b == 0)) / sum(b == 0)
    object <- step(glm(status ~  ., data = wpbc[,-2], subset = b > 0, 
                       family = binomial(), weights = b),
                       trace = 0)
    p <- predict(object, newdata = wpbc[,-2], type = "link")
    soob <- c(soob, foo(flink(p)))
}



save(bs, goob, soob, boob, file = "wpbc_benchmarks.rda")

### WPBC survival benchmarks

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

mod <- glmboost(log(time) ~ ., data = wpbc3)

grid <- seq(from = 5, to = 1000, by = 10)

boob <- cvrisk(mod, folds = bs, grid = grid)

save(boob, grid, file = "wpbc_survivalbenchmarks.rda")

