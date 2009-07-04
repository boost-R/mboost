
library("splines")
library("Matrix")
library("mboost")

source("helpers.R")
source("bl.R")
source("mboost.R")


x <- gl(50, 1009) ###rpois(10, lambda = 10)
x[sample(1:length(x), 100)] <- NA
y <- rnorm(length(x))
x[sample(1:length(x), 100)] <- NA
w <- rpois(length(x), lambda = 1)

system.time(c1 <- bols3(x)$dpp(w)$fit(y)$model)
system.time(c2 <- coef(lm(y ~ x, weights = w)))
max(abs(c1 - c2))

set.seed(29)
x <- rpois(100, lambda = 7)
y <- rnorm(length(x))
w <- rpois(length(x), lambda = 1)
system.time(a1 <- bbs3(x)$dpp(w)$fit(y)$model)
system.time(a2 <- as.vector(attr(mboost:::bbs(x), "dpp")(w)$fit(y)$model))
max(abs(a1 - a2))

b1 <- bbs3(x)$dpp(w)
b1$predict(list(b1$fit(y), b1$fit(y + 2)), newdata = NULL, Sum = FALSE)
fitted(b1$fit(y))

a1 <- bbs3(x, y)$dpp(w)$fit(y)$model
a2 <- as.vector(attr(mboost:::bspatial(x, y, df = 4), "dpp")(w)$fit(y)$model)
max(abs(a1 - a2))

data("bodyfat", package = "mboost")

attach(bodyfat)
b <- list(blage = bbs3(age),
          blhih = bbs3(hipcirc))
a <- mboost_fit(b, DEXfat)

cc <- gamboost(DEXfat ~ age + hipcirc)

max(abs(a$predict() - cc$fit))

plot(model.frame(a, which = 2)[[c(1,1)]], predict(a, which = 2))

coef(a)

x <- names(bodyfat)
x <- x[x != "DEXfat"]
fm1 <- paste("DEXfat ~ ", paste("bbs3(", x, ")", collapse = "+"), sep = "")
fm1 <- as.formula(fm1)
system.time(a1 <- mboost(fm, data = bodyfat))

fm2 <- paste("DEXfat ~ ", paste("bbs(", x, ")", collapse = "+"), sep = "")
fm2 <- as.formula(fm2)
system.time(a2 <- gamboost(fm2, data = bodyfat))

(max(abs(drop(predict(a2)) - predict(a1))))
