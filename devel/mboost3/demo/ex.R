
library("splines")
library("Matrix")
library("mboost")

f <- list.files(path = "../R", pattern = "R$", full = TRUE)
sapply(f, source)
library("MASS")
library("Matrix")

x <- gl(50, 1009) ###rpois(10, lambda = 10)
x[sample(1:length(x), 100)] <- NA
y <- rnorm(length(x))
x[sample(1:length(x), 100)] <- NA
w <- rpois(length(x), lambda = 1)

system.time(c1 <- bols3(x)$dpp(w)$fit(y)$model)
system.time(c2 <- coef(lm(y ~ x - 1, weights = w)))
max(abs(c1 - c2))

set.seed(29)
x <- rpois(100, lambda = 7)
y <- rnorm(length(x))
w <- rpois(length(x), lambda = 1)
system.time(a1 <- bbs3(x)$dpp(w)$fit(y)$model)
system.time(a2 <- as.vector(attr(mboost:::bbs(x), "dpp")(w)$fit(y)$model))
max(abs(a1 - a2))

b1 <- bbs3(x)$dpp(w)
b1$predict(list(b1$fit(y), b1$fit(y + 2)), newdata = NULL, aggre = "none")
fitted(b1$fit(y))

a1 <- bbs3(x, y)$dpp(w)$fit(y)$model
a2 <- as.vector(attr(mboost:::bspatial(x, y, df = 4), "dpp")(w)$fit(y)$model)
max(abs(a1 - a2))

data("bodyfat", package = "mboost")

ctrl <- boost_control(mstop = 500)
x <- names(bodyfat)
x <- x[x != "DEXfat"]
fm1 <- paste("DEXfat ~ ", paste("bbs3(", x, ")", collapse = "+"), sep = "")
fm1 <- as.formula(fm1)
system.time(a1 <- mboost(fm1, data = bodyfat))

fm2 <- paste("DEXfat ~ ", paste("bbs(", x, ")", collapse = "+"), sep = "")
fm2 <- as.formula(fm2)
system.time(a2 <- gamboost(fm2, data = bodyfat))

p1 <- predict(a1, newdata = bodyfat)
p2 <- predict(a2, newdata = bodyfat)

(max(abs(drop(p2) - p1)))


Rprof("a1")
a1 <- mboost(fm1, data = bodyfat, control = ctrl)
Rprof(NULL)
Rprof("a2")
a1[200]
Rprof(NULL)

f <- a1$basemodel[[3]]$fit
Rprof("a3")
for (i in 1:1000)
    z <- f(bodyfat$DEXfat)
Rprof(NULL)




x <- rnorm(10)
w <- rpois(length(x), lambda = 1)
y <- rnorm(10)

a <- fit(dpp(bolscw(x), w), y)

b1 <- mboost(y ~ bolscw(x), weights = w)
coef(b1)
predict(b1)

b2 <- glmboost(y ~ x, weights = w)
coef(b2)
predict(b2)




data("bodyfat", package = "mboost")

w <- rpois(nrow(bodyfat), lambda = 2)

system.time(b1 <- glmboost(DEXfat ~ ., data = bodyfat, weights = w))
system.time(b2 <- mboost(DEXfat ~ bolscw(bodyfat[, colnames(bodyfat) != "DEXfat"]), 
             data = bodyfat, weights = w))
b3 <- Glmboost(DEXfat ~ ., data = bodyfat, weights = w)



b1 <- glmboost(DEXfat ~ ., data = bodyfat, 
               control = boost_control(center = TRUE))


b3 <- Glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)

X <- bodyfat[,-2]

a <- dpp(bolscw(X, center = TRUE), rep(1, nrow(X)))$fit
a(bodyfat$DEXfat)


max(abs(coef(b1) - coef(b2)))
max(abs(predict(b1) - predict(b2)))
max(abs(b1$risk - b2$risk()))

Rprof("a5")
b1 <- glmboost(DEXfat ~ ., data = bodyfat, weights = w, control = boost_control(mstop = 1000))
Rprof(NULL)


Rprof("a4")
b2 <- mboost(DEXfat ~ bolscw(bodyfat[, colnames(bodyfat) != "DEXfat"]),
             data = bodyfat, weights = w, control = boost_control(mstop = 1000))
Rprof(NULL)


b1 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w)
x <- names(bodyfat)
x <- x[x != "DEXfat"]
fm2 <- paste("DEXfat ~ ", paste("bbs3(", x, ")", collapse = "+"), sep = "")
fm2 <- as.formula(fm2)
b2 <- mboost(fm2, data = bodyfat, weights = w)


sapply(1:length(coef(b1)), function(i) max(abs(coef(b1)[[i]] - coef(b2)[[i]])))
max(abs(predict(b1) - predict(b2)))
max(abs(b1$risk - b2$risk()))

Rprof("a6")
z <- mboost(DEXfat ~ btree(age) + btree(waistcirc), data = bodyfat)
Rprof(NULL)

library("gbm")

Rprof("a7")
z2 <- gbm(DEXfat ~ age + waistcirc, data = bodyfat, distr = "gaussian")
Rprof(NULL)


n <- 500
df <- data.frame(y = rnorm(n), x1 = round(runif(n), 2), 
                 x2 = round(runif(n), 2),
                 z1 = round(runif(n), 2), 
                 z2 = round(runif(n),2),
                 id = gl(100, n / 100))

system.time(a <- mboost(y ~ bbs3(x1) + bbs3(x2) + bspatial3(z1, z2, knots = 6) + 
                        brandom3(id), data = df))

system.time(b <- predict(a, components = TRUE))
system.time(b1 <- predict(a, newdata = df, components = TRUE))

