library(mboost)

## library(mboost, lib.loc="~/R/experimental")
## attach(asNamespace("mboost"))
## library(Matrix)
## source("~/_Dateien_/R_Sources/mboost/pkg/R/brad.R")

########################################
#        Bivariate kriging             #
########################################

set.seed(1907)
x1 <- runif(300, -2, 2)
x2 <- runif(300, -x1, 2)
f <- function(x,y){
    x^2 + y^2
}
y <- rnorm(300, mean=f(x1,x2), sd = 0.1)
DF <- data.frame(y = y, x1=x1, x2=x2)
x1 <- x2 <- seq(-2,2, by=0.1)
X <- expand.grid(x1,x2)
nD <-  data.frame(x1=X[,1], x2=X[,2])
z <- matrix(f(X[,1], X[,2]), nrow=length(x1), ncol=length(x2))

mstop <- 1000 ## brad seems to "converge" much quicker in this situation
# mstop <- 10000 ## now bspatial looks ok, too.

###
# Fit model with brad
model_1 <- gamboost(y ~ brad(x1, x2, knots=100),
                data=DF, control=boost_control(mstop=mstop))
pr <- matrix(predict(model_1, nD), nrow=length(x1), ncol=length(x2))

par(mfrow=c(1,3))

contour(x1, x2, z)
contour(x1, x2, pr, add=TRUE, lty="dashed", col="red")
with(DF, points(x1, x2, cex=0.5))

###
# Fit model with bspatial
model_2 <- gamboost(y ~ bspatial(x1, x2, knots=10),
                 data=DF, control=boost_control(mstop=mstop))
pr2 <- matrix(predict(model_2, nD), nrow=length(x1), ncol=length(x2))

contour(x1, x2, z)
contour(x1, x2, pr2, add=TRUE, lty="dashed", col="red", levels=c(-10:10))
with(DF, points(x1, x2, cex=0.5))

###
# Fit model with Krig (package fields)
require(fields)
TMP <- mboost:::hyper_brad(DF[,2:3], vary= "", knots=100, cov.function = stationary.cov, args=list(Covariance="Matern", smoothness = 1.5, theta=NULL))
a <- Krig(x = as.matrix(DF[,2:3]), Y = DF[,1], df=20, knots = TMP$knots, cov.args=TMP$args)
pr3 <- matrix(predict(a, nD), nrow=length(x1), ncol=length(x2))

contour(x1, x2, z)
contour(x1, x2, pr3, add=TRUE, lty="dashed", col="red")
with(DF, points(x1, x2, cex=0.5))

########################################
#        Univariate kriging            #
########################################

set.seed(1907)
x1 <- rnorm(300)
y <- rnorm(300, mean=x1^2, sd = 0.5)
DF <- data.frame(y = y, x1=x1, int=rep(1, length(x1)))

nknots <- 20
model_uni <- gamboost(y ~ bols(int, intercept=FALSE) + brad(x1, knots=nknots),
                data=DF, control=boost_control(mstop=mstop))

#cv_uni <- cvrisk(model_uni) ## no overfitting until mstop = 1000
#plot(cv_uni)

par(mfrow=c(1,1))
plot(model_uni, which=2)
lines(sort(x1), sort(x1)^2 - mean(sort(x1)^2), col="red")

########################################
#        Trivariate kriging            #
########################################

#.. todo

########################################
# Check effective range computation    #
########################################

x = DF[1:300,2:3]
require(fields)
theta_hat <- mboost:::effective_range(x, eps = 0.001, cov.function=stationary.cov,
                             args=list(Covariance = "Matern", smoothness = 1.5, theta = NULL))
stopifnot( (attr(theta_hat, "c_value") - 9.23339238) < sqrt(.Machine$double.eps) )
