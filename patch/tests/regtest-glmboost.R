
require("mboost")

set.seed(290875)

dgp <- function(n = 100, beta = rep(0, 10), sd = 1) {
    p <- length(beta) - 1
    x <- cbind(1, matrix(runif(n * p), ncol = p))
    lp <- x %*% beta
    y <- lp + rnorm(n, sd = sd)
    ls <- data.frame(y = y, x[,-1])
    attr(ls, "lp") <- lp
    ls
}

### well-defined problem
mydf <- dgp(beta = c(1, 2.5, rep(0, 2)))

### for easy comparison with lm
fm <- Gaussian()
fm@offset <- function(y, w) 0

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm,
                    control = boost_control(mstop = 1000, nu = 1))
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
stopifnot(max(abs(coef(mydf.gb) - coef(mydf.lm))) < 1e-10)

### a little bit more difficult
mydf <- dgp(beta = c(1, 2.5, rep(0, 38)))

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm,
                    control = boost_control(mstop = 1000, nu = 1))
aic <- AIC(mydf.gb, method = "corrected")
ht <- hatvalues(mydf.gb)
mstop(aic)
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
which(abs(coef(mydf.lm)) < abs(coef(mydf.gb[mstop(aic)])))

#### check boosting hat matrix and subsetting / predict
stopifnot(isTRUE(all.equal(drop(attr(ht, "hatmatrix") %*% mydf$y),
                           as.vector(predict(mydf.gb[1000])))))
ht25 <- hatvalues(mydf.gb[25])
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% mydf$y),
                           as.vector(predict(mydf.gb[25])))))
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% mydf$y),
                           as.vector(fitted(mydf.gb[25])))))

### a simple two-dimensional example from `glmboost.Rd'
data("cars")
cars.gb <- glmboost(dist ~ speed, data = cars, family = fm,
                    control = boost_control(mstop = 1000, nu = 1))
cars.gb

### coefficients should coincide
cf <- coef(cars.gb)
attr(cf, "offset") <- NULL
stopifnot(all.equal(cf, coef(lm(dist ~ speed, data = cars))))

### logistic regression
mydf <- data.frame(x = runif(100), z = rnorm(100),
                   y = factor(c(rep(0, 30), rep(1, 70))))
bmod <- glmboost(y ~ x + z, data = mydf, family = Binomial(),
                 control = boost_control(mstop = 1000, nu = 1))
gmod <- glm(y ~ x + z, data = mydf, family = binomial())
llg <- logLik(gmod)
attributes(llg) <- NULL
stopifnot(all.equal(logLik(bmod), llg))
stopifnot(max(abs(predict(gmod, type = "link")/2 - fitted(bmod))) <
                  sqrt(.Machine$double.eps))
cfb <- coef(bmod, off2int = TRUE) * 2
stopifnot(all.equal(cfb, coef(gmod)))
aic <- AIC(bmod, "classical")
stopifnot(abs(AIC(gmod) - attr(aic, "AIC")[mstop(bmod)]) < 1e-5)

### weighted least squares problem

x <- runif(100)
df <- data.frame(y = 2 + 3 * x + rnorm(length(x)),
                 x = x, z = runif(length(x)),
                 w = runif(length(x)) * 10)

### linear model, classical fit
lmmod <- lm(y ~ x + z, data = df, weights = w)

### linear model, boosting fit
lmb <- glmboost(y ~ x + z, data = df, weights = df$w,
                control = boost_control(mstop = 5000, nu = 1))

### compare fitted values
stopifnot(max(abs(fitted(lmmod) -fitted(lmb))) < sqrt(.Machine$double.eps))

### compare hat matrices
stopifnot(max(abs(hatvalues(lmmod) - hatvalues(lmb))) < sqrt(.Machine$double.eps))

### compare boosting hat matrix with fitted values
stopifnot(max(abs(attr(hatvalues(lmb), "hatmatrix") %*% (df$y - lmb$offset) + lmb$offset -
        fitted(lmb))) < sqrt(.Machine$double.eps))

### Cox model (check for CoxPH family)
if (require("survival")) {

    test <- data.frame(time = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                       event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                       x     = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0))

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test, method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test,
                       family = CoxPH(),
                       control = boost_control(mstop = 2000, nu = 1)), which = 1:2)[2]))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))

    indx <- c(1, 1, 1, 2:10)
    w <- tabulate(indx)

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test, weights = w,
                                   method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test, weights = w,
                       family = CoxPH(),
                       control = boost_control(mstop = 200, nu = 1)), which = 1:2)[2]))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))

    indx <- c(1, 1, 1, 3:10)
    w <- tabulate(indx)

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test[indx,],
                                   method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test, weights = w,
                       family = CoxPH(),
                       control = boost_control(mstop = 1000)), which = 1:2)[2], tolerance = .Machine$double.eps ^ 0.125))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))
}


## Cox model with predictions obtained from survFit function

fm <- Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps - 1
fit <- coxph(fm, data = ovarian)
fit2 <- glmboost(fm, data = ovarian, family = CoxPH(),
    control=boost_control(mstop = 1000, center = TRUE))
fit3 <- glmboost(fm, data = ovarian, family = CoxPH(),
    control=boost_control(mstop = 1000, center = FALSE))

A1 <- survfit(fit, censor = FALSE)
A2 <- survFit(fit2)
A3 <- survFit(fit3)

max(A1$surv-A2$surv)
max(A1$surv-A3$surv)

newdata <- ovarian[c(1,3,12),]
A1 <- survfit(fit, newdata = newdata, censor = FALSE)
A2 <- survFit(fit2, newdata = newdata)
A3 <- survFit(fit3, newdata = newdata)

max(A1$surv-A2$surv)
max(A1$surv-A3$surv)

### Poisson models

df <- data.frame(x1 = runif(100), x2 = runif(100))
f <- -1 + 3 * df$x2
df$y <- round(exp(f) )
ctrl <- boost_control(mstop = 2000, nu = 0.1)

gmod <- glm(y ~ x1 + x2, data = df, family = poisson())
gbmod <- glmboost(y ~ x1 + x2, data = df, family = Poisson(), control = ctrl)

llg <- logLik(gmod)
attributes(llg) <- NULL
stopifnot(all.equal(logLik(gbmod), llg))

### hat matrix is only approximate!
stopifnot(abs(AIC(gmod) - attr(AIC(gbmod, "classical"), "AIC")[mstop(gbmod)]) < 1)

stopifnot(max(abs(predict(gmod) -  predict(gbmod))) < 1e-4)

### predictions:
set.seed(1907)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
y <- rnorm(100, mean = 3 * x1, sd = 2)
DF <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

amod <- glmboost(y ~ -1 + x1 + x2, data = DF)
agg <- c("none", "sum", "cumsum")
whi <- list(NULL, 1, 2, c(1,2))
for (i in 1:4){
    pred <- vector("list", length=3)
    for (j in 1:3){
        pred[[j]] <- predict(amod, aggregate=agg[j], which = whi[[i]])
    }
    if (i == 1){
        stopifnot(max(abs(pred[[2]] - pred[[3]][,ncol(pred[[3]])]))  < sqrt(.Machine$double.eps))
        if ((pred[[2]] - rowSums(pred[[1]]))[1] - attr(coef(amod), "offset") < sqrt(.Machine$double.eps))
            warning(sQuote("aggregate = sum"), " adds the offset, ", sQuote("aggregate = none"), " doesn't.")
        stopifnot(max(abs(pred[[2]] - rowSums(pred[[1]]) - attr(coef(amod), "offset")))   < sqrt(.Machine$double.eps))
    } else {
        stopifnot(max(abs(pred[[2]] - sapply(pred[[3]], function(obj) obj[,ncol(obj)])))  < sqrt(.Machine$double.eps))
        stopifnot(max(abs(pred[[2]] - sapply(pred[[1]], function(obj) rowSums(obj))))  < sqrt(.Machine$double.eps))
    }
}

agg <- c("none", "sum", "cumsum")
whi <- list(NULL, "x1", "x2", c("x1","x2"))
for (i in 1:4){
    pred <- vector("list", length=3)
    for (j in 1:3){
        pred[[j]] <- predict(amod, aggregate=agg[j], which = whi[[i]])
    }
    if (i == 1){
        stopifnot(max(abs(pred[[2]] - pred[[3]][,ncol(pred[[3]])]))  < sqrt(.Machine$double.eps))
        if ((pred[[2]] - rowSums(pred[[1]]))[1] - attr(coef(amod), "offset") < sqrt(.Machine$double.eps))
            warning(sQuote("aggregate = sum"), " adds the offset, ", sQuote("aggregate = none"), " doesn't.")
        stopifnot(max(abs(pred[[2]] - rowSums(pred[[1]]) - attr(coef(amod), "offset")))   < sqrt(.Machine$double.eps))
    } else {
        stopifnot(max(abs(pred[[2]] - sapply(pred[[3]], function(obj) obj[,ncol(obj)])))  < sqrt(.Machine$double.eps))
        stopifnot(max(abs(pred[[2]] - sapply(pred[[1]], function(obj) rowSums(obj))))  < sqrt(.Machine$double.eps))
    }
}

y <- rnorm(100, mean = 3 * x1^2, sd = 2)
DF2 <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
amod <- glmboost(y ~ -1 + x1 + I(x1^2), data = DF2)
stopifnot(ncol(predict(amod, which="x1")) == 2 && all(rowSums(predict(amod, which="x1")) + attr(coef(amod), "offset") - predict(amod) < sqrt(.Machine$double.eps)))


amod <- glmboost(y ~ 1+ x1 + x2, data = DF)
pr1 <- predict(amod, aggre = "sum", which= 1:2)
foo <- DF
foo$x2 <- 0
pr2 <- predict(amod, aggre = "sum", newdata=foo)
stopifnot(rowSums(pr1) + attr(coef(amod),"offset") - pr2 < sqrt(.Machine$double.eps))
newData <- as.data.frame(rbind(mean(DF)[-1], mean(DF)[-2]+1*sd(DF)[-1]))
if (!is.list(pr <- predict(amod, newdata=newData, which=1:2)))
    warning("predict(amod, newdata=newData, which=1:2) does not return a list") # no list but a matrix is returned!
stopifnot(is.list(pr <- predict(amod, newdata=newData, aggregate="cumsum", which=1:2)))
amod[10]
pr <- predict(amod, which=1:3)
stopifnot(ncol(pr) == 3 || all(pr[,c(1,ncol)] == 0))
amod[100]

# compare predictions with gamboost
mod1 <- glmboost(y ~ -1 + x1 + x2 + x3, data = DF)
mod2 <- gamboost(y ~ x1 + x2 + x3, data = DF, baselearner= function(x) bols(x, intercept=FALSE))
pr1_2 <- predict(mod2, aggre = "cumsum")
pr2_2 <- predict(mod2, aggre = "none")
pr3_2 <- predict(mod2, aggre = "sum")

stopifnot(max(abs(predict(mod1) - predict(mod2))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, aggre = "none") - predict(mod2, aggre = "none"))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, aggre = "cumsum") - predict(mod2, aggre = "cumsum"))) < sqrt(.Machine$double.eps))

# check type argument
set.seed(1907)
x1 <- rnorm(100)
p <- 1/(1 + exp(- 3 * x1))
y <- as.factor(runif(100) < p)
DF <- data.frame(y = y, x1 = x1)

logitBoost <- glmboost(y ~ x1, family = Binomial(),
                 data = DF,  control=boost_control(mstop=5000))
logit <- glm(y ~ x1, data=DF, family=binomial)
stopifnot(coef(logitBoost)[2]*2 - coef(logit)[2] < 1e-5) # * 2 as we use y = {-1, 1}

pr <- predict(logitBoost)
pr2 <- predict(logit)
stopifnot(pr * 2 - pr2 < 1e-5)  # * 2 as we use y = {-1, 1}

pr <- predict(logitBoost, type="class")
pr2 <- predict(logit, type="response") > 0.5
foo <- table(pr, pr2)
stopifnot(foo[1,2] + foo[2,1] == 0)

pr <- predict(logitBoost, type="response")
pr2 <- predict(logit, type="response")
stopifnot(pr - pr2 < sqrt(.Machine$double.eps))

### coefficients:
set.seed(1907)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
y <- rnorm(100, mean = 3 * x1, sd = 2)
DF <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
amod <- glmboost(y ~ x1 + x2 + x3, data = DF)

stopifnot(length(coef(amod)) == 4)
amod[10]
stopifnot(length(coef(amod)) == 1)
stopifnot(length(coef(amod, which=1:3)) == 3)
