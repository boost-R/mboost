
library("mboostDevel")

set.seed(290875)
x <- rnorm(100)
y <- rnorm(100)
w <- drop(rmultinom(1, 100, rep(1 / 100, 100)))

G <- Gaussian()
fm <- Family(ngradient = G@ngradient, risk = G@risk)

glmboost(y ~ x, family = G)

glmboost(y ~ x, family = fm)


x <- rnorm(100)
y <- rnbinom(length(x), size = 2, mu = exp(x * 2))
mod <- glmboost(y ~ x, family = NBinomial(), center = FALSE)
mod[1000]
coef(mod)
nuisance(mod)

### QuantReg and ExpectReg
gamboost(y ~ x, family = QuantReg())
gamboost(y ~ x, family = ExpectReg())


if (require("MASS")) {

summary(glm.nb(y ~ x))

y <- cut(x, breaks = c(-Inf, quantile(x, prob = c(0.25, 0.5, 0.75)), Inf), ordered = TRUE)
x <- rnorm(100)
polr(y ~ x)
mod <- glmboost(y ~ x, family = PropOdds(), center = FALSE)
nuisance(mod) - attr(coef(mod), "offset")
coef(mod)

}


### Weibull model

if (require("survival")) {

rextrval <- function(x) log( -log(1-x) )
sigma <- 0.5
u <- runif(100)
u.1 <- runif(100)
w <- rextrval(u)
w.1 <- rextrval(u.1)

x1 <- rnorm(100,sd=1)
x2 <- x1 + rnorm(100,sd=1)
x1.1 <- rnorm(100,sd=1)
x2.1 <- x1.1 + rnorm(100,sd=1)
X <- cbind(x1,x2)
X.1 <- cbind(x1.1,x2.1)
beta <- c(1,0.5)
survtime <- exp(X%*%beta + sigma*w)
censtime <- exp(X.1%*%beta + sigma*w.1)
event <- survtime < censtime
stime <- pmin(survtime,censtime)

model1 <- glmboost(Surv(stime,event)~x1+x2, family=Weibull(),
    control = boost_control(mstop=100), center = FALSE)
coef(model1)
nuisance(model1)
model2 <- survreg(Surv(stime,event)~x1+x2)
coef(model2)
model2$scale

}


### Log logistic model

if (require("survival")) {

sigma <- 0.5
w <- rlogis(100)
w.1 <- rlogis(100)

x1 <- rnorm(100,sd=1)
x2 <- x1 + rnorm(100,sd=1)
x1.1 <- rnorm(100,sd=1)
x2.1 <- x1.1 + rnorm(100,sd=1)
X <- cbind(x1,x2)
X.1 <- cbind(x1.1,x2.1)
beta <- c(1,0.5)
survtime <- exp(X%*%beta + sigma*w)
censtime <- exp(X.1%*%beta + sigma*w.1)
event <- survtime < censtime
stime <- pmin(survtime,censtime)

model1 <- glmboost(Surv(stime,event)~x1+x2, family=Loglog(),
    control = boost_control(mstop=200), center = FALSE)
coef(model1)
nuisance(model1)
model2 <- survreg(Surv(stime,event)~x1+x2, dist="loglogistic")
coef(model2)
model2$scale

}


### Log normal model

if (require("survival")) {

sigma <- 0.5
w <- rnorm(100)
w.1 <- rnorm(100)

x1 <- rnorm(100,sd=1)
x2 <- x1 + rnorm(100,sd=1)
x1.1 <- rnorm(100,sd=1)
x2.1 <- x1.1 + rnorm(100,sd=1)
X <- cbind(x1,x2)
X.1 <- cbind(x1.1,x2.1)
beta <- c(1,0.5)
survtime <- exp(X%*%beta + sigma*w)
censtime <- exp(X.1%*%beta + sigma*w.1)
event <- survtime < censtime
stime <- pmin(survtime,censtime)

model1 <- glmboost(Surv(stime,event)~x1+x2, family=Lognormal(),
    control = boost_control(mstop=200), center = FALSE)
coef(model1)
nuisance(model1)
model2 <- survreg(Surv(stime,event)~x1+x2, dist="lognormal")
coef(model2)
model2$scale

}


### AUC
data("wpbc", package = "mboostDevel")
wpbc[,colnames(wpbc) != "status"] <- scale(wpbc[,colnames(wpbc) != "status"])
wpbc <- wpbc[complete.cases(wpbc), colnames(wpbc) != "time"]
mAUC <- gamboost(status ~ ., data = wpbc, family = AUC())
1 - mAUC$risk()


# rank-based boosting

if (require("survival")) {

set.seed(1907)
n <- 100
beta <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
p <- length(beta)
x <- matrix(rnorm(n*p), n, p)
yt <- x %*% beta + rnorm(n)
x <- cbind(1, x)
colnames(x) <- c("intercept", paste("x", 1:8, sep=""))
cens <- runif(n, 0, 6)
y <- exp(pmin(yt, cens))
del <- yt <= cens
## check fitting and cvrisk
mod <- glmboost(x = x, y = Surv(y, del),
                control = boost_control(mstop = 500, nu = 0.1),
                center = TRUE,
                family = Gehan())
coef(mod)
plot(mod$risk())
cvr <- cvrisk(mod, folds = cv(model.weights(mod), "kfold"), papply=lapply)
plot(cvr)
## check weighting:
wMat <- cv(rep(1, n), type = "kfold",B = 2)
modWeighted <- glmboost(x = x, y = Surv(y, del), weights = wMat[, 1],
                        control = boost_control(mstop = 300, nu = 0.20),
                        family = Gehan())
# same model with data set subseted:
modSubset <- glmboost(x = x[as.logical(wMat[, 1]),],
                      y = Surv(y, del)[as.logical(wMat[, 1]),],
                      control = boost_control(mstop = 300, nu = 0.20),
                      family = Gehan())
## <FIXME> there are still some minor discrepancies. Perhaps this is due to
## different pre-processing? </FIXME>
round(coef(modWeighted) - coef(modSubset), 3)
}
