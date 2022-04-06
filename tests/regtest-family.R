
library("mboost")

set.seed(290875)
x <- rnorm(100)
y <- rnorm(100)
w <- drop(rmultinom(1, 100, rep(1 / 100, 100)))

G <- Gaussian()
fm <- Family(ngradient = G@ngradient, risk = G@risk)

glmboost(y ~ x, family = G)

glmboost(y ~ x, family = fm)

mG <- glmboost(y ~ x, family = Gaussian())
mH <- glmboost(y ~ x, family = Huber(10))
all.equal(coef(mG), coef(mH))


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

## Count models as before but with highly skewed data 
## (used to lead to an error as offset was not properly computed)
x <- rnorm(100)
y <- rnbinom(length(x), size = 2, mu = exp(x * 5))
mod <- glmboost(y ~ x, family = NBinomial())
## check same for Hurdle model
x <- x[y > 0]
y <- y[y > 0]
mod <- glmboost(y ~ x, family = Hurdle())

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
data("wpbc", package = "TH.data")
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
round(coef(modWeighted, which = "") - coef(modSubset, which = ""), 3)
}


## Binomial
y <- as.factor(sample(0:1, 100, replace = TRUE))
x1 <- rnorm(100)
x2 <- rnorm(100)

mod <- glmboost(y ~ x1 + x2, family = Binomial())
mod[500]
2 * coef(mod, off2int = TRUE)

glmMod <- glm(y ~ x1 + x2, family = 'binomial')
coef(glmMod)
stopifnot(all(abs((coef(glmMod) - coef(mod, off2int = TRUE) * 2)) < sqrt(.Machine$double.eps)))

## C-index boosting
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
  dat <- data.frame(time = stime, event = event, x1 = x1, x2 = x2)
  
  # compute ipcweights outside the family
  ipcw <- IPCweights(x = Surv(dat$time, dat$event))
  model1 <- glmboost(Surv(time,event)~x1+x2, family=Cindex(ipcw = ipcw),
                     control = boost_control(mstop=50), data = dat)
  
  
  # compute ipcweights inside
  model2 <- glmboost(Surv(time,event)~x1+x2, family=Cindex(ipcw = 1),
                     control = boost_control(mstop=50), data = dat)
  rbind(coef(model1, off2int = TRUE, which = ""), 
        coef(model2, off2int = TRUE, which = ""))
  
  stopifnot(identical(coef(model1), coef(model2)))
  
  # change sigma
  model1 <- glmboost(Surv(time,event)~x1+x2, family=Cindex(sigma = 0.01),
                     control = boost_control(mstop=20), data = dat)
  model2 <- glmboost(Surv(time,event)~x1+x2, family=Cindex(sigma = 0.2),
                     control = boost_control(mstop=20), data = dat)
  rbind(coef(model1, off2int = TRUE, which = ""), 
        coef(model2, off2int = TRUE, which = ""))
  stopifnot(!identical(coef(model1), coef(model2)))
  
}

## check Binomial(type = "glm")
glmModboost <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm"))
glmModboost[1000] 
round(rbind(coef(glmMod), coef(glmModboost, off2int =TRUE), 2*coef(mod, off2int = TRUE)),3)
## use different link
glmMod <- glm(y ~ x1 + x2, family = binomial(link = "probit"))
coef(glmMod)
glmModboost <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm", link = "probit"))
glmModboost[500]
round(rbind(coef(glmMod), coef(glmModboost, off2int =TRUE)),3)
## use matrix of successes and failures
y <- matrix(ncol = 2, nrow = length(x1), data = rpois(lambda = 30, n =2*length(x1)) )
glmMod <- glm(y ~ x1 + x2, family = binomial())
coef(glmMod)
glmModboost <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm"))
glmModboost[500]
round(rbind(coef(glmMod), coef(glmModboost, off2int =TRUE)),3)
## use binary vector
y <- rbinom(prob = plogis(x1 + x2), size = 1, n = length(x1))
glmMod <- glm(y ~ x1 + x2, family = binomial())
coef(glmMod)
glmModboost <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm"),
                        control = boost_control(nu = 0.2))
glmModboost[1000]
round(rbind(coef(glmMod), coef(glmModboost, off2int =TRUE)),2)


## 
## Binomial with other links
## and interface type = "adaboost" or "glm"
set.seed(123)
y <- as.factor(sample(0:1, 100, replace = TRUE))
x1 <- rnorm(100)
x2 <- rnorm(100)

mod <- glmboost(y ~ x1 + x2, family = Binomial(type = "adaboost", link = "cauchit"))
mod[500]
mod2 <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm", link = "cauchit"))
mod2[500]
glmMod <- glm(y ~ x1 + x2 , family = binomial(link = "cauchit"))
stopifnot(all(round(coef(glmMod),2) == round(coef(mod, off2int =TRUE),2)))
rbind(coef(glmMod), coef(mod, off2int = TRUE), coef(mod2, off2int = TRUE))

mod <- glmboost(y ~ x1 + x2, family = Binomial(type = "adaboost", link = "probit"))
mod[500]
glmMod <- glm(y ~ x1 + x2 , family = binomial(link = "probit"))
stopifnot(all(round(coef(glmMod),2) == round(coef(mod, off2int =TRUE),2)))
mod2 <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm", link = "probit"))
mod2[500]
rbind(coef(glmMod), coef(mod, off2int = TRUE), coef(mod2, off2int = TRUE))


mod <- glmboost(y ~ x1 + x2, family = Binomial(type = "adaboost", link = "log"))
mod[500]
glmMod <- glm(y ~ x1 + x2 , family = binomial(link = "log"))
stopifnot(all(round(coef(glmMod),2) == round(coef(mod, off2int =TRUE),2)))
mod2 <- glmboost(y ~ x1 + x2, family = Binomial(type = "glm", link = "log"))
mod2[500]
rbind(coef(glmMod), coef(mod, off2int = TRUE), coef(mod2, off2int = TRUE))

