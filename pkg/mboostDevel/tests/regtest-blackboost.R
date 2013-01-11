
require("mboostDevel")
if (require("party")) {

set.seed(290875)

tst <- try(data("BostonHousing", package = "mlbench"))
if (!inherits(tst, "try-error")) {

    system.time(a <- blackboost(medv ~ ., data = BostonHousing,
        tree_controls = ctree_control(teststat = "max",
                                testtype = "Teststatistic",
                                mincriterion = 0,
                                maxdepth = 2),
        control = boost_control(mstop = 500)))

    print(ae <- mean((predict(a) - BostonHousing$medv)^2))

    pdiffs <- max(abs(predict(update(a, model.weights(a))) - predict(a)))
    stopifnot(pdiffs < sqrt(.Machine$double.eps))


    ### attach `gbm', quietly
    sink("tmpfile")
    if (require("gbm")) cat()
    sink()
    file.remove("tmpfile")

    if (require("gbm")) {
        system.time(b <- gbm(medv ~ ., data = BostonHousing,
            n.trees = 500, interaction = 2, distribution = "gaussian",
            shrinkage = 0.1, bag = 1))
    print(be <- mean((predict(b, newdata = BostonHousing, n.trees = 500) -
                BostonHousing$medv)^2))
    plot(BostonHousing$medv, predict(a), col = "red", pch = "+")
    points(BostonHousing$medv,
           predict(b, newdata = BostonHousing, n.trees = 500),
           col = "blue", pch = "+")
    stopifnot(ae < be)
    }
}

### check different interfaces
x <- as.matrix(BostonHousing[,colnames(BostonHousing) != "medv"])
y <- BostonHousing$medv
p2 <- predict(blackboost(medv ~ ., data = BostonHousing, family = Laplace()),
              newdata = BostonHousing)

## Cox model
library("survival")

fit2 <- blackboost(Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps,
    data = ovarian, family = CoxPH(), control = boost_control(mstop = 1000))

A2 <- survFit(fit2)
A2

newdata <- ovarian[c(1,3,12),]
A2 <- survFit(fit2, newdata = newdata)
A2

### predictions:
set.seed(1907)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
y <- rnorm(100, mean = 3 * x1, sd = 2)
DF <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

amod <- blackboost(y ~ -1 + x1 + x2, data = DF)
agg <- c("none", "sum", "cumsum")
whi <- list(NULL, 1)
for (i in 1:2){
    pred <- vector("list", length=3)
    for (j in 1:3){
        pred[[j]] <- predict(amod, aggregate=agg[j], which = whi[[i]])
    }
    if (i == 1){
        stopifnot(max(abs(pred[[2]] - pred[[3]][,ncol(pred[[3]])]))  < sqrt(.Machine$double.eps))
        if ((pred[[2]] - rowSums(pred[[1]]))[1] - amod$offset < sqrt(.Machine$double.eps))
            warning(sQuote("aggregate = sum"), " adds the offset, ", sQuote("aggregate = none"), " doesn't.")
        stopifnot(max(abs(pred[[2]] - rowSums(pred[[1]]) - amod$offset))   < sqrt(.Machine$double.eps))
    } else {
        stopifnot(max(abs(pred[[2]] - sapply(pred[[3]], function(obj) obj[,ncol(obj)])))  < sqrt(.Machine$double.eps))
        stopifnot(max(abs(pred[[2]] - sapply(pred[[1]], function(obj) rowSums(obj))))  < sqrt(.Machine$double.eps))
    }
}

stopifnot(all(predict(amod, which=1) + amod$offset  - predict(amod) < sqrt(.Machine$double.eps)))


# check type argument
set.seed(1907)
x1 <- rnorm(100)
p <- 1/(1 + exp(- 3 * x1))
y <- as.factor(runif(100) < p)
DF <- data.frame(y = y, x1 = x1)

mod <- blackboost(y ~ x1, family = Binomial(),
                  data = DF,  control=boost_control(mstop=5000))

pr <- predict(mod)
pr <- predict(mod, type="class")
foo <- table(pr, y)
stopifnot(foo[1,2] + foo[2,1] == 0)
pr <- predict(mod, type="response")
# <FIXME> How do we check "correctness" of results?</FIXME>

}
