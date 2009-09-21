
require("mboost")

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

    pdiffs <- max(abs(predict(a$update(a$data, a$control, a$weights)) - predict(a)))
    stopifnot(pdiffs < sqrt(.Machine$double.eps))


    ### attach `gbm', quietly
    sink("tmpfile")
    require("gbm")
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
p1 <- predict(blackboost(x = x, y = y, family = Laplace()), newdata = x)
p2 <- predict(blackboost(medv ~ ., data = BostonHousing, family = Laplace()),
              newdata = BostonHousing)
stopifnot(identical(abs(max(p1 - p2)), 0))

## Cox model

fit2 <- blackboost(Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps, 
    data = ovarian, family = CoxPH(), control = boost_control(mstop = 1000, 
    center = TRUE))

A2 <- survFit(fit2)
A2

newdata <- ovarian[c(1,3,12),]
A2 <- survFit(fit2, newdata = newdata)
A2
