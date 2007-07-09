
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
data("iris")
x <- as.matrix(iris[,colnames(iris) != "Species"])
y <- iris$Species
p1 <- predict(blackboost(x = x, y = y), newdata = x)
p2 <- predict(blackboost(Species ~ ., data = iris), newdata = iris)
stopifnot(identical(abs(max(p1 - p2)), 0))


