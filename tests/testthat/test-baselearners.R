library(mboost)
set.seed(1234)

test_that("rankMatrix does not give a warning", {
    ## setting with 'p>n' should not give warning in rankMatrix()
    n <- 10
    x1 <- rnorm(n)
    y <- 2*x1
    expect_that(mod1 <- mboost(y ~ bbs(x1, df = 4)),
                not(gives_warning()))
    expect_equal(dim(extract(mod1, "design")[[1]]),
                 c(10, 24))
    ## for comparison: setting with 'p<n'
    n <- 30
    x1 <- rnorm(n)
    y <- 2*x1
    expect_that(mod2 <- mboost(y ~ bbs(x1, df = 4)),
                not(gives_warning()))
    expect_equal(dim(extract(mod2, "design")[[1]]),
                 c(30, 24))
})


test_that("extrapolation is possible for kronecker product, kronecker product and sum", {
    x1 <- 1:90
    x2 <- 1:60
    y <- rnorm(length(x1) * length(x2))
    mod1 <- mboost(y ~ bbs(x1, df = 3, knots = 10) %O%
                       bbs(x2, df = 3, knots = 10),
                   control = boost_control(nu = 0.25, mstop = 10))
    expect_warning(predict(mod1, newdata=list(x1=x1, x2=1:61)),
                   "Linear extrapolation used")
    x <- expand.grid(x1, x2)
    mod2 <- mboost(y ~ bbs(Var2, df = 3, knots = 10) %X%
                       bbs(Var1, df = 3, knots = 10), data = x,
                   control = boost_control(nu = 0.25, mstop = 10))
    expect_warning(predict(mod2, newdata = data.frame(Var1=1:61, Var2=1:61)),
                   "Linear extrapolation used")
    mod3 <- mboost(y ~ bbs(Var2, df = 3, knots = 10) %+%
                       bbs(Var1, df = 3, knots = 10), data = x,
                   control = boost_control(nu = 0.25, mstop = 10))
    expect_warning(predict(mod3, newdata = data.frame(Var1=1:61, Var2=1:61)),
                   "Linear extrapolation used")
})


test_that("center = TRUE/FALSE throws an error in bols", {
    data("bodyfat", package = "TH.data")
    
    expect_error(mboost(DEXfat ~ bols(waistcirc, center=TRUE), data = bodyfat), regexp = ".*deprecated.*")
    expect_error(mboost(DEXfat ~ bols(waistcirc, center=FALSE), data = bodyfat), regexp = ".*deprecated.*")
    expect_error(mboost(DEXfat ~ bols(waistcirc, center=T), data = bodyfat), regexp = ".*deprecated.*")
    expect_error(mboost(DEXfat ~ bols(waistcirc, center=F), data = bodyfat), regexp = ".*deprecated.*")
    
    ## not throw errors if variable is called center
    bodyfat$center <- rnorm(nrow(bodyfat))
    mod <- mboost(DEXfat ~ bols(waistcirc, center), data = bodyfat)
    expect_length(coef(mod)[[1]], 3)
    expect_equal(names(coef(mod)[[1]])[3], "center")
    
    ## not throw errors if logical variable is called center
    bodyfat$center <- sample(c(TRUE, FALSE), nrow(bodyfat), replace = TRUE)
    mod <- mboost(DEXfat ~ bols(waistcirc, center), data = bodyfat)
    expect_length(coef(mboost(DEXfat ~ bols(waistcirc, center), data = bodyfat))[[1]], 3)
    expect_equal(names(coef(mod)[[1]])[3], "centerTRUE")
})
