library(mboost)
set.seed(1234)

testthat("rankMatrix does not give a warning", {
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
