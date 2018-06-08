require("mboost")

### check confidence intervals
test_that("predict with mstop = 0 works", {
    m <- glmboost(speed ~ dist, data = cars, control = boost_control(mstop = 0))
    expect_warning(p <- predict(m, type = "response", newdata = cars), regexp  = NA)
    expect_named(p, as.character(1:50))
})
