require("mboost")

### check confidence intervals
test_that("predict with mstop = 0 works", {
    glm <- glmboost(speed ~ dist, data = cars, control = boost_control(mstop = 0))
    expect_silent(p <- predict(m, type = "response", newdata = cars))
    expect_equal(names(p), as.character(1:50))
})
