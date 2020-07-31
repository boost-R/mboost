require("mboost")

set.seed(1907)
data("bodyfat", package = "TH.data")
bodyfat$ID <- factor(sample(1:5, size = nrow(bodyfat), replace = TRUE))

### fit models
glm <- glmboost(DEXfat ~ ., data = bodyfat)
gam <- gamboost(DEXfat ~ age + elbowbreadth + kneebreadth + anthro3a +
                    anthro3b + anthro3c + anthro4 + bbs(hipcirc, waistcirc, df = 6) + brandom(ID), data = bodyfat)

### check confidence intervals
test_that("confint.glmboost works", {
    
    refit <- glm$update(weights = model.weights(glm), risk = "inbag")
    ## must work for confint to work
    expect_equal(coef(refit), coef(glm))
    mstop(glm) <- 200
    confint.glm <- confint(glm, B = 100, B.mstop = 2)

    expect_output(print(confint.glm), 
                  ".*Bootstrap Confidence Intervals.*2.5%.*97.5%.*age.*-0.002.*0.05.*ID5.*")
    expect_output(print(confint.glm, which = 2), 
                  ".*Bootstrap Confidence Interval\n.*2.5%.*97.5%.*age.*-0.002.*0.05.*")
    expect_error(print(confint.glm, which = 20), "is wrongly specified")
    expect_output(print(confint.glm, level = 0.8, pe = TRUE), 
                  ".*Bootstrap Confidence Intervals.*beta.*10%.*90%.*age.*0.009.*0.0000.*0.035.*ID5.*")
})

test_that("confint.gamboost works", {
    
    expect_warning(confint.gam <- confint(gam, B = 10, B.mstop = 1), "zero weights", all = TRUE)
    expect_output(confint.gam2 <- confint(gam, B = 10, B.mstop = 0), "Start computing bootstrap confidence intervals")
    
    ## plot some effects
    plot(confint.gam, which = 1)
    lines(confint.gam, which = 1, level = 0.9)
    plot(confint.gam, which = 4)
    
    ## plot raw data
    plot(confint.gam, which = 5, raw = TRUE)
    lines(confint.gam, which = 5)
    lines(confint.gam, which = 5, raw = TRUE)
    
    ## level plots for interaction effects
    expect_warning(plot(confint.gam, which = 8), "The scale is not the same")
    ## return plots without printing
    res <- plot(confint.gam, which = 8, print_levelplot = FALSE)
    expect_equal(res$mean$main, "Mean surface")
    expect_equal(res$lowerCI$main, "2.5% CI surface")
    expect_equal(res$upperCI$main, "97.5% CI surface")
    
    ## plots for factors
    plot(confint.gam, which = 9)
    lines(confint.gam, which = 9, level = 0.8)

    ## B.mstop = 0 ad B.mstop = 1 almost identical in this CI
    plot(confint.gam, which = 4)
    lines(confint.gam2, which = 4)
    
    expect_error(plot(confint.gam), 
                 ".*Specify a single base-learner.*")
    expect_error(lines(confint.gam), 
                 ".*Specify a single base-learner.*")
    
})


