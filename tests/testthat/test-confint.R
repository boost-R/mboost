require("mboost")

set.seed(1907)
data("bodyfat", package = "TH.data")
bodyfat$ID <- factor(sample(1:5, size = nrow(bodyfat), replace = TRUE))

### fit models
glm <- glmboost(DEXfat ~ ., data = bodyfat)
gam <- gamboost(DEXfat ~ age + elbowbreadth + kneebreadth + anthro3a +
                    anthro3b + anthro3c + anthro4 + bbs(hipcirc, waistcirc, df = 6), data = bodyfat)

### check confidence intervals
test_that("confint.glmboost works", {
    
    refit <- glm$update(weights = model.weights(glm), risk = "inbag")
    ## must work for confint to work
    expect_equal(coef(refit), coef(glm))
    mstop(glm) <- 200
    confint.glm <- confint(glm, B = 100, B.mstop = 2)

    expect_output(print(confint.glm), 
                  ".*Bootstrap Confidence Intervals.*2.5%.*97.5%.*age.*-0.0100.*0.0478.*ID5.*")
    expect_output(print(confint.glm, which = 2), 
                  ".*Bootstrap Confidence Interval\n.*2.5%.*97.5%.*age.*-0.0100.*0.0478.*")
    expect_error(print(confint.glm, which = 20), "is wrongly specified")
    expect_output(print(confint.glm, level = 0.8, pe = TRUE), 
                  ".*Bootstrap Confidence Intervals.*beta.*10%.*90%.*age.*0.0077.*0.0000.*0.0354.*ID5.*")
})

test_that("confint.gamboost works", {
    
    confint.gam <- confint(gam, B = 100, B.mstop = 1)
    confint.gam2 <- confint(gam, B = 100, B.mstop = 0)
    
    ## no effect (but some variability)
    plot(confint.gam, which = 3)
    lines(confint.gam, which = 3, level = 0.9)
    ## effects
    plot(confint.gam, which = 4)
    plot(confint.gam, which = 5)
    ## level plots for interaction effects
    plot(confint.gam, which = 8)
    
    ## B.mstop = 0 ad B.mstop = 1 almost identical in this CI
    plot(confint.gam, which = 4)
    lines(confint.gam2, which = 4)
    
    expect_error(plot(confint.gam), 
                 ".*Specify a single base-learner.*")
})


