require("mboost")

set.seed(1907)
data("bodyfat", package = "TH.data")
bodyfat$ID <- factor(sample(1:5, size = nrow(bodyfat), replace = TRUE))

### fit models
glm <- glmboost(DEXfat ~ ., data = bodyfat)
gam <- gamboost(DEXfat ~ age + elbowbreadth + kneebreadth + anthro3a +
                    anthro3b + anthro3c + anthro4 + bbs(hipcirc, waistcirc, df = 6), data = bodyfat)

### check confidence intervals
test_that("cvrisk works with 'broken' folds", {
    folds <- cv(model.weights(glm), type = "kfold")
    folds[1, 1] <- NA
    expect_warning(cvrisk(glm, folds = folds, papply = lapply),
                   ".*1 fold.* encountered an error.*Results are based on 9 folds only.*")
    expect_warning(cvrisk(glm, folds = folds, papply = mclapply),
                   ".*1 fold.* encountered an error.*Results are based on 9 folds only.*")

})

test_that("cvrisk starts at 0 and provides a sensible model", {
    
    data <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
    glm <- glmboost(y ~ ., data = data)
    gam <- gamboost(y ~ ., data = data)
    
    cvr.glm <- cvrisk(glm)
    cvr.gam <- cvrisk(gam)
    
    expect_equal(dim(cvr.glm), c(25,101))
    expect_equal(colnames(cvr.glm), as.character(0:100))
    
    expect_equal(dim(cvr.gam), c(25,101))
    expect_equal(colnames(cvr.gam), as.character(0:100))
    
    expect_equal(mstop(cvr.glm), 9)
    expect_equal(mstop(cvr.gam), 1)
})