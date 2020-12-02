require("mboost")

set.seed(1907)

data <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
glm <- glmboost(y ~ ., data = data)
gam <- gamboost(y ~ ., data = data)
cvr.glm <- cvrisk(glm)
cvr.gam <- cvrisk(gam)

### check confidence intervals
test_that("cvrisk works with 'broken' folds", {
    folds <- cv(model.weights(glm), type = "kfold")
    folds[1, 1] <- NA
    expect_warning(cvrisk(glm, folds = folds, papply = lapply),
                   ".*1 fold.* encountered an error.*Results are based on 9 folds only.*")
    expect_warning(cvrisk(gam, folds = folds, papply = mclapply),
                   ".*1 fold.* encountered an error.*Results are based on 9 folds only.*")

})

test_that("cvrisk starts at 0 and provides a sensible model", {
    expect_equal(dim(cvr.glm), c(25,101))
    expect_equal(colnames(cvr.glm), as.character(0:100))
    
    expect_equal(dim(cvr.gam), c(25,101))
    expect_equal(colnames(cvr.gam), as.character(0:100))
    
    expect_equal(mstop(cvr.glm), 2)
    expect_equal(mstop(cvr.gam), 14)
})

test_that("print.cvrisk works", {
    expect_output(print(cvr.glm), ".*Cross-validated Squared Error.*glmboost.*Optimal number of boosting iterations.*")
    expect_output(print(cvr.gam), ".*Cross-validated Squared Error.*gamboost.*Optimal number of boosting iterations.*")    
})

test_that("sampling types work correctly", {

    set.seed(1234)
    folds_bootstrap <- cv(weights = rep(1, 100))
    folds_subsampling <- cv(weights = rep(1, 100), type = "subsampling")
    folds_kfold <- cv(weights = rep(1, 100), type = "kfold")
    
    ## draw 100 observations with replacement
    expect_equal(colSums(folds_bootstrap), rep(100, 25))    
    expect_gt(sum(folds_bootstrap > 1), 0) ## some weights must be > 1
    
    ## draw 50% of observations randomly without replacement
    expect_false(all(rowSums(folds_subsampling) == rep(9, 100)))
    expect_equal(colSums(folds_subsampling), rep(50, 25))
    expect_true(all(folds_subsampling %in% c(0,1)))
    
    ## leave each observation out in one of 10 folds
    expect_equal(rowSums(folds_kfold), rep(9, 100))
    expect_equal(colSums(folds_kfold), rep(90, 10))
    expect_true(all(folds_kfold %in% c(0,1)))
    
    ## check that folds = NULL works as well
    set.seed(1234)
    weights <- model.weights(glm)
    folds <- rmultinom(25, length(weights), weights/sum(weights))
    expect_equivalent(folds, folds_bootstrap)
    set.seed(1234)
    expect_equal(cvrisk(glm, folds = NULL), cvrisk(glm, folds = folds_bootstrap))
    
    ## check user defined folds
    attr(folds_bootstrap, "type") <- NULL
    expect_equal(attr(cvrisk(glm, folds = folds_bootstrap), "type"), "user-defined")
})

if (require("survival")) {
    data("ovarian", package = "survival")
    fm <- Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps
    fit <- glmboost(fm, data = ovarian, family = CoxPH())
    
    test_that("crossvalidation works for CoxPH models", {
        expect_error(cvrisk(fit, folds = cv(weights = model.weights(fit), type = "kfold", B = nrow(ovarian)),
                            grid = 0:10), "Leave-one-out cross-validation cannot be used with .*family = CoxPH().*")
        
        expect_silent(cvr_uncor <- cvrisk(fit, grid = seq(0, 10, by = 2)))
        expect_equal(dim(cvr_uncor), c(25, 6))
        expect_gt(mstop(cvr_uncor), 0)
    })
}


