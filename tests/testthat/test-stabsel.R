library(mboost)
set.seed(1234)
data("bodyfat", package = "TH.data")

test_that("stabsel works", {
    ### low-dimensional example
    mod <- glmboost(DEXfat ~ ., data = bodyfat)
    
    ## check if parameters are ocrrectly pre-computed
    s1 <- stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                       sampling.type = "MB")
    s2 <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)
    s2b <- stabsel_parameters(mod, q = 3, PFER = 1, sampling.type = "MB")
    expect_equal(s1, s2)
    expect_equal(s1, s2b)
    
    s3 <- stabsel(mod, q = 3, PFER = 1, folds = NULL, 
                  sampling.type = "SS", eval = FALSE)
    expect_equal(s3$B, 50)
    
    ## dimension of folds does not match B
    expect_warning(s4 <- stabsel(mod, q = 3, PFER = 1, folds = subsample(model.weights(mod), B = 50),
                               B = 10, sampling.type = "SS", eval = FALSE),
                   ".*B should be equal to number of folds.*")
    expect_equal(s4$B, 50)
    
    ## check if computations work correctly
    s5 <- stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20)
    ## expect that s5 uses 40 subsamples
    expect_gt(length(setdiff(sort(unique(c(s5$phat))), (0:20)/20)), 0) 
    expect_equal(length(setdiff(sort(unique(c(s5$phat))), (0:40)/40)), 0)
    
    s6 <- stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", B = 20)
    ## expect that s6 uses 20 subsamples
    expect_equal(length(setdiff(sort(unique(c(s6$phat))), (0:20)/20)), 0) 
    expect_equal(length(setdiff(sort(unique(c(s6$phat))), (0:40)/40)), 0)
    
    ## test verbosity:
    expect_warning(stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20, grid = 0:3),
                   ".*mstop.*too small in .* of the 40 subsampling replicates.*")
    
    ## test if mstop / grid can be specified correctly
    expect_error(stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20, grid = 1:50),
                 ".*grid must be of the form")
    expect_error(stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20, grid = 0:50, mstop = 10),
                 ".*Please specify only one.*")
    set.seed(1234)
    s7 <- stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20, grid = 0:20, mstop = 20)
    set.seed(1234)
    s8 <- stabsel(mod, q = 3, PFER = 1, sampling.type = "SS", B = 20, mstop = 20)
    expect_equal(s7$max, s8$max)
})
