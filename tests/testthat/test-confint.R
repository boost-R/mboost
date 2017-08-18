require("mboost")

set.seed(1907)

### check confidence intervals
data("bodyfat", package = "TH.data")
bodyfat$ID <- factor(sample(1:5, size = nrow(bodyfat), replace = TRUE))
glm <- glmboost(DEXfat ~ ., data = bodyfat)
gam <- gamboost(DEXfat ~ ., data = bodyfat)

refit <- glm$update(weights = model.weights(glm), risk = "inbag")
stopifnot(all.equal(coef(refit), coef(glm)))

glm[200]
confint.glm <- confint(glm, B = 100, B.mstop = 2)
confint.glm
print(confint.glm, which = 2)
res <- try(print(confint.glm, which = 20), silent = TRUE)
stopifnot(inherits(res, "try-error"))
print(confint.glm, level = 0.8, pe = TRUE)

confint.gam <- confint(gam, B = 100, B.mstop = 1)
plot(confint.gam, which = 1)
plot(confint.gam, which = 2)
plot(confint.gam, which = 3)


### check cvrisk (it should run even if a fold leads to an error)
folds <- cv(model.weights(glm), type = "kfold")

folds[1, 1] <- NA
cvrisk(glm, folds = folds, papply = lapply)
cvrisk(glm, folds = folds, papply = mclapply)

## test if cvrisk starts at 0 and provides a sensible model

data <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
glm <- glmboost(y ~ ., data = data)
gam <- gamboost(y ~ ., data = data)

cvr.glm <- cvrisk(glm)
cvr.gam <- cvrisk(gam)
stopifnot(mstop(cvr.glm) == 0)
stopifnot(mstop(cvr.gam) == 0)
