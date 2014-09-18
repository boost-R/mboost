require("mboostDevel")
attach(asNamespace("mboostDevel"))

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

confint.gam <- confint(gam, B = 100, B.mstop = 1)
plot(confint.gam, which = 1)
plot(confint.gam, which = 2)
plot(confint.gam, which = 3)
