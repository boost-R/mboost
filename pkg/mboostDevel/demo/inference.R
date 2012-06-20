
library("mboostDevel")

data("bodyfat", package = "mboostDevel")

object <- glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)
a <- stabsel(object, cutoff = 0.9)
a
plot(a)

object <- gamboost(DEXfat ~ ., data = bodyfat)
a <- stabsel(object, cutoff = 0.9)
a
plot(a)
