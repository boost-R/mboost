
library("mboost")
library("glmnet")

url <- "http://www-stat.stanford.edu/~hastie/glmnet/glmnetData"
wd <- getwd()
setwd(system.file("demo", package = "mboost"))
if (!file.exists("Leukemia.RData"))
    download.file(paste(url, "Leukemia.RData", sep = "/"), "Leukemia.RData")
load("Leukemia.RData")

if (!file.exists("InternetAd.RData"))
    download.file(paste(url, "InternetAd.RData", sep = "/"), "InternetAd.RData")
load("InternetAd.RData")

if (!file.exists("NewsGroup.RData"))
    download.file(paste(url, "NewsGroup.RData", sep = "/"), "NewsGroup.RData")
load("NewsGroup.RData")

dim(Leukemia$x)
system.time(a <- glmboost(x = Leukemia$x, y = as.factor(Leukemia$y), family = Binomial()))
system.time(b <- glmnet(x = Leukemia$x, y = Leukemia$y, family = "binomial"))

dim(InternetAd$x)
system.time(a <- glmboost(x = InternetAd$x, y = as.factor(InternetAd$y), family = Binomial()))
system.time(b <- glmnet(x = InternetAd$x, y = InternetAd$y, family = "binomial"))

dim(NewsGroup$x)
system.time(a <- glmboost(x = NewsGroup$x, y = as.factor(NewsGroup$y), family = Binomial()))
system.time(b <- glmnet(x = NewsGroup$x, y = NewsGroup$y, family = "binomial"))

setwd(wd)
