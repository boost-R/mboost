###################################################
### chunk number 1: pkg-attach
###################################################
source("setup.R")
library("Matrix")
library("mboost")
attach(asNamespace("mboost"))

f <- list.files(path = "../R", pattern = "R$", full = TRUE)
sapply(f, source)


###################################################
### chunk number 2: bodyfat-lm-fit
###################################################
bf_lm <- lm(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(bf_lm)  


###################################################
### chunk number 3: bodyfat-glmboost-fit
###################################################
bf_glm <- glmboost(DEXfat ~ ., data = bodyfat, control = boost_control(center = TRUE))
bf_glm3 <- Glmboost(DEXfat ~ ., data = bodyfat, control = boost_control(center = TRUE))

max(abs(attr(hatvalues(bf_glm), "trace") - attr(hatvalues(bf_glm3), "trace")))
AIC(bf_glm)
AIC(bf_glm3)

cv <- Cvrisk(bf_glm3)

###################################################
### chunk number 4: bodyfat-glmboost-coef
###################################################
max(abs(coef(bf_glm) - coef(bf_glm3)))

max(abs(predict(bf_glm, newdata = bodyfat[1:10]) - 
        predict(bf_glm3, newdata = bodyfat[1:10])))

###################################################
### chunk number 9: bodyfat-pkg-attach
###################################################
source("setup.R")


###################################################
### chunk number 10: bodyfat-gamboost-fit
###################################################
bf_gam <- gamboost(DEXfat ~ ., data = bodyfat)
bf_gam3 <- Gamboost(DEXfat ~ ., data = bodyfat)
sapply(1:length(coef(bf_gam)), function(i) max(abs(coef(bf_gam)[[i]] - coef(bf_gam3)[[i]])))

max(abs(attr(hatvalues(bf_gam), "trace") - attr(hatvalues(bf_gam3), "trace")))
AIC(bf_gam)
AIC(bf_gam3)


max(abs(predict(bf_gam, newdata = bodyfat[1:10]) -
        predict(bf_gam3, newdata = bodyfat[1:10])))


###################################################
### chunk number 13: bodyfat-pkg-attach
###################################################
source("setup.R")
library("splines")
indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
bsfm <- as.formula(paste("DEXfat ~ ", paste("bs(", indep, ")", collapse = " + ", sep = ""), sep = ""))


###################################################
### chunk number 14: bodyfat-bs
###################################################
bsfm


###################################################
### chunk number 15: bodyfat-fpboost-fit
###################################################
ctrl <- boost_control(mstop = 5000)
bf_bs <- glmboost(bsfm, data = bodyfat, control = ctrl)
bf_bs3 <- Glmboost(bsfm, data = bodyfat, control = ctrl)
max(abs(coef(bf_bs) - coef(bf_bs3)))

max(abs(predict(bf_bs, newdata = bodyfat[1:10]) -
        predict(bf_bs3, newdata = bodyfat[1:10])))

max(abs(attr(hatvalues(bf_bs), "trace") - attr(hatvalues(bf_bs3), "trace")))
AIC(bf_bs)
AIC(bf_bs3)


###################################################
### chunk number 17: pkg-attach
###################################################
source("setup.R")
n <- sum(complete.cases(wpbc))
p <- ncol(wpbc) - 2


###################################################
### chunk number 18: wpbc-glm-fit
###################################################
### remove missing values and time variable
cc <- complete.cases(wpbc)
wpbc2 <- wpbc[cc, colnames(wpbc) != "time"]
### fit logistic regression model
wpbc_step <- step(glm(status ~ ., data = wpbc2, family = binomial()), trace = 0)


###################################################
### chunk number 19: wpbc-glm-aic
###################################################
logLik(wpbc_step)
AIC(wpbc_step)


###################################################
### chunk number 20: wpbc-glmboost-fit
###################################################
### fit logistic regression model via gradient boosting
ctrl <- boost_control(mstop = 500, center = TRUE)
wpbc_glm <- glmboost(status ~ ., data = wpbc2, family = Binomial(),
                     control = ctrl)
wpbc_glm3 <- Glmboost(status ~ ., data = wpbc2, family = Binomial(),
                     control = ctrl)
max(abs(coef(wpbc_glm) - coef(wpbc_glm3)))

max(abs(predict(wpbc_glm, newdata = wpbc2[1:10,]) -
        predict(wpbc_glm3, newdata = wpbc2[1:10,])))


###################################################
### chunk number 23: wpbc-gamboost-fit
###################################################
wpbc_gam <- gamboost(status ~ ., data = wpbc2, family = Binomial())
wpbc_gam3 <- Gamboost(status ~ ., data = wpbc2, family = Binomial())
sapply(1:length(coef(wpbc_gam)), function(i) max(abs(coef(wpbc_gam)[[i]] - coef(wpbc_gam3)[[i]])))

max(abs(predict(wpbc_gam, newdata = wpbc2[1:10,]) -
        predict(wpbc_gam3, newdata = wpbc2[1:10,])))


###################################################
### chunk number 25: pkg-attach
###################################################
source("setup.R")


###################################################
### chunk number 26: wpbc-glmboost-PIC
###################################################
### calculate IPC weights
censored <- wpbc$status == "R"
iw <- IPCweights(Surv(wpbc$time, censored))
wpbc3 <- wpbc[,names(wpbc) != "status"]


###################################################
### chunk number 27: wpbc-glmboost-censored-fit
###################################################
ctrl <- boost_control(mstop = 500, center = TRUE)
wpbc_surv <- glmboost(log(time) ~ ., data = wpbc3,
                  control = ctrl, weights = iw)
wpbc_surv3 <- Glmboost(log(time) ~ ., data = wpbc3,
                  control = ctrl, weights = iw)
max(abs(coef(wpbc_surv) - coef(wpbc_surv3)))

max(abs(predict(wpbc_surv, newdata = wpbc3[1:10,]) -
        predict(wpbc_surv3, newdata = wpbc3[1:10,])))


###################################################
### chunk number 28: wpbc-glmboost-coef
###################################################
names(coef(wpbc_surv)[abs(coef(wpbc_surv)) > 0])

