
library("splines")
library("Matrix")
library("mboost")

f <- list.files(path = "../R", pattern = "R$", full = TRUE)
sapply(f, source)


###################################################
### chunk number 1: pkg-attach
###################################################
source("setup.R")


###################################################
### chunk number 2: bodyfat-lm-fit
###################################################
bf_lm <- lm(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(bf_lm)  


###################################################
### chunk number 3: bodyfat-glmboost-fit
###################################################
bf_glm <- Glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)


###################################################
### chunk number 4: bodyfat-glmboost-coef
###################################################
coef(bf_glm)


bf_gam <- Gamboost(DEXfat ~ ., data = bodyfat)
coef(bf_gam)


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
bf_bs <- Glmboost(bsfm, data = bodyfat, control = ctrl)
coef(bf_bs)


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
ctrl <- boost_control(mstop = 500)
wpbc_glm <- Glmboost(status ~ ., data = wpbc2, family = Binomial(), center = TRUE,
                     control = ctrl)
coef(wpbc_glm)

###################################################
### chunk number 23: wpbc-gamboost-fit
###################################################
wpbc_gam <- Gamboost(status ~ ., data = wpbc2, family = Binomial())
coef(wpbc_gam)

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
ctrl <- boost_control(mstop = 500)
wpbc_surv <- Glmboost(log(time) ~ ., data = wpbc3,
                  control = ctrl, weights = iw, center = TRUE)
coef(wpbc_surv)
names(coef(wpbc_surv)[abs(coef(wpbc_surv)) > 0])   
