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
bf_glm <- glmboost(DEXfat ~ ., data = bodyfat, control = boost_control(center = TRUE))


###################################################
### chunk number 4: bodyfat-glmboost-coef
###################################################
coef(bf_glm)


###################################################
### chunk number 5: bodyfat-oob-plot
###################################################
load(system.file("cache/bodyfat_benchmarks.rda", package = "mboost"))
aic <- AIC(bf_glm)
pdf("figures/bodyfat_glmboost-bodyfat-oob-plot.pdf", version = "1.4", width = 6, height = 10)
par(mai = par("mai") * c(1, 1, 0.5, 1))
mopt <- grid[which.min(colMeans(boob))]
layout(matrix(1:2, nrow = 2))
perfplot(boob, grid, ylab = "Out-of-bootstrap squared error", 
    xlab = "Number of boosting iterations", alpha = 0.05)
abline(h = mean(boobrest), lty = 2)
lines(c(which.min(colMeans(boob)), which.min(colMeans(boob))), 
      c(0, min(colMeans(boob))), lty = 2)
points(which.min(colMeans(boob)), min(colMeans(boob)))
plot(aic, ylim = c(3, 5.5))
dev.off()


###################################################
### chunk number 6: bodyfat-glmboost-AIC
###################################################
mstop(aic <- AIC(bf_glm))


###################################################
### chunk number 7: bodyfat-glmboost-coef
###################################################
coef(bf_glm[mstop(aic)])


###################################################
### chunk number 8: bodyfat-glmboost-coef-count
###################################################
cf <- coef(bf_glm[mopt])
nsel <- sum(abs(cf) > 0)


###################################################
### chunk number 9: bodyfat-pkg-attach
###################################################
source("setup.R")


###################################################
### chunk number 10: bodyfat-gamboost-fit
###################################################
bf_gam <- gamboost(DEXfat ~ ., data = bodyfat, baselearner = "bss")


###################################################
### chunk number 11: bodyfat-gamboost-prune
###################################################
mstop(aic <- AIC(bf_gam))


###################################################
### chunk number 12: bodyfat-gamboost-plot
###################################################
bf_gam <- bf_gam[mstop(aic)]
fpartial <- mboost:::gamplot(bf_gam)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
par(mai = par("mai") * c(1, 1, 0.5, 1))
x <- bf_gam$data$input
varorder <- rev(order(colMeans(abs(fpartial))))[1:4]

out <- sapply(varorder, function(i) {
    plot(x[,i], fpartial[,i],  main = "",
         xlab = colnames(x)[i], ylab = expression(f[partial]),
         ylim = max(abs(fpartial))*c(-1, 1))
    abline(h = 0, lty = 2, lwd = 0.5)
    })


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
mstop(aic <- AIC(bf_bs))


###################################################
### chunk number 16: bodyfat-fpboost-plot
###################################################
layout(matrix(1:4, ncol = 2, byrow = TRUE))
par(mai = par("mai") * c(1, 1, 0.5, 1))
cf <- coef(bf_bs[mstop(aic)])
x <- bf_bs$data$x
varorder <- c("hipcirc", "waistcirc", "kneebreadth", "anthro3b")
fpartial <- sapply(varorder, function(v) {
    indx <- grep(v, names(cf))
    x[,indx] %*% cf[indx]
})
out <- sapply(varorder, function(i) {
    plot(bodyfat[,i], fpartial[,i],  main = "",
         xlab = i, ylab = expression(f[partial]), ylim = max(abs(fpartial)) * c(-1, 1))
    abline(h = 0, lty = 2, lwd = 0.5)
    })


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


###################################################
### chunk number 21: wpbc-glmboost-AIC
###################################################
aic <- AIC(wpbc_glm, "classical") 
aic


###################################################
### chunk number 22: wpbc-glmboost-fit2
###################################################
### fit with new mstop
wpbc_glm <- wpbc_glm[mstop(aic)]
coef(wpbc_glm)[abs(coef(wpbc_glm)) > 0]


###################################################
### chunk number 23: wpbc-gamboost-fit
###################################################
wpbc_gam <- gamboost(status ~ ., data = wpbc2, family = Binomial(), baselearner = "bss")
mopt <- mstop(aic <- AIC(wpbc_gam, "classical"))
aic


###################################################
### chunk number 24: wpbc-gamboost-plot
###################################################
fpartial <- mboost:::gamplot(wpbc_gam[mopt])
x <- wpbc_gam$data$input
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mai = par("mai") * c(1, 1, 0.5, 1))
out <- sapply(rev(order(colMeans(abs(fpartial))))[1:4], function(i) {
    plot(x[,i], fpartial[,i], xlab = colnames(x)[i], main = "",
         ylab = expression(f[partial]), ylim = c(-0.5, 0.5), type = "p")
    abline(h = 0, lty = 2, lwd = 0.5)
        })


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
mstop(aic <- AIC(wpbc_surv))
wpbc_surv <- wpbc_surv[mstop(aic)]


###################################################
### chunk number 28: wpbc-glmboost-coef
###################################################
names(coef(wpbc_surv)[abs(coef(wpbc_surv)) > 0])


###################################################
### chunk number 29: wpbc-glmboost-censored-fit
###################################################
plot(log(wpbc3$time), predict(wpbc_surv),
     cex = iw, ylim = c(0, 5), xlim = c(0, 5), 
     xlab = "Time to recurrence (log-scale)", 
     ylab = "Predicted time to recurrence")
abline(a = 0, b = 1, lty = 2, lwd = 0.5)


