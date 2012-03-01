
library("mboost")
source("india_helpfunc.R")

QuantileReg <- function(q = 0.25) {
   Family(
   ngradient = function(y, f, w = 1) q*((y-f)>0) - (1-q)*((y-f)<0) + 0*((y-f)==0) ,
   loss = function(y, f) q*(y-f)*((y-f)>=0) - (1-q)*(y-f)*((y-f)<0) ,
   offset = function(y, w = rep(1, length(y))) quantile(y[w==1], q),
   check_y = function(y) {
     if (!is.numeric(y) || !is.null(dim(y)))
         stop("response is not a numeric vector but ", sQuote("family = QuantileReg()"))
     TRUE
     },
   name = "Quantile Regression",
   )
 }

fm <- stunting ~ bbs(cage, center=FALSE, df=5, knots = 20) +
  bbs(breastfeeding, center=FALSE, df=5, knots=20) +
  bols(csex, center=FALSE) +
  bols(ctwin, center=FALSE) +
  bols(cbirthorder, center=FALSE) +
  bbs(mbmi, center=FALSE, df=5, knots=20) +
  bbs(mage, center=FALSE, df=5, knots=20) +
  bbs(medu, center=FALSE, df=5, knots=20) +
  bbs(edupartner, center=FALSE, df=5, knots=20) +
  bols(munemployed, center=FALSE) +
  bols(mreligion, center=FALSE) +
  bols(mresidence, center=FALSE) +
  bols(deadchildren, center=FALSE) +
  bols(wealth, center=FALSE) +
  bols(electricity, center=FALSE) +
  bols(radio, center=FALSE) +
  bols(television, center=FALSE) +
  bols(refrigerator, center=FALSE) +
  bols(bicycle, center=FALSE) +
  bols(motorcycle, center=FALSE) +
  bols(car, center=FALSE)

inb1 <- gamboost(fm,
                    data = india,
                    control = boost_control(mstop = 200000, nu=0.25, trace = TRUE, save_ensembless=FALSE, risk = "oob"),
                    weights = india$cv,
                    family = QuantileReg(0.05))
mstop1 <- which.min(inb1$risk)
its1 <- round(seq(from=1000, to=mstop1, length=20))
fit1 <- extractfit(inb1, fm, its1)

inb2 <- gamboost(fm,
                    data = india,
                    control = boost_control(mstop = 200000, nu=0.25, trace = TRUE, save_ensembless=FALSE, risk = "oob"),
                    weights = india$cv,
                    family = QuantileReg(0.1))
mstop2 <- which.min(inb2$risk)
its2 <- round(seq(from=1000, to=mstop2, length=20))
fit2 <- extractfit(inb2, fm, its2)

inb3 <- gamboost(fm,
                    data = india,
                    control = boost_control(mstop = 200000, nu=0.25, trace = TRUE, save_ensembless=FALSE, risk = "oob"),
                    weights = india$cv,
                    family = QuantileReg(0.5))
mstop3 <- which.min(inb3$risk)
its3 <- round(seq(from=1000, to=mstop3, length=20))
fit3 <- extractfit(inb3, fm, its3)

save(inb1, fit1, inb2, fit2, inb3, fit3, file="india_results.Rdata")

pdf("india_results1.pdf", paper="a4", width=8, height=12)
par(mfrow=c(3,2))
plotnonpar("cage", 1, "age of the child", fit1)
plotnonpar("breastfeeding", 2, "duration of breastfeeding", fit1)
plotnonpar("mbmi", 6, "body mass index of the mother", fit1)
plotnonpar("mage", 7, "age of the mother", fit1)
plotnonpar("medu", 8, "years of education (mother)", fit1)
plotnonpar("edupartner", 9, "years of education (partner)", fit1)

plotpar("csex", 3, "sex of the child", fit1)
plotpar("ctwin", 4, "twin", fit1)
plotpar("cbirthorder", 5, "birth order", fit1)
plotpar("munemployed", 10, "unemployment of the mother", fit1)
plotpar("mreligion", 11, "religion", fit1)
plotpar("mresidence", 12, "place of residence", fit1)
plotpar("deadchildren", 13, "no. of dead children", fit1)
plotpar("wealth", 14, "wealth index", fit1)
plotpar("electricity", 15, "electricity", fit1)
plotpar("radio", 16, "radio", fit1)
plotpar("television", 17, "television", fit1)
plotpar("refrigerator", 18, "refrigerator", fit1)
plotpar("bicycle", 19, "bicycle", fit1)
plotpar("motorcycle", 20, "motorcycle", fit1)
plotpar("car", 21, "car", fit1)
dev.off()

pdf("india_results2.pdf", paper="a4", width=8, height=12)
par(mfrow=c(3,2))
plotnonpar("cage", 1, "age of the child", fit2)
plotnonpar("breastfeeding", 2, "duration of breastfeeding", fit2)
plotnonpar("mbmi", 6, "body mass index of the mother", fit2)
plotnonpar("mage", 7, "age of the mother", fit2)
plotnonpar("medu", 8, "years of education (mother)", fit2)
plotnonpar("edupartner", 9, "years of education (partner)", fit2)

plotpar("csex", 3, "sex of the child", fit2)
plotpar("ctwin", 4, "twin", fit2)
plotpar("cbirthorder", 5, "birth order", fit2)
plotpar("munemployed", 10, "unemployment of the mother", fit2)
plotpar("mreligion", 11, "religion", fit2)
plotpar("mresidence", 12, "place of residence", fit2)
plotpar("deadchildren", 13, "no. of dead children", fit2)
plotpar("wealth", 14, "wealth index", fit2)
plotpar("electricity", 15, "electricity", fit2)
plotpar("radio", 16, "radio", fit2)
plotpar("television", 17, "television", fit2)
plotpar("refrigerator", 18, "refrigerator", fit2)
plotpar("bicycle", 19, "bicycle", fit2)
plotpar("motorcycle", 20, "motorcycle", fit2)
plotpar("car", 21, "car", fit2)
dev.off()

pdf("india_results3.pdf", paper="a4", width=8, height=12)
par(mfrow=c(3,2))
plotnonpar("cage", 1, "age of the child", fit3)
plotnonpar("breastfeeding", 2, "duration of breastfeeding", fit3)
plotnonpar("mbmi", 6, "body mass index of the mother", fit3)
plotnonpar("mage", 7, "age of the mother", fit3)
plotnonpar("medu", 8, "years of education (mother)", fit3)
plotnonpar("edupartner", 9, "years of education (partner)", fit3)

plotpar("csex", 3, "sex of the child", fit3)
plotpar("ctwin", 4, "twin", fit3)
plotpar("cbirthorder", 5, "birth order", fit3)
plotpar("munemployed", 10, "unemployment of the mother", fit3)
plotpar("mreligion", 11, "religion", fit3)
plotpar("mresidence", 12, "place of residence", fit3)
plotpar("deadchildren", 13, "no. of dead children", fit3)
plotpar("wealth", 14, "wealth index", fit3)
plotpar("electricity", 15, "electricity", fit3)
plotpar("radio", 16, "radio", fit3)
plotpar("television", 17, "television", fit3)
plotpar("refrigerator", 18, "refrigerator", fit3)
plotpar("bicycle", 19, "bicycle", fit3)
plotpar("motorcycle", 20, "motorcycle", fit3)
plotpar("car", 21, "car", fit3)
dev.off()