###################################################
### chunk number 1: setup
###################################################
options(width = 50)
require("mboost")
require("Biobase")
### West et al data: binary classification
data("Westbc", package = "mboost")

westbc <- new("ExpressionSet", 
              phenoData = new("AnnotatedDataFrame", data = Westbc$pheno),
              assayData = assayDataNew(exprs = Westbc$assay))

if (!require("kidpack")) {
    install.packages("kidpack", repos = "http://bioconductor.org/packages/monograph")
    require("kidpack")
}
library("survival")
data("westbc", package = "mboost")
data("eset", package = "kidpack")
remove <- is.na(pData(phenoData(eset))$survival.time)
eset <- eset[,!remove]
library("party")


###################################################
### chunk number 2: west-datapp
###################################################
x <- exprs(westbc)
x <- t(x - rowMeans(x))
y <- pData(westbc)$nodal.y


###################################################
### chunk number 3: west-logistic
###################################################
westglm <- glmboost(x, y, family = Binomial(), control = boost_control(mstop = 500))


###################################################
### chunk number 4: west-AIC
###################################################
plot(westAIC <- AIC(westglm, "classical"))


###################################################
### chunk number 5: west-prune
###################################################
westglm <- westglm[mstop(westAIC)]


###################################################
### chunk number 6: west-coef eval=FALSE
###################################################
## coef(westglm)


###################################################
### chunk number 7: west-table
###################################################
table(y, predict(westglm, type = "response"))


###################################################
### chunk number 8: kidpack-datapp
###################################################
x <- exprs(eset)
x <- t(x - rowMeans(x))
y <- Surv(eset$survival.time, eset$died)


###################################################
### chunk number 9: kidpack-Cox
###################################################
kidpackCox <- glmboost(x, y, family = CoxPH())


###################################################
### chunk number 10: kidpack-varsel
###################################################
sum(abs(coef(kidpackCox)) > 0)


###################################################
### chunk number 11: kidpack-predict eval=FALSE
###################################################
## predict(kidpackCox, type = "lp")


###################################################
### chunk number 12: coxgam eval=FALSE
###################################################
## gamboost(x, y, family = CoxPH())


###################################################
### chunk number 13: vignette eval=FALSE
###################################################
## system.file("mboost_Bioinf.R", package = "mboost")


