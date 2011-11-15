
library("mboost")
data("birds", package = "mboost")

# define characteristics of the boosting algorithm
bcr <- boost_control(mstop=200, trace=TRUE)

# estimation of a purely linear GLM
fm <- SG5 ~ bols(GST) + bols(DBH) + bols(AOT) + bols(AFS) + bols(DWC) +
            bols(LOG)
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)

# extract and plot AIC curve against iteration index and determine stopping
# iteration
birdsaic <- AIC(sp, "classical")
plot(birdsaic)
ms <- mstop(birdsaic)

# selection frequencies of the model terms
table(sp$xselect()[1:ms])

# estimated coefficients
coef(sp[ms])

# re-define boosting iterations
bcr <- boost_control(mstop=500, trace=TRUE)

# Variable selection in a GLM without spatial component
fm <- SG4 ~ bols(GST) + bols(DBH) + bols(AOT) + bols(AFS) + bols(DWC) +
            bols(LOG)
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
table(sp$xselect())
coef(sp, which=1:6)

# Variable selection in a GLM with high df spatial component
fm <- SG4 ~ bols(GST) + bols(DBH) + bols(AOT) + bols(AFS) + bols(DWC) +
            bols(LOG) + bspatial(x_gk, y_gk, df=5, differences=1, knots=c(12,12))
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
table(sp$xselect())
coef(sp, which=1:6)

# Variable selection in a GLM with small df spatial component
fm <- SG4 ~ bols(GST) + bols(DBH) + bols(AOT) + bols(AFS) + bols(DWC) +
            bols(LOG) + bspatial(x_gk, y_gk, df=1, differences=1, knots=c(12,12), center=TRUE)
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
table(sp$xselect())
coef(sp, which=1:6)


# Geoadditive regression model without centering
fm <- SG5 ~ bbs(GST) + bbs(DBH) + bbs(AOT) + bbs(AFS) + bbs(DWC) +
            bbs(LOG) + bspatial(x_gk, y_gk, df=4, differences=1, knots=c(12,12))
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
plot(sp)

# Geoadditive regression model with centering

fm <- SG5 ~ bols(GST) + bbs(GST, df=1, center=TRUE) +
            bols(AOT) + bbs(AOT, df=1, center=TRUE) +
            bols(AFS) + bbs(AFS, df=1, center=TRUE) +
            bols(DWC) + bbs(DWC, df=1, center=TRUE) +
            bols(LOG) + bbs(LOG, df=1, center=TRUE) +
            bspatial(x_gk, y_gk, df=1, differences=1, knots=c(12,12),
              center=TRUE)
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
plot(sp)


# re-define boosting iterations
bcr <- boost_control(mstop=200, trace=TRUE)

# transform covariates to [0,1]
birds$GST <- (birds$GST-min(birds$GST))/(max(birds$GST)-min(birds$GST))
birds$AOT <- (birds$AOT-min(birds$AOT))/(max(birds$AOT)-min(birds$AOT))
birds$AFS <- (birds$AFS-min(birds$AFS))/(max(birds$AFS)-min(birds$AFS))
birds$DWC <- (birds$DWC-min(birds$DWC))/(max(birds$DWC)-min(birds$DWC))
birds$LOG <- (birds$LOG-min(birds$LOG))/(max(birds$LOG)-min(birds$LOG))

# Space-varying coefficient models (with centered spatial effects)
fm <- SG5 ~ bols(GST) + bspatial(x_gk, y_gk, by = GST, df=1, differences=1,
              knots=c(12, 12), center=TRUE) +
            bols(AOT) + bspatial(x_gk, y_gk, by = AOT, df=1, differences=1,
              knots=c(12, 12), center=TRUE) +
            bols(AFS) + bspatial(x_gk, y_gk, by = AFS, df=1, differences=1,
              knots=c(12, 12), center=TRUE) +
            bols(DWC) + bspatial(x_gk, y_gk, by = DWC, df=1, differences=1,
              knots=c(12, 12), center=TRUE) +
            bols(LOG) + bspatial(x_gk, y_gk, by = LOG, df=1, differences=1,
              knots=c(12, 12), center=TRUE) +
            bspatial(x_gk, y_gk, df=1, differences=1, knots=c(12, 12),
              center=TRUE)
sp <- gamboost(fm, data = birds, family = Poisson(), control = bcr)
plot(sp, which = "GST")
plot(sp, which = "AOT")
plot(sp, which = "AFS")
plot(sp, which = "DWC")
plot(sp, which = "LOG")

