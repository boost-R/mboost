
require("mboost")

set.seed(290875)

### predict did not return factor levels for blackboost models
dummy <- data.frame(y = gl(2, 100), x = runif(200))
pr <- predict(blackboost(y ~ x, data = dummy, family = Binomial()),
              newdata = dummy, type = "class")
stopifnot(is.factor(pr) && all(levels(pr) %in% levels(dummy$y)))

### predict for g{al}mboost.matrix did not work
ctrl <- boost_control(mstop = 10)
X <- cbind(int = 1, x = dummy$x)
gb <- glmboost(x = X, y = dummy$y, family = Binomial(),
               control = ctrl)
stopifnot(all.equal(predict(gb), predict(gb, newdata = X)))

if (FALSE) {
gb <- gamboost(x = X, y = dummy$y, family = Binomial(),
               control = ctrl)
stopifnot(all.equal(predict(gb), predict(gb, newdata = X)))
}

### blackboost _did_ touch the response, arg!

data("bodyfat", package = "mboost")
ctrl <- boost_control(mstop = 500, nu = 0.01)
bb <- blackboost(DEXfat ~ ., data = bodyfat, control = ctrl)
n <- nrow(bodyfat)
bs <- rmultinom(3, n, rep(1, n) / n)
x <- seq(from = 10, to = 500, by = 10)
cv <- cvrisk(bb, bs, grid = x)
ctrl$risk <- "oobag"
tmp <- blackboost(DEXfat ~ ., data = bodyfat, control = ctrl,
                 weights = bs[,3])

stopifnot(identical(max(abs(tmp$risk()[x] / sum(bs[,3] == 0)  - cv[3,])), 0))

### center = TRUE and cvrisk were broken; same issue with masking original data

gb <- glmboost(DEXfat ~ ., data = bodyfat, control = boost_control(center = TRUE))
cv1 <- cvrisk(gb, folds = bs)
tmp <- glmboost(DEXfat ~ ., data = bodyfat,
                control = boost_control(center = TRUE, risk = "oobag"),
                weights = bs[,3])
stopifnot(identical(max(tmp$risk()[attr(cv1, "mstop")] / sum(bs[,3] == 0) - cv1[3,]), 0))

### same problem, just another check

indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
cbodyfat <- bodyfat
cbodyfat[indep] <- lapply(cbodyfat[indep], function(x) x - mean(x))
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4

bf_glm_1 <- glmboost(bffm, data = cbodyfat)
cv1 <- cvrisk(bf_glm_1, folds = bs)
bf_glm_2 <- glmboost(bffm, data = bodyfat, control = boost_control(center = TRUE))
cv2 <- cvrisk(bf_glm_2, folds = bs)

stopifnot(mstop(cv1) == mstop(cv2))

### dfbase=1 was not working correctly for ssp
### spotted by Matthias Schmid <Matthias.Schmid@imbe.imed.uni-erlangen.de>
data("bodyfat", package = "mboost")
ctrl <- boost_control(mstop = 100, center = TRUE)
ga <- gamboost(DEXfat ~ ., data = bodyfat, dfbase = 1, control = ctrl)
gl <- glmboost(DEXfat ~ ., data = bodyfat, control = ctrl)
stopifnot(max(abs(predict(ga) - predict(gl))) < 1e-8)
AIC(gl)

if (FALSE) {
### prediction with matrices was broken for gamboost,
### spotted by <Max.Kuhn@pfizer.com>
x <- matrix(rnorm(1000), ncol = 10)
colnames(x) <- 1:ncol(x)
y <- rnorm(100)
fit <- gamboost(x = x, y = y, control = boost_control(mstop = 10))
a <- predict(fit, newdata = x[1:10,])
}

### AIC for centered covariates didn't work
y <- gl(2, 10)
xn <- rnorm(20)
xnm <- xn - mean(xn)
xf <- gl(2, 10)
gc <- glmboost(y ~ xn + xf, control = boost_control(center = TRUE),
               family = Binomial())
g <- glmboost(y ~ xnm + xf, family = Binomial())
cgc <- coef(gc)
cg <- coef(g)
names(cgc) <- NULL
names(cg) <- NULL
stopifnot(all.equal(cgc, cg))
stopifnot(all.equal(mstop(AIC(gc, "classical")), mstop(AIC(g, "classical"))))

gc <- gamboost(y ~ xn + bols(xf), control = boost_control(center = TRUE),
               family = Binomial())
g <- gamboost(y ~ xnm + bols(xf), family = Binomial())
stopifnot(all.equal(mstop(AIC(gc, "classical")), mstop(AIC(g, "classical"))))

y <- rnorm(20)
gc <- gamboost(y ~ xn + bols(xf), control = boost_control(center = TRUE))
g <- gamboost(y ~ xnm + bols(xf))
stopifnot(all.equal(mstop(AIC(gc, "corrected")), mstop(AIC(g, "corrected"))))

### check gamboost with weights (use weighted some of residuals
### for variable selection)
data("bodyfat", package = "mboost")

set.seed(290875)
n <- nrow(bodyfat)
w <- numeric(n)
w <- rmultinom(1, n, rep(1, n) / n)[,1]
ctrl <- boost_control(mstop = 20)

mod1 <- glmboost(DEXfat ~ ., data = bodyfat, weights = w)
aic1 <- AIC(mod1, "corrected")
attributes(aic1) <- NULL

mod2 <- glmboost(DEXfat ~ ., data = bodyfat[rep(1:n, w),])
aic2 <- AIC(mod2, "corrected")
attributes(aic2) <- NULL

stopifnot(all.equal(round(aic1, 3), round(aic2, 3)))

mod1 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl)
aic1 <- AIC(mod1, "corrected")
attributes(aic1) <- NULL

mod2 <- gamboost(DEXfat ~ ., data = bodyfat[rep(1:n, w),], con = ctrl)
aic2 <- AIC(mod2, "corrected")
attributes(aic2) <- NULL

stopifnot(all.equal(round(aic1, 1), round(aic2, 1)))

mod1 <- blackboost(DEXfat ~ ., data = bodyfat, weights = w)
mod2 <- blackboost(DEXfat ~ ., data = bodyfat[rep(1:n, w),])

ratio <- mod1$risk() / mod2$risk()
stopifnot(ratio[1] > 0.95 && ratio[2] < 1.05)

### df <= 2
ctrl$risk <- "oobag"
mod1 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl, base = "bss", dfbase = 2)
mod2 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl, base = "bbs", dfbase = 2)
mod3 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl, base = "bols")
stopifnot(max(abs(predict(mod1) - predict(mod2))) < .Machine$double.eps)
# stopifnot(max(abs(predict(mod1) - predict(mod3))) < .Machine$double.eps)

### check predictions of zero-weight observations (via out-of-bag risk)
mod1 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl, base = "bss")
mod2 <- gamboost(DEXfat ~ ., data = bodyfat, weights = w, con = ctrl, base = "bbs")
stopifnot(abs(coef(lm(mod1$risk() ~ mod2$risk() - 1)) - 1) < 0.01)

### not really a bug, a new feature; test fastp
df <- data.frame(y = rnorm(100), x = runif(100), z = runif(100))
eps <- sqrt(.Machine$double.eps)
s <- seq(from = 1, to = 100, by = 3)

x <- glmboost(y ~ ., data = df)
for (i in s)
    stopifnot(max(abs(predict(x[i]) - predict(x[max(s)], agg = "cumsum")[,i])) < eps)

x <- gamboost(y ~ ., data = df)
for (i in s)
    stopifnot(max(abs(predict(x[i]) - predict(x[max(s)], agg = "cumsum")[,i])) < eps)

x <- blackboost(y ~ ., data = df)
for (i in s)
    stopifnot(max(abs(predict(x[i]) - predict(x[max(s)], agg = "cumsum")[,i])) < eps)

### negative gradient of GaussClass was incorrectly specified
### negative gradient of GaussClass was incorrectly specified
data("BreastCancer", package = "mlbench")
tmp <- BreastCancer[complete.cases(BreastCancer), -1]
learn <- sample(1:nrow(tmp), ceiling(nrow(tmp) * 0.7))

stump <- blackboost(Class ~ ., data = tmp[learn,],
    tree_controls = ctree_control(teststat = "max",
        testtype = "Teststatistic", mincriterion = 0, stump = TRUE),
    control = boost_control(mstop = 176), family = GaussClass())
mean(predict(stump, newdata = tmp[-learn,], type = "class") != tmp[-learn, "Class"])

stump <- blackboost(Class ~ ., data = tmp[learn,],
    tree_controls = ctree_control(teststat = "max",
        testtype = "Teststatistic", mincriterion = 0, stump = TRUE),
    control = boost_control(mstop = 275), family = GaussClass())
mean(predict(stump, newdata = tmp[-learn,], type = "class") != tmp[-learn, "Class"])

cspline <- gamboost(Class ~ ., data = tmp[learn,],
    control = boost_control(mstop = 126), family = GaussClass())
mean(predict(cspline, newdata = tmp[-learn,], type = "class") != tmp[-learn, "Class"])

cspline <- gamboost(Class ~ ., data = tmp[learn,],
    control = boost_control(mstop = 73), family = GaussClass())
mean(predict(cspline, newdata = tmp[-learn,], type = "class") != tmp[-learn, "Class"])

### make sure environment(formula) is used for evaluation
data("cars")
ctl  <- boost_control(mstop = 100, trace = TRUE)
tctl <- ctree_control(teststat = "max", testtype = "Teststat",
                      mincrit = 0, maxdepth = 5, savesplitstat = FALSE)
myfun <- function(cars, xx, zz){
  mboost(dist ~ btree(speed, tree_controls = zz),
         data = cars, control = xx)
}
try(mod <- myfun(cars, xx = ctl, zz = tctl))

### bbs with weights and expanded data
x <- runif(100)
y <- rnorm(length(x))
knots <- seq(from = 0.1, to = 0.9, by = 0.1)
w <- rmultinom(1, length(x), rep(1, length(x)) / length(x))[,1]
iw <- rep(1:length(x), w)

m1 <- bbs(x, knots = knots)$dpp(w)$fit(y)$model
m2 <- bbs(x[iw], knots = knots)$dpp(rep(1, length(iw)))$fit(y[iw])$model

stopifnot(max(abs(m1 - m2)) < sqrt(.Machine$double.eps))

### base learner handling
stopifnot(max(abs(fitted(gamboost(DEXfat ~ age, data = bodyfat)) -
                  fitted(gamboost(DEXfat ~ bbs(age), data = bodyfat)))) <
          sqrt(.Machine$double.eps))

### predict for matrix interface to glmboost
x <- matrix(runif(1000), ncol = 10)
y <- rowMeans(x) + rnorm(nrow(x))
mod <- glmboost(x = x, y = y)
stopifnot(length(predict(mod, newdata = x[1:2,])) == 2)
try(predict(mod, newdata = as.data.frame(x[1:2,])))


### predict for varying coefficient models with categorical z
set.seed(1907)
x <- rnorm(100)
z <- gl(2,50)
y <- rnorm(100, sd = 0.1)
data <- data.frame(y=y, x=x, z=z)

model <- gamboost(y ~ bols(x, by=z), data = data)
stopifnot(!is.na(predict(model,data[1,-1])))

model <- gamboost(y ~ bbs(x, by=z), data = data)
stopifnot(!is.na(predict(model,data[1,-1])))

x <- as.factor(sample(1:10, 100, replace=TRUE))
data <- data.frame(y=y, x=x, z=z)
model <- gamboost(y ~ brandom(x, by=z), data = data)
stopifnot(!is.na(predict(model,data[1,-1])))

x1 <- rnorm(100)
x2 <- rnorm(100)
data <- data.frame(y=y, x1=x1, x2=x2, z=z)
model <- gamboost(y ~ bspatial(x1,x2, by=z), data = data)
stopifnot(!is.na(predict(model,data[1,-1])))

### bols with intercept = FALSE for categorical covariates was broken
x <- gl(2, 50)
y <- c(rnorm(50, mean = -1), rnorm(50, mean = 1))
stopifnot(length(coef(gamboost(y ~ bols(x, intercept=FALSE)))) == 1)

# check also for conntinuous x
x <- rnorm(100)
y <- c(rnorm(100, mean = 1 * x))
stopifnot(length(coef(gamboost(y ~ bols(x, intercept=FALSE)))) == 1)

### check interface of coef
set.seed(1907)
x1 <- rnorm(100)
int <- rep(1, 100)
y <- 3 * x1 + rnorm(100, sd=0.1)
dummy <- data.frame(y = y, int = int, x1 = x1)

gbm <- gamboost(y ~ bols(int, intercept=FALSE) +  bols(x1, intercept=FALSE) + bbs(x1, center=TRUE, df=1), data = dummy)

stopifnot(names(coef(gbm, which=1:3)) == c("bols(int,intercept=FALSE)", "bols(x1,intercept=FALSE)", "bbs(x1,center=TRUE,df=1)"))
stopifnot(names(coef(gbm)) == c("bols(x1,intercept=FALSE)", "bbs(x1,center=TRUE,df=1)"))
stopifnot(names(coef(gbm, "x1")) == c("bols(x1,intercept=FALSE)", "bbs(x1,center=TRUE,df=1)"))
stopifnot(names(coef(gbm, "bbs")) ==  "bbs(x1,center=TRUE,df=1)")
stopifnot(names(coef(gbm, "center=TRUE")) == "bbs(x1,center=TRUE,df=1)")

### check prediction if intercept=FALSE
gbm <- gamboost(y ~ bols(x1, intercept=FALSE), data = dummy)
stopifnot(!is.na(predict(gbm)) & max(abs(predict(gbm) - fitted(gbm))) < sqrt(.Machine$double.eps))
