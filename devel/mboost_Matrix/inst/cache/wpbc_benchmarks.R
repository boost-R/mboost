
library("mboost")

### Wisconsin prognostic breast cancer data, without missing values
data("wpbc")
wpbc <- wpbc[complete.cases(wpbc),]

indep <- names(wpbc)[!(names(wpbc) %in% c("time", "status"))]
wpbc[indep] <- lapply(wpbc[indep], function(x) x - mean(x))

### setup benchmark experiments with B = 100 bootstrap samples
n <- nrow(wpbc)
set.seed(290875)
bs <- rmultinom(100, n, rep(1, n)/n)
obj <- boost_dpp(status ~ ., data = wpbc[,-2])

boob <- c()
grid <- seq(from = 5, to = 500, by = 5)

### negative binom. log-lik / n for different number of boosting iterations
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f) 
        Binomial()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- glmboost_fit(obj, family = Binomial(),
        control = boost_control(mstop = 500), weights = b)

    gm <- mboost:::gm.glmboost(object)
    fm <- t(apply(gm, 1, cumsum))[,grid]
    boob <- rbind(boob, apply(fm, 2, function(f) foo(f)))
}

goob <- c()

flink <- function(x)
    pmin(abs(x), 18) * sign(x)

### negative binom. log-lik / n for step(glm(...))
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f) 
        Binomial()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- glm(status ~  ., data = wpbc[,-2], subset = b > 0, 
                  family = binomial(), weights = b)
    p <- predict(object, newdata = wpbc[,-2], type = "link")
    goob <- c(goob, foo(flink(p)))
}

soob <- c()

flink <- function(x)
    pmin(abs(x), 18) * sign(x)

### negative binom. log-lik / n for step(glm(...))
for (j in 1:ncol(bs)) {
    b <- bs[,j]
    print(j)
    foo <- function(f) 
        Binomial()@risk(obj$yfit, f, as.numeric(b == 0)) / sum(b == 0)
    object <- step(glm(status ~  ., data = wpbc[,-2], subset = b > 0, 
                       family = binomial(), weights = b),
                       trace = 0)
    p <- predict(object, newdata = wpbc[,-2], type = "link")
    soob <- c(soob, foo(flink(p)))
}



save(obj, bs, goob, soob, boob, file = "wpbc_benchmarks.rda")

plot(grid, colMeans(boob), xlim = range(grid), type = "l")
print(mean(goob))

bmod <- glmboost_fit(obj, family = Binomial(), control = boost_control(mstop = 500))
### AIC(bmod, "classical")
gmod <- step(glm(status ~  ., data = wpbc[,-2], family = binomial()), trace = FALSE)

boxplot(as.data.frame(boob), pars = list(boxwex = 0.2, staplewex = 0.1, outwex = 0.1), axes = FALSE,
        ylim = c(-0.65, -0.44), xlab = "Number of Boosting Iterations", 
        ylab = "Negative Binomial Log-Likelihood / n")
axis(1, at = 1:ncol(boob), labels = grid)
axis(2)
box()
abline(h = logLik(bmod[415]) / nrow(wpbc), lty = 2)
abline(h = logLik(gmod) / nrow(wpbc), lty = 3)
legend(0, -0.44, legend = c("Logistic Regression (step)", "AIC-stopped boosting"), 
       lty = c(3, 2), bty = "n")
