
library("mboost")

x <- rnorm(100)
y <- rnorm(100)
w <- drop(rmultinom(1, 100, rep(1 / 100, 100)))

G <- Gaussian()
fm <- Family(ngradient = G@ngradient, risk = G@risk)

glmboost(y ~ x, family = G)

glmboost(y ~ x, family = fm)


x <- rnorm(100)
y <- rnbinom(length(x), size = 2, mu = exp(x * 2))
mod <- glmboost(y ~ x, family = NBinomial())
mod[1000]
coef(mod)
nuisance(mod)

if (require("MASS")) {

summary(glm.nb(y ~ x))

y <- cut(x, breaks = c(-Inf, quantile(x, prob = c(0.25, 0.5, 0.75)), Inf), ordered = TRUE)
x <- rnorm(100)
polr(y ~ x)
mod <- glmboost(y ~ x, family = PropOdds())
nuisance(mod) - attr(coef(mod), "offset")
coef(mod)

}

