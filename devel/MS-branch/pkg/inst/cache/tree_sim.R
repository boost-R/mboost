
library("mboost")

### number of simulation runs
B <- 50

### number of learning samples
n <- 100

### number of test samples
m <- 2000

### boosting control parameters
iterations <- 1000
ctrl <- boost_control(mstop = iterations)

### data generating process
dgp <- function(n, d = 20, J = 5) {
    x <- matrix(runif(n * d), ncol = d)
    eta <- 4 * apply(x[,1:J] - 0.5, 1, sum)  
    pp <- exp(eta)/(exp(-eta) + exp(eta))
    y <- as.factor(yraw <- rbinom(n, 1, pp))
    list(x = x, y = y, yraw = yraw, pp = pp)
}

### fast computation of all linear predictors
fastlp <- function(object, newdata) {

    newinp <- party:::newinputs(object$data, newdata)
    lp <- matrix(0, nrow = nrow(newdata), ncol = mstop(object))
    for (m in 1:mstop(object)) {
        wh <- .Call("R_get_nodeID", object$ensemble[[m]], newinp, 0.0,
                    PACKAGE = "party")
        if (m == 1) {
            tmp <- object$offset
        } else {
            tmp <- lp[,m-1]
        }
        lp[,m] <- tmp + object$control$nu * unlist(.Call("R_getpredictions",
            object$ensemble[[m]], wh, PACKAGE = "party"))
    }
    lp
}

### results
misclassblack.logit <- matrix(0, nrow = B, ncol = iterations)
misclassblack.logit.prob <- matrix(0, nrow = B, ncol = iterations)
misclassblack.logit.suloss <- matrix(0, nrow = B, ncol = iterations)

misclass8black.logit <- matrix(0, nrow = B, ncol = iterations)
misclass8black.logit.prob <- matrix(0, nrow = B, ncol = iterations)
misclass8black.logit.suloss <- matrix(0, nrow = B, ncol = iterations)


set.seed(22)

for (b in 1:B) {

    print(b)
    learn <- dgp(n)
    test <- dgp(m)
   
    fitlogit <- blackboost(x = learn$x, y = learn$y, family = Binomial(),
        tree_controls = ctree_control(mincriterion = 0, maxdepth = 1),
        control = boost_control(mstop = iterations))

    fit8logit <- blackboost(x = learn$x, y = learn$y, family = Binomial(),
        tree_controls = ctree_control(mincriterion = 0, maxdepth = 3),
        control = boost_control(mstop = iterations))

    hh <- fastlp(fitlogit, newdata = test$x)
    hh8 <- fastlp(fit8logit, newdata = test$x)

    for (m in 1:iterations) {

        h <- hh[,m]
        misclassblack.logit[b,m] <- mean(as.numeric(h > 0) != test$yraw)
        misclassblack.logit.prob[b,m] <- mean(abs(exp(h)/(exp(h) + exp(-h)) - test$pp))
        misclassblack.logit.suloss[b,m] <- mean(log(1+ exp(-h*2*(2*test$yraw -1)), base=2))

        h <- hh8[,m]
        misclass8black.logit[b,m] <- mean(as.numeric(h > 0) != test$yraw)
        misclass8black.logit.prob[b,m] <- mean(abs(exp(h)/(exp(h) + exp(-h)) - test$pp))
        misclass8black.logit.suloss[b,m] <- mean(log(1+ exp(-h*2*(2*test$yraw -1)), base=2))
    }
}

save(misclassblack.logit, misclassblack.logit.prob,
     misclassblack.logit.suloss, misclass8black.logit, misclass8black.logit.prob,
     misclass8black.logit.suloss, file = "tree_sim.Rda")
