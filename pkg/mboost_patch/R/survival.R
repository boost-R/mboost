
survFit <- function(object, ...)
    UseMethod("survFit")

survFit.mboost <- function(object, newdata = NULL, ...)
{

    n <- length(w <- model.weights(object))
    if (!all.equal(w,rep(1,n)))
        stop("survFit cannot (yet) deal with weights")
    y <- object$response
    stopifnot(inherits(y, "Surv"))

    ord <- order(y[,1])
    time <- y[ord,1]
    event <- y[ord,2]
    n.event <- aggregate(event, by = list(time), sum)[,2]

    pr <- predict(object)
    linpred <- (pr - mean(pr))[ord]
    R <- rev(cumsum(rev(exp(linpred))))
    R[R<1e-6] <- Inf

    d <- c(1,diff(time)) != 0
    R <- R[d]
    H <- cumsum(1/R*n.event)
    time <- time[d]

    devent <- n.event > 0
    if (!is.null(newdata)){
        S <- exp( tcrossprod( -H, exp(as.numeric(predict(object,
        newdata=newdata)- mean(predict(object)))) ))[devent,]
        colnames(S) <- rownames(newdata)
    } else S <- matrix(exp(-H)[devent], ncol = 1)

    ret <- list(surv = S, time = time[devent], n.event = n.event[devent])
    class(ret) <- "survFit"
    ret
}

plot.survFit <- function(x, xlab = "Time", ylab = "Probability", ...) {
    plot(x$time, rep(1, length(x$time)), ylim = c(0, 1), type = "n",
         ylab = ylab, xlab = xlab, ...)
    tmp <- apply(rbind(1,x$surv), 2, function(s) lines(c(0,x$time), s,
                 type = "s"))
}


### inverse probability of censoring weights
### see van der Laan & Robins (2003)
IPCweights <- function(x, maxweight = 5) {

    if (!extends(class(x), "Surv"))
        stop(sQuote("x"), " is not a Surv object")

    event <- x[,"status"]
    x[,"status"] <- 1 - event
    km <- survfit(x ~ 1)
    Ghat <- getsurv(km, times = x[,"time"])
    Ghat[event == 0] <- 1
    w <- event / Ghat
    w[w > maxweight] <- maxweight
    w
}

### extract survival probabilities
### taken from ipred:::getsurv
### DO NOT TOUCH HERE
getsurv <- function(obj, times)
{
    # get the survival probability for times from KM curve j'

    if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
    # <FIXME: methods may have problems with that>
    class(obj) <- NULL
    # </FIXME>
    lt <- length(times)
    nsurv <- times

    # if the times are the same, return the km-curve

    if(length(times) == length(obj$time)) {
        if (all(times == obj$time)) return(obj$surv)
    }

    # otherwise get the km-value for every element of times separatly

    inside <- times %in% obj$time
    for (i in (1:lt)) {
        if (inside[i])
            nsurv[i] <- obj$surv[obj$time == times[i]]
        else  {
            less <- obj$time[obj$time < times[i]]
            if (length(less) == 0)
                nsurv[i] <- 1
            else
                nsurv[i] <- obj$surv[obj$time == max(less)]
        }
    }
    nsurv
}
