
basesel <- function(object, q = floor(sqrt(length(variable.names(object)) / 2)), ...) {
    ibase <- 1:length(variable.names(object))
    fun <- function(model) {
        xs <- selected(model)
        qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
        xs[qq > q] <- xs[1]
        xs
    }
    ss <- cvrisk(object, fun  = selected, ...)
    ret <- matrix(0, nrow = length(ibase), ncol = m <- mstop(object))
    for (i in 1:length(ss)) {
        tmp <- sapply(ibase, function(x) 
            ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
        ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
    }
    phat <- ret / length(ss)
    rownames(phat) <- names(variable.names(object))
    if (extends(class(object), "glmboost"))
        rownames(phat) <- variable.names(object)
    p <- nrow(phat)
    pi <- (q^2 / p + 1) / 2
    attr(phat, "selected") <- which(apply(phat, 1, max) >= pi)
    phat
}

fitsel <- function(object, newdata = NULL, which = NULL, ...) {
    fun <- function(model) {
        tmp <- predict(model, newdata = newdata, 
                       which = which, agg = "cumsum")
        ret <- c()
        for (i in 1:length(tmp))
            ret <- rbind(ret, tmp[[i]])
        ret
    }
    ss <- cvrisk(object, fun = fun, ...)
    ret <- matrix(0, nrow = nrow(ss[[1]]), ncol = ncol(ss[[1]]))
    for (i in 1:length(ss))
        ret <- ret + sign(ss[[i]])
    ret <- abs(ret) / length(ss)
    ret
}
