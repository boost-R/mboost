### just a try
plot.mboost <- function(x, which = NULL, newdata = NULL,
                        type = "b", rug = TRUE, rugcol = "black", ylim = NULL,
                        xlab = NULL, ylab = expression(f[partial]), add = FALSE,
                        ...) {

    if (inherits(x, "blackboost"))
        stop("partial dependency plots for ", sQuote("blackboost"), " not yet implemented.",
             "See ?blackboost.")

    which <- x$which(which, usedonly = is.null(which))

    if (is.null(xlab)){
        userspec <- FALSE
        xlab <- variable.names(x)
    } else {
        userspec <- TRUE
        if (length(which) == length(xlab) & length(which) != length(variable.names(x))){
            ## if labels are user specified but not all base-learners are selected
            foo <- rep(NA, length(variable.names(x)))
            foo[which] <- xlab ## xlab is supposed to be in the same ordering as which
            xlab <- foo
        }
    }
    if (length(which) == length(ylab) & length(which) != length(variable.names(x))){
        ## if labels are user specified but not all base-learners are selected
        foo <- rep(NA, length(variable.names(x)))
        foo[which] <- ylab ## ylab is supposed to be in the same ordering as which
        ylab <- foo
    }

    ## labs must have either length 1 or length(variable.names(x))
    stopifnot(length(xlab) %in% c(1, length(variable.names(x))))
    stopifnot(length(ylab) %in% c(1, length(variable.names(x))))
    ## FIXED?

    RET <- vector("list", length = length(which))

    for (w in which) {

        data <- model.frame(x, which = w)[[1]]
        get_vary <- x$baselearner[[w]]$get_vary
        vary <- ""
        if (!is.null(get_vary))
            vary <- get_vary()
        if (!is.null(newdata)) {
            data <- newdata[colnames(data)]
            if (is.list(data))
                data <- as.data.frame(data)
        }
        if (vary != "") {
            v <- data[[vary]]
            if (is.factor(v)) v <- factor(levels(v)[-1], levels = levels(v))
            if (is.numeric(v)) v <- 1
        }

        plot_helper <- function(xl, yl){
            pr <- predict(x, newdata = data, which = w)
            if (is.null(ylim)) ylim <- range(pr, na.rm = TRUE)

            if (vary != "") {
                datavary <- data[, colnames(data) == vary, drop = FALSE]
                data <- data[, colnames(data) != vary, drop = FALSE]
            }
            if (ncol(data) == 1) {
                if (!add) {
                    if (is.factor(data[[1]])) {
                        xVals <- unique(sort(data[[1]]))
                        xValsN <- as.numeric(xVals)
                        ## make sure that we get the same number of values as in
                        ## x; this is only a problem if pr is equal for
                        ## different x values.
                        yVals <- unique(cbind(pr[order(data[[1]], na.last = NA)],
                                              sort(data[[1]])))[, 1]

                        if (length(pr) == 1 && pr == 0) {
                            yVals <- rep(0, length(xVals))
                        }
                        plot(xValsN, yVals,
                             type = "n", xaxt = "n",
                             xlim = range(as.numeric(xVals)) + c(-0.5, 0.5),
                             xlab = xl, ylab = yl, ylim = ylim)
                        axis(1, at = xValsN, labels = levels(xVals))
                        for (i in 1:length(xVals)) {
                            lines(x = rep(xValsN[i], 2) + c(-0.35, 0.35),
                                  y = rep(yVals[i], 2), ...)
                        }
                    } else {
                        plot(sort(data[[1]]), pr[order(data[[1]], na.last = NA)], type = type,
                             xlab = xl, ylab = yl, ylim = ylim, ...)
                    }
                    if (rug) rug(data[[1]], col = rugcol)
                } else {
                    if (is.factor(data[[1]])){
                        xVals <- unique(sort(data[[1]]))
                        xValsN <- as.numeric(xVals)
                        yVals <- unique(pr[order(data[[1]], na.last = NA)])
                        if (length(pr) == 1 && pr == 0) {
                            yVals <- rep(0, length(xVals))
                        }
                        axis(1, at = xValsN, labels = levels(xVals))
                        for (i in 1:length(xVals)) {
                            lines(x = rep(xValsN[i], 2) + c(-0.35, 0.35),
                                  y = rep(yVals[i], 2), ...)
                        }
                    } else {
                        lines(sort(data[[1]]), pr[order(data[[1]], na.last = NA)], type =
                              type, ...)
                        if (rug){
                            rug(data[[1]], col = rugcol)
                            warning(sQuote("rug  =TRUE"),
                                    " should be used with care if ",
                                    sQuote("add = TRUE"))
                        }
                    }
                }
            }
            if (ncol(data) == 2) {
                if (is.null(newdata)){
                    tmp <- expand.grid(unique(data[,1]), unique(data[,2]))
                    colnames(tmp) <- colnames(data)
                    data <- tmp
                    if (vary != "") {
                        ## datavary contains only 1 value, thus use only first
                        ## entry and use recycling to get appropriate vector
                        tmp <- cbind(datavary[1,], tmp)
                        colnames(tmp)[1] <- vary
                    }
                    pr <- predict(x, newdata = tmp, which = w)
                }
                fm <- as.formula(paste("pr ~ ", paste(colnames(data), collapse = "*"), sep = ""))
                RET[[w]] <<- levelplot(fm, data = data, ...)
            }
            if (ncol(data) > 2) {
                for (v in colnames(data)) {
                    tmp <- data
                    pardep <- sapply(data[[v]], function(vv) {
                        tmp[[v]] <- vv
                        mean(predict(x, newdata = tmp, which = w))
                    })
                    plot(sort(data[[v]]), pardep[order(data[[v]], na.last = NA)], type = type,
                         xlab = v, ylab = "Partial Dependency", ylim = ylim, ...)
                }
            }
        } # END plot_helper()

        xl <- ifelse(length(xlab) > 1, xlab[w], xlab[1])
        yl <- ifelse(length(ylab) > 1, ylab[w], ylab[1])

        if (add && length(which) > 1)
            warning(sQuote("add = TRUE"),
                    " should be used for single plots only")

        if (!isMATRIX(data)){
            if (vary == "") plot_helper(xl, yl)
            if (vary != ""){
                for (i in 1:length(v)){
                    data[[vary]] <- v[rep(i, nrow(data))]
                    if (!userspec){
                        ## xlab not user specified
                        plot_helper(paste(xl, "=", v[i]), yl)
                    } else {
                        plot_helper(xl, yl)
                    }
                }
            }
        } else {
            warning(paste(extract(x, what = "bnames", which = w),
                          ": automated plot not reasonable for base-learners",
                          " of matrices", sep=""))
        }
    }
    if (any(foo <- !sapply(RET, is.null))){
        if (sum(foo) == 1)
            return(RET[foo][[1]])
        stop(paste("Cannot plot multiple surfaces via levelplot;",
             "Plot base-learners seperately via plot(model, which = ...)!"))
    }
}


lines.mboost <- function(x, which = NULL, type = "l", rug = FALSE, ...){
    plot(x, which = which, type = type, add = TRUE, rug = rug, ...)
}


plot.glmboost <- function(x, main = deparse(x$call), col = NULL,
                          off2int = FALSE, ...) {

    cp <- coef(x, aggregate = "cumsum", off2int = off2int)
    ncp <- names(cp)
    cp <- matrix(unlist(cp), nrow = length(cp), byrow = TRUE)
    cf <- cp[, ncol(cp)]
    if (is.null(col))
        col <- hcl(h = 40, l = 50, c= abs(cf) / max(abs(cf)) * 490)
    matplot(t(cp), type = "l", lty = 1, xlab = "Number of boosting iterations",
            ylab = "Coefficients", main = main, col = col, ...)
    abline(h = 0, lty = 1, col = "lightgray")
    axis(4, at = cf, labels = ncp, las = 1)
}
