
### just a try
plot.mboost <- function(x, which = NULL, newdata = NULL,
                        type = "b", rug = TRUE, rugcol = "black", ylim = NULL,
                        xlab = NULL, ylab = expression(f[partial]), ...) {

    which <- x$which(which, usedonly = is.null(which))

    pr <- predict(x, which = which)
    if (is.null(ylim)) ylim <- range(pr)

    ## FIXED?
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

    for (w in which) {

        data <- model.frame(x, which = w)[[1]]
        get_vary <- x$baselearner[[w]]$get_vary
        vary <- ""
        if (!is.null(get_vary)) vary <- get_vary()
        if (!is.null(newdata)) data <- newdata[, colnames(data), drop = FALSE]
        if (vary != "") {
            v <- data[[vary]]
            if (is.factor(v)) v <- factor(levels(v)[-1], levels = levels(v))
            if (is.numeric(v)) v <- 1
        }

        plot_helper <- function(xl, yl){
            pr <- predict(x, newdata = data, which = w)
            if (vary != "")
                data <- data[, colnames(data) != vary, drop = FALSE]

            if (ncol(data) == 1) {
                plot(sort(data[[1]]), pr[order(data[[1]])], type = type,
                     xlab = xl, ylab = yl, ylim = ylim, ...)
                if (rug) rug(data[[1]], col = rugcol)
            }
            if (ncol(data) == 2) {
                fm <- as.formula(paste("pr ~ ", paste(colnames(data), collapse = "*"), sep = ""))
                print(levelplot(fm, data = data, ...))
            }
            if (ncol(data) > 2) {
                for (v in colnames(data)) {
                    tmp <- data
                    pardep <- sapply(data[[v]], function(vv) {
                        tmp[[v]] <- vv
                        mean(predict(x, newdata = tmp, which = w))
                    })
                    plot(sort(data[[v]]), pardep[order(data[[v]])], type = type,
                         xlab = v, ylab = "Partial Dependency", ylim = ylim, ...)
                }
            }
        } # END plot_helper()

        xl <- ifelse(length(xlab) > 1, xlab[w], xlab[1])
        yl <- ifelse(length(ylab) > 1, ylab[w], ylab[1])

        if (!isMATRIX(data)){
            if (vary == "") plot_helper(xl, yl)
            if (vary != ""){
                for (i in 1:length(v)){
                    data[[vary]] <- v[rep(i, nrow(data))]
                    if (!userspec){
                        ## xlab not user specified
                        plot_helper(paste(xl, "=", v[i]), yl)
                    } else {
                        plot_helper(paste(xl, "(", vary, "=", v[i], ")"), yl)
                    }
                }
            }
        } else {
            warning(paste(extract(x, what = "bnames", which = w),
                          ": automated plot not reasonable for base-learners",
                          " of matrices", sep=""))
        }
    }
}

plot.glmboost <- function(x, main = deparse(x$call), col = NULL, ...) {

    cp <- coef(x, aggregate = "cumsum")
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
