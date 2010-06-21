
### component-wise linear models baselearner
bolscw <- function(X) {

    ### formula parsing, data handling is done elsewhere
    stopifnot(isMATRIX(X))

    ### deal with missing values
    index <- NULL
    if (!all(cc <- Complete.cases(X))) {
        nd <- which(cc)
        index <- match(1:nrow(X), nd)
        X <- X[cc, , drop = FALSE]
    }

    ### we actually don't have raw data available
    ret <- list(model.frame = function() return(NULL),
                get_data = function() X,
                get_index = function() index,
                get_vary = function() NULL,
                get_names = function() colnames(X),
                set_names = function(value) attr(X, "colnames") <<- value)
    class(ret) <- "blg"

    ret$dpp <- function(weights) {

        weights <- weights[cc]
        xtx <- colSums(X^2 * weights, na.rm = TRUE)
        ### some columns may be zero
        xtx[xtx < .Machine$double.eps] <- 1
        sxtx <- sqrt(xtx)
        p <- ncol(X)
        ### use crossprod for estimation but not for Matrix objects
        ### since S4 dispatch is way to slow
        if (is(X, "Matrix")) {
            MPinvS <- t(X * weights) / sxtx
            est <- function(y) MPinvS %*% y
        } else {
            MPinvS <- NULL
            storage.mode(X) <- "double" ### crossprod is doing this anyway
            est <- function(y) crossprod(X, y * weights) / sxtx
        }

        fit <- function(y) {
            if (!is.null(index))
                y <- y[cc]

            mmu <- max(amu <- abs(mu <- est(y)))
            xselect <- which(as.logical(mmu == amu))[1]
            coef <- mu[xselect] / sxtx[xselect]
            ret <- list(model = c(coef = coef, xselect = xselect, p = p),
                        fitted = function() {
                            if (is.null(index)) return(coef * X[,xselect, drop = FALSE])
                            return(coef * X[index, xselect,drop = FALSE])
                        })
            class(ret) <- c("bm_cwlin", "bm_lin", "bm")
            return(ret)
        }

        predict <- function(bm, newdata = NULL, 
                            aggregate = c("sum", "cumsum", "none")) {

            aggregate <- match.arg(aggregate)
            cf <- switch(aggregate, "sum" = {
                cf <- rep(0, bm[[1]]$model["p"])
                for (i in 1:length(bm)) {
                    m <- bm[[i]]$model
                    cf[m[2]] <- cf[m[2]] + m[1]
                }
                cf
            },
            "cumsum" = {
                cf <- matrix(coef(bm[[1]], all = TRUE), ncol = 1)
                if (length(bm) > 2) {
                    for (i in 2:length(bm))
                        cf <- cbind(cf, cf[[i-1]] + coef(bm[[i]], all = TRUE))
                }
                cf
            },
            "none" = {
                cf <- matrix(coef(bm[[1]], all = TRUE), ncol = 1)
                if (length(bm) > 2) {
                    for (i in 2:length(bm))
                        cf <- cbind(cf, coef(bm[[i]], all = TRUE))
                }
                cf
            })

            if (!is.null(newdata)) {
                stopifnot(all(class(newdata) == class(X)))
                stopifnot(all(colnames(newdata) == colnames(X)))
                X <- newdata
                return(X %*% cf)
            }
            if (!is.null(index))
                return(X[index,,drop = FALSE] %*% cf)
            return(X %*% cf)
        }

        ret <- list(fit = fit, predict = predict, Xnames = colnames(X), 
                    MPinv = function() {
                        if (is.null(MPinvS)) MPinvS <<- t(X * weights) / sxtx
                        return(MPinvS / sxtx)
                    },
                    hatvalues = function() {
                        if (is.null(MPinvS)) MPinvS <<- t(X * weights) / sxtx
                        X %*% (MPinvS / sxtx)
                    })
        class(ret) <- c("bl_cwlin", "bl")
        return(ret)
    }
    return(ret)
}

coef.bm_cwlin <- function(object, all = FALSE, ...) {
    if (!all) return(object$model[1])
    cf <- numeric(object$model[3])
    cf[object$model[2]] <- object$model[1]
    cf
}

### extract hatmatrix
hatvalues.bl_cwlin <- function(model)
    model$hatvalues()
