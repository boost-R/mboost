
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
        xw <- t(X * weights)
        xtx <- colSums(X^2 * weights, na.rm = TRUE)
        ### some columns may be zero
        xtx[xtx < .Machine$double.eps] <- 1
        sxtx <- sqrt(xtx)
        MPinvS <- (1 / sxtx) * xw
        p <- ncol(X)

        fit <- function(y) {
            if (!is.null(index))
                y <- y[cc]

            mmu <- max(amu <- abs(mu <- MPinvS %*% y))
            xselect <- which(as.logical(mmu == amu))[1]
            coef <- mu[xselect] / sxtx[xselect]
            ret <- list(model = c(coef, xselect, p),
                        fitted = function() {
                            if (is.null(index)) return(coef * X[,xselect])
                            return(coef * X[index, xselect])
                        })
            class(ret) <- c("bm_cwlin", "bm_lin", "bm")
            return(ret)
        }

        predict <- function(bm, newdata = NULL, 
                            aggregate = c("sum", "cumsum", "none")) {

            aggregate <- match.arg(aggregate)
            cf <- switch(aggregate, "sum" = {
                cf <- rep(0, bm[[1]]$model[3])
                for (i in 1:length(bm)) {
                    m <- bm[[i]]$model
                    cf[m[2]] <- cf[m[2]] + m[1]
                }
                cf
            },
            "cumsum" = {
                cf <- coef(bm[[1]])
                for (i in 2:length(bm))
                    cf <- cbind(cf, cf[[i-1]] + coef(bm[[i]]))
                cf
            },
            "none" = {
                cf <- coef(bm[[1]])
                for (i in 2:length(bm))
                    cf <- cbind(cf, coef(bm[[i]]))
                cf
            })

            if (!is.null(newdata)) {
                stopifnot(all(class(newdata) == class(X)))
                stopifnot(all(colnames(newdata) == colnames(X)))
                X <- newdata
            }
            return(X %*% cf)
        }

        ret <- list(fit = fit, predict = predict, Xnames = colnames(X), 
                    MPinv = function() MPinvS / sxtx,
                    hatvalues = function() X %*% (MPinvS / sxtx))
        class(ret) <- c("bl_cwlin", "bl")
        return(ret)
    }
    return(ret)
}

coef.bm_cwlin <- function(object, ...) {
    cf <- numeric(object$model[3])
    cf[object$model[2]] <- object$model[1]
    cf
}

### extract hatmatrix
hatvalues.bl_cwlin <- function(model)
    model$hatvalues()
