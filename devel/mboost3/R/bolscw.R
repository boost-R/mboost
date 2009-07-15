
bolscw <- function(..., z = NULL, center = FALSE, intercept = TRUE, contrast.arg = NULL) {

    mf <- list(...)
    if (length(mf) == 1 && (isMATRIX(mf[[1]]) || is.data.frame(mf[[1]]))) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    vary <- ""
    if (!is.null(z)) {
        stopifnot(is.data.frame(mf))
        mf <- cbind(mf, z)
        colnames(mf) <- c(colnames(mf), deparse(substitute(z)))
        vary <- colnames(mf)[ncol(mf)]
    }

    index <- NULL
    if (!all(cc <- Complete.cases(mf))) {
        nd <- which(cc)
        index <- match(1:nrow(mf), nd)
        mf <- mf[cc, , drop = FALSE]
    }

    ret <- list(model.frame = function() {
                    if (is.null(index)) return(mf)
                    return(mf[index, , drop = FALSE])
                },
                get_names = function() colnames(mf),
                set_names = function(value) attr(mf, "names") <<- value)
    class(ret) <- "blg"

    newX <- function(newdata) {
        if (isMATRIX(newdata)) return(newdata)
        fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary],
                    collapse = "+"), sep = "")
        if (!intercept)
            fm <- paste(fm, "-1", collapse = "")
        X <- model.matrix(as.formula(fm), data = newdata, contrasts.arg = contrast.arg)
        if (vary != "") {
            ### <FIXME> see X_ols
            z <- model.matrix(as.formula(paste("~", vary, collapse = "")), data = mf)[,2]
            X <- X * z
            ### </FIXME>
        }
        return(X)
    }

    X <- newX(mf)

    ### <FIXME> centering with or without weights?
    if (center) {
        cm <- colSums(X) / nrow(X)
        cls <- sapply(mf, class)[colnames(mf) != vary]
        num <- which(cls == "numeric")
        cm[!attr(X, "assign") %in% num] <- 0
        X <- scale(X, center = cm, scale = FALSE)
    }
    ### </FIXME>

    ret$dpp <- function(weights) {

        weights <- weights[cc]
        xw <- t(X * weights)
        xtx <- colSums(X^2 * weights)
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
                X <- newX(newdata)
                X <- scale(X, center = cm, scale = FALSE)
            }
            return(X %*% cf)
        }
    
        ret <- list(fit = fit, predict = predict)
        class(ret) <- c("bl_lin", "bl")
        return(ret)
    }
    return(ret)
}

coef.bm_cwlin <- function(object, ...) {
    cf <- numeric(object$model[3])
    cf[object$model[2]] <- object$model[1]
    cf
}
