
hatMatTH <- function(x, w = NULL, df = 4) {
    n <- NROW(x)
    indx <- diag(n)
    x <- signif(x, 10)
    apply(indx, 2, function(y)
        predict(smoothbase(x = x, ux = unique(sort(x)), y = y, w = w, df = df),
                x = x)$y)
}

complete_cases <- function(x, y = NULL, z = NULL) {

    tmp <- list(x = x, y = y, z = z)
    tmp <- tmp[!sapply(tmp, is.null)]
    rowSums(sapply(tmp, is.na)) == 0
}

predict.baselist <- function(object, ...) {

    pr <- predict(object[[1]], ...)
    if (length(object) == 1) return(pr)
    for (i in 2:length(object)) {
        if (!any(is.na(pr)))
            break
        tmp <- predict(object[[i]], ...)
        pr[is.na(pr)] <- tmp[is.na(pr)]
    }
    return(pr)
}

fitted.baselist <- function(object) {

    pr <- fitted(object[[1]])
    if (length(object) == 1) return(pr)
    for (i in 2:length(object)) {
        if (!any(is.na(pr)))
            break
        tmp <- fitted(object[[i]])
        pr[is.na(pr)] <- tmp[is.na(pr)]
    }
    return(pr)
}

### what happens to weights
### when calculating knots etc?
bbs <- function(x, z = NULL, df = 4, knots = NULL, degree = 3, differences = 2,
                center = FALSE, xname = NULL, zname = NULL) {

    cc <- complete_cases(x = x, z = z)

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (is.factor(x) || (df <= 2 && !center)) 
        return(bols(x = x, z = z, xname = xname, zname = zname))

    if (differences < 1 || differences > 3) 
        stop(sQuote("differences"), " are not in 1:3")
    if ((!center) && (df < differences))
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")
    
    n.kn <- function(n) {
        ## Number of inner knots
        if(n < 50) n
        else trunc({
            a1 <- log( 50, 2)
	    a2 <- log(100, 2)
            a3 <- log(140, 2)
	    a4 <- log(200, 2)
	    if	(n < 200) 2^(a1+(a2-a1)*(n-50)/150)
	    else if (n < 800) 2^(a2+(a3-a2)*(n-200)/600)
	    else if (n < 3200)2^(a3+(a4-a3)*(n-800)/2400)
	    else  200 + (n-3200)^0.2
        })
    }
    
    if(is.null(knots)) {
        n <- length(x)
        nk <- n.kn(n)    
        knots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = nk)
        knots <- knots[2:(length(knots) - 1)]   
    }

    if (length(knots) == 1) {
        knots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = knots+2) 
        knots <- knots[2:(length(knots) - 1)]
    }
    boundary.knots <- range(x, na.rm = TRUE)
    newX <- function(x, z = NULL, na.rm = TRUE) {
        if (na.rm) {
            x <- x[cc]
            if (!is.null(z))
                z <- z[cc]
        }
        X <- bs(x, knots = knots, degree = degree, intercept = TRUE, 
                Boundary.knots = boundary.knots)
        if (!is.null(z))
            X <- X * z
        if (center) {
            K <- diff(diag(ncol(X)), differences = differences)
            X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
        }
        return(X)
    }
    X <- newX(x, z)
    Xna <- X
    if (any(!cc))
        Xna <- newX(x, z, na.rm = FALSE)

    if (center) {
        K <- diag(ncol(X))
    } else {
        K <- diff(diag(ncol(X)), differences = differences)
        K <- crossprod(K, K) 
    }

    dpp <- function(weights) {

        if (any(!cc)) weights <- weights[cc]
        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            
            if (any(!cc)) y <- y[cc]
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(Xna %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]], na.rm = FALSE)
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = Xna %*% coef)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret	
    }
    attr(X, "dpp") <- dpp
    return(X)
}

predict.basefit <- function(object, newdata = NULL)
    object$predict(newdata)

fitted.basefit <- function(object)
    object$fitted

df2lambda <- function(X, df = 4, dmat = NULL, weights) {

#    if (df <= 2) stop(sQuote("df"), " must be greater than two")

    if (is.null(dmat)) {
        dmat <- diff(diag(ncol(X)), differences = 2)
        dmat <- crossprod(dmat, dmat)
    }

    # singular value decomposition
    A <- crossprod(X * weights, X)
    decomp <- svd(A+0.0001)
    # A is equal to decomp$u %*% diag(decomp$d) %*% t(decomp$v)
    u <- decomp$u
    v <- decomp$v
    d <- decomp$d
    K <- crossprod(u,dmat)%*%v

    # df2lambda
    dd <- diag(d)
    df2l <- function(lambda)
        (sum(diag(d * solve( (dd+lambda*K) ) ))-df)^2

    lower.l <- 0
    upper.l <- 5000
    lambda <- upper.l

    while (lambda >= upper.l - 200 ) {
        upper.l <- upper.l * 1.5
        lambda <- optimize(df2l, interval=c(lower.l,upper.l))$minimum
        lower.l <- upper.l-200
        if (lower.l > 1e+06)
            stop("Maximum size of lambda reached")
    }

    tmp <- sum(diag(X %*% solve(crossprod(X * weights, X) + 
                    lambda*dmat) %*% t(X * weights))) - df
    if (abs(tmp) > sqrt(.Machine$double.eps)) 
        warning("trace of hat matrix is not equal df with difference", tmp)

    lambda
}

bns <- function(x, z = NULL, df = 4, knots = NULL, differences = 2,
                xname = NULL, zname = NULL) {

    if (is.null(xname)) xname <- deparse(substitute(x))
    if (is.null(zname)) zname <- deparse(substitute(z))

    if (is.factor(x) || df <= 2) 
        return(bols(x = x, z = z, xname = xname, zname = zname))

    if (differences < 1 || differences > 3) 
        stop(sQuote("differences"), " are not in 1:3")
    if (df < differences)
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")

    
    n.kn <- function(n) {
        ## Number of inner knots
        if(n < 50) n
        else trunc({
            a1 <- log( 50, 2)
	    a2 <- log(100, 2)
            a3 <- log(140, 2)
	    a4 <- log(200, 2)
	    if	(n < 200) 2^(a1+(a2-a1)*(n-50)/150)
	    else if (n < 800) 2^(a2+(a3-a2)*(n-200)/600)
	    else if (n < 3200)2^(a3+(a4-a3)*(n-800)/2400)
	    else  200 + (n-3200)^0.2
        })
    }
        
    
    if(is.null(knots)) {
        n <- length(x)
        nk <- n.kn(n)    
        knots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = nk)
        knots <- knots[2:(length(knots) - 1)]   
    }
    
    
    
    if (length(knots) == 1) {
        knots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = knots + 2) 
        knots <- knots[2:(length(knots) - 1)]
        #knots <- c(min(x)-sd(x),knots,max(x)+sd(x))
    }
    newX <- function(x, z = NULL) {
        X <- ns(x, knots = knots, intercept = TRUE, Boundary.knots = c(min(x)-sd(x),max(x)+sd(x)) )
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, z)

    K <- diff(diag(ncol(X)), differences = differences)
    K <- crossprod(K, K) 

    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]])
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = X %*% coef)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

bss <- function(x, df = 4, xname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.factor(x) || df <= 2) 
        return(bols(x = x, xname = xname))

    xs <- signif(x, 10)   
    ux <- unique(sort(xs))

    dpp <- function(weights) {

        fitfun <- function(y) {
            object <- smoothbase(x = xs, ux = ux, y = y, w = weights, df = df)

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(object$yfit)
                stats:::predict.smooth.spline.fit(object, x = newdata[[xname]])$y
            }

            ret <- list(basemodel = object, predict = predictfun, 
                        fitted = object$yfit)
            class(ret) <- "basefit"
            ret
        }  

        ret <- list(fit = fitfun, hatmatrix = function() 
                                      hatMatTH(x = x, w = weights, df = df))
        class(ret) <- "basisdpp"
        ret
    }
    x <- matrix(x, nc = 1)
    attr(x, "dpp") <- dpp
    return(x)
}

bspatial <- function(x, y, z = NULL, df = 5, xknots = NULL, yknots = NULL, 
                     degree = 3, differences = 2, center = FALSE, xname = NULL,
                     yname = NULL, zname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(yname)) yname = deparse(substitute(y))
    if (is.null(zname)) zname = deparse(substitute(z))

#    if (df <= 2) stop(sQuote("df"), " must be greater two")
    if (differences < 1 || differences > 3) 
        stop(sQuote("differences"), " are not in 1:3")
    if ((!center) && (df < differences^2))
        stop(sQuote("df"), " is less than ", sQuote("differences^2"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")
    if (length(unique(y)) < 6)
        stop(sQuote(yname), " has less than 6 unique values")
        
        
    n.kn <- function(n) {
        ## Number of inner knots
        if(n < 50) n
        else trunc({
            a1 <- log( 50, 2)
	    a2 <- log(100, 2)
            a3 <- log(140, 2)
	    a4 <- log(200, 2)
	    if	(n < 200) 2^(a1+(a2-a1)*(n-50)/150)
	    else if (n < 800) 2^(a2+(a3-a2)*(n-200)/600)
	    else if (n < 3200)2^(a3+(a4-a3)*(n-800)/2400)
	    else  200 + (n-3200)^0.2
        })
    }
        
    
    if(is.null(xknots)) {
        n <- length(x)
        nk <- n.kn(n)    
        xknots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = nk + 2)
        xknots <- xknots[2:(length(xknots) - 1)]   
    }    
    if(is.null(yknots)) {
        n <- length(y)
        nk <- n.kn(n)    
        yknots <- seq(from = min(y, na.rm = TRUE), 
                     to = max(y, na.rm = TRUE), length = nk + 2)
        yknots <- yknots[2:(length(yknots) - 1)]   
    }        
        

    if (length(xknots) == 1) {
        xknots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = xknots + 2) 
        xknots <- xknots[2:(length(xknots) - 1)]
    }
    if (length(yknots) == 1) {
        yknots <- seq(from = min(y, na.rm = TRUE), 
                     to = max(y, na.rm = TRUE), length = yknots + 2) 
        yknots <- yknots[2:(length(yknots) - 1)]
    }
    newX <- function(x, y, z) {
        Xx <- bs(x, knots = xknots, degree = degree, intercept = TRUE)
        Xy <- bs(y, knots = yknots, degree = degree, intercept = TRUE)
        X <- kronecker(Xx, matrix(1, nc = ncol(Xy))) * kronecker(matrix(1, nc = ncol(Xx)), Xy)
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, y, z)

    xd <- length(xknots) + degree + 1
    yd <- length(yknots) + degree + 1

    Kx <- diff(diag(xd), differences = differences)
    Kx <- crossprod(Kx, Kx) 
    Ky <- diff(diag(yd), differences = differences)
    Ky <- crossprod(Ky, Ky) 
    K <- kronecker(Kx, diag(yd)) + kronecker(diag(xd), Ky)

    L <- 0
    if(center) {
        L <- eigen(K, symmetric=TRUE, EISPACK=TRUE)
        L$vectors <- L$vectors[,1:(ncol(X)-differences^2)]
        L$values <- sqrt(L$values[1:(ncol(X)-differences^2)])
        L <- L$vectors%*%diag(1/L$values)
        X <- X%*%L
        K <- diag(ncol(X))
    }
 
    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)
        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], y = newdata[[yname]],
                           z = newdata[[zname]])
                if(center) {
                    nX <- nX%*%L
                }
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, 
                        fitted = X %*% coef)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

bols <- function(x, z = NULL, xname = NULL, zname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(zname)) zname = deparse(substitute(z))

    cc <- complete_cases(x = x, z = z)

    newX <- function(x, z = NULL, na.rm = TRUE) {
        if (na.rm) {
            x <- x[cc]
            if (!is.null(z))
                z <- z[cc]
        }
        X <- model.matrix(~ x)
        if (any(!cc) & !na.rm) {
            Xtmp <- matrix(NA, ncol = ncol(X), nrow = length(cc))
            Xtmp[cc,] <- X
            X <- Xtmp
        }
        if (!is.null(z)) X <- X * z
        X
    }
    X <- newX(x, z)
    Xna <- X
    if (any(!cc))
        Xna <- newX(x, z, na.rm = FALSE)

    dpp <- function(weights) {

        if (any(!cc)) weights <- weights[cc]
        Xw <- X * weights
        Xsolve <- tcrossprod(solve(crossprod(Xw, X)), Xw)

        fitfun <- function(y) {
           
            if (any(!cc)) y <- y[cc]
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(Xna %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]], 
                           na.rm = FALSE)
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, 
                        fitted = Xna %*% coef)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}

brandom <- function(x, z = NULL, df = 4, xname = NULL, 
                    zname = NULL) {

    if (is.null(xname)) xname = deparse(substitute(x))
    if (is.null(zname)) zname = deparse(substitute(z))

    newX <- function(x, z = NULL) {
        if (!is.factor(x)) stop(sQuote("x"), " is not a factor")
        X <- model.matrix(~ x - 1)
        if (!is.null(z))
            X <- X * z
        return(X)
    }
    X <- newX(x, z)

    K <- diag(ncol(X))

    dpp <- function(weights) {

        lambda <- df2lambda(X, df = df, dmat = K, weights = weights)

        Xw <- X * weights
        XtX <- crossprod(Xw, X)
        ### XtX and K are diagonal matrices
        Xsolve <- tcrossprod(solve(XtX + lambda * K), Xw)

        fitfun <- function(y) {
            coef <- Xsolve %*% y

            predictfun <- function(newdata = NULL) {
                if (is.null(newdata)) return(X %*% coef)
                nX <- newX(x = newdata[[xname]], z = newdata[[zname]])
                nX %*% coef
            }
            ret <- list(model = coef, predict = predictfun, fitted = X %*% coef)
            class(ret) <- "basefit"
            ret
        }
        ret <- list(fit = fitfun, hatmatrix = function() X %*% Xsolve)
        class(ret) <- "basisdpp"
        ret
    }
    attr(X, "dpp") <- dpp
    return(X)
}
