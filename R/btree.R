
### the classical tree-based baselearner; stumps by default
### (also fits an additive model)
btree <- function(..., by = NULL, nmax = Inf,
    tree_controls = partykit::ctree_control(
        teststat = "quad",
        testtype = "Teststatistic",
#        splittest = TRUE,
        mincriterion = 0,
        minsplit = 10,
        minbucket = 4,
        maxdepth = 1, 
        saveinfo = FALSE))
{

    cll <- match.call()
    cll[[1]] <- as.name("btree")

    ctrl <- tree_controls
    mf <- list(...)
    if (length(mf) == 1 && is.data.frame(mf[[1]])) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }

    vary <- ""
    if (!is.null(by)){
        vary <- deparse(substitute(by))
        stopifnot(is.factor(by) || all(by %in% 0:1))
        if (is.factor(by)) {
            stopifnot(nlevels(by) == 2L)
            xlev <- levels(by)
            iby <- (0:1)[by]
        } else {
            iby <- by
        }
    } else {
        iby <- rep(1, NROW(mf))
    }

    mffun <- function() {
        if (vary == "") return(mf)
        ret <- cbind(mf, by)
        colnames(ret)[ncol(ret)] <- vary
        ret
    }

    ret <- list(model.frame = mffun,
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_names = function() colnames(mf),
                get_vary = function() vary,
                set_names = function(value) {
                    if(length(value) != length(colnames(mf)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf, "names") <<- value
                })
    class(ret) <- "blg"


    ret$dpp <- function(weights) {

        ### construct design matrix etc.
        y <- vector(length = nrow(mf), mode = "numeric")
        ### name for working response (different from any x)
        rname <- paste(sample(LETTERS, 25, replace = TRUE), collapse = "")
        fm <- as.formula(paste(rname, " ~ ", paste(colnames(mf), collapse = "+")))
        df <- mffun()
        df[[rname]] <- y
        d <- extree_data(fm, data = df, yx = "none", 
                         nmax = c(yx = Inf, z = nmax))
        Y <- NULL
        ytrafo <- function(subset, weights, info, estfun, object, ...) 
            list(estfun = Y, unweighted = TRUE) 
        mymf <- model.frame(d)
        tree_controls$update <- FALSE
        subset <- which(weights > 0 & iby == 1)
        weights <- weights * iby

        fitfun <- function(y) {
            if (!is.matrix(y)) y <- matrix(y, ncol = 1)
 
            assign("Y", y, envir = environment(ytrafo))
            tree <- extree_fit(data = d, trafo = ytrafo, converged = TRUE, 
                               partyvars = d$variables$z, subset = subset, 
                               weights = weights, ctrl = tree_controls, 
                               doFit = TRUE)$node
            where <- factor(fitted_node(tree, mymf))
            coef <- do.call("rbind", tapply(1:NROW(y), where, function(i)
                colSums(y[i,,drop = FALSE] * weights[i]) / sum(weights[i]), 
                simplify = FALSE))

            fitted <- function()
                return(coef[unclass(where),,drop = FALSE] * iby)

            predict <- function(newdata = NULL) {
                if (is.null(newdata)) newdata <- mffun()
                vmatch <- match(names(mymf), names(newdata))
                wh <- factor(fitted_node(tree, newdata, vmatch = vmatch), 
                             levels = levels(where), labels = levels(where))
                ret <- coef[unclass(wh),,drop = FALSE]
                if (vary != "") {
                    by <- newdata[, vary]
                    stopifnot(is.factor(by) || all(by %in% 0:1))
                    if (is.factor(by)) {
                        stopifnot(isTRUE(all.equal(levels(by), xlev)))
                        by <- (0:1)[by]
                    }
                    ret <- ret * by
                }
                return(ret)
            }

            ret <- list(model = tree, fitted = fitted, predict = predict)
            class(ret) <- c("bm_tree", "bm")
            ret
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            aggregate <- match.arg(aggregate)

            if (is.null(newdata)) 
                newdata <- mffun()
            if (!is.data.frame(newdata))
                newdata <- as.data.frame(newdata)

            pr <- 0
            for (i in 1:length(bm)) {
                pri <- bm[[i]]$predict(newdata)
                if (aggregate == "sum") {
                    pr <- pr + pri
                } else {
                    if (i > 1) {
                        pr <- cbind(pr, pri)
                    } else {
                        pr <- pri
                    }
                    if (aggregate == "cumsum")
                        if (i > 1) pr[,i] <- pr[,i] + pr[,i-1]
                }
            }
            return(pr)
        }

        ret <- list(fit = fitfun, predict = predict)
        class(ret) <- c("bl_tree", "bl")
        ret
    }
    return(ret)
}
