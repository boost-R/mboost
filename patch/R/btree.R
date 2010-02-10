
### the classical tree-based baselearner; stumps by default
### (also fits an additive model)
btree <- function(..., tree_controls = ctree_control(stump = TRUE, mincriterion = 0)) {

    if (!require("party"))
        stop("cannot load ", sQuote("party"))

    cll <- match.call()
    cll[[1]] <- as.name("btree")
    cll <- deparse(cll, width.cutoff=500L)
    if (length(cll) > 1)
        cll <- paste(cll, collapse="")

    ctrl <- tree_controls
    mf <- list(...)
    if (length(mf) == 1 && is.data.frame(mf[[1]])) {
        mf <- mf[[1]]
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }

    ret <- list(model.frame = function() return(mf),
                get_call = function() cll,
                get_names = function() colnames(mf),
                set_names = function(value) {
                    attr(mf, "names") <<- value
                    cll <<- paste("btree", "(", paste(colnames(mf),
                        collapse = ", "), ")", sep = "")
                })
    class(ret) <- "blg"


    ret$dpp <- function(weights) {

        ### construct design matrix etc.
        y <- vector(length = nrow(mf), mode = "numeric")
        ### name for working response (different from any x)
        rname <- paste(sample(LETTERS, 25, replace = TRUE), collapse = "")
        fm <- as.formula(paste(rname, " ~ ", paste(colnames(mf), collapse = "+")))
        df <- mf
        df[[rname]] <- y
        object <- party:::ctreedpp(fm, data = df)
        fitmem <- ctree_memory(object, TRUE)
        where <- rep.int(0, nrow(mf))
        storage.mode(where) <- "integer"
        storage.mode(weights) <- "double"

        fitfun <- function(y) {

            .Call("R_modify_response", as.double(y), object@responses,
                 PACKAGE = "party")
            tree <- .Call("R_TreeGrow", object, weights, fitmem, ctrl,
                          where, PACKAGE = "party")
            .Call("R_remove_weights", tree, package = "party")

            fitted <- function() {
                wh <- .Call("R_get_nodeID", tree, object@inputs, 0.0, PACKAGE = "party")
                return(unlist(.Call("R_getpredictions", tree, wh, PACKAGE = "party")))
            }

            ret <- list(model = tree, fitted = fitted)
            class(ret) <- c("bm_tree", "bm")
            ret
        }

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            aggregate <- match.arg(aggregate)

            if (is.null(newdata)) {
                newinp <- object@inputs
            } else {
                newinp <- party:::newinputs(object, newdata)
            }

            pr <- 0
            for (i in 1:length(bm)) {
                wh <- .Call("R_get_nodeID", bm[[i]]$model, newinp, 0.0,
                         PACKAGE = "party")
                pri <- unlist(.Call("R_getpredictions", bm[[i]]$model, wh, PACKAGE = "party"))
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
