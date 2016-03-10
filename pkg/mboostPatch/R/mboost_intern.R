mboost_intern <- function(...,
    fun = c("df2lambda", "hyper_bbs", "hyper_ols",
            "bl_lin", "bl_lin_matrix",
            "Complete.cases", "get_index", "isMATRIX",
            "cbs", "bsplines", "model.frame.blg", 
            "check_newdata", "do_trace", "rescale_weights", "nnls2D")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
