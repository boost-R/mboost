varimp <- function(object, ...)
  UseMethod("varimp")

varimp.mboost <- function(object) {
  
  ### check arguments
  if( !(inherits(object, "mboost")) )
    stop(paste(deparse(substitute(object)), "has to be a mboost object."))
  
  ### which baselearners were selected in boosting steps
  # differentiation between gamboost- and glmboost-object
  blearner_names <- if( inherits(object, "mboost") && !(inherits(object, "glmboost")) ) 
    names(object$baselearner) else names(object$baselearner[[1]])  
  blearner_selected <- object$xselect()
   
  ### compute risks for each step
  # initial risk for the intercept model
  y <- object$response
  if(is.factor(y)) y <- 2*as.numeric(y) - 3
  risk0 <- object$family@risk( y = y, f = object$offset ) 
  # risk after each boosting-steps
  riskstep <- object$risk()
  # risk reduction per step
  riskdiff <- c(risk0, riskstep[-length(riskstep)]) - riskstep  
  
  ### compute empirical inbag risk (according to output in cvrisk)
  riskdiff <- riskdiff / length(object$response)
  
  ### explained Risk attributed to baselearners
  explained <- sapply(seq_along(blearner_names), FUN = function(i) {
    sum( riskdiff[which(blearner_selected == i)] ) 
  })
  names(explained) <- blearner_names
  
  class(explained) <- "varimp"
  attr(explained, "selprobs") <- sapply(seq_along(blearner_names), function(i) mean(blearner_selected == i))
  
  # add variable names per baselearner
  variable_names = sapply(names(object$baselearner), 
    function(x) names(object$baselearner[[x]]))
  
  attr(explained, "variables") <- unlist( lapply(variable_names, function(x) {
    do.call(function(...) paste(..., sep = ","), as.list(x)) 
  }) )
  
  return(explained)
}


plot.varimp <- function(x, percent = TRUE, type = "blearner", 
  nbars = 20L, maxchar = 20L, xlim, ...) {
  
  args <- as.list(match.call())
  
  ### check arguments
  if( !(class(x) == "varimp") ) 
    stop(paste(deparse(substitute(x)), "is not of class varimp."))
  
  if( !(is.numeric(nbars) && nbars > 0 && length(nbars) == 1) )
    stop("Parameter nbars has to be a positive integer.")
  
  if( !(is.numeric(maxchar) && maxchar > 0 && length(maxchar) == 1) )
    stop("Parameter maxchar has to be a positive integer.")
  
  if( hasArg(xlab) || hasArg(ylab) )
    stop("xlab, ylab already defined by default.")  
  
  ### set x-axis label
  if( percent ) xlab = "In-bag Risk Reduction (%)" else 
    xlab = "In-bag Risk Reduction"
  
  ### set range for x-axis
  if( !("xlim" %in% names(args)) ) {
    if( percent ) xlim = c(0,1)
    else xlim = c(0,sum(x))
  }  
  
  ### transform to pecent if corresponding argument is true
  if( percent ) { x <- x / sum(x) }
  
  
#   # add up values for baselearners that are based on same variable
#   # e.g. bols(var1) and bbs(var1)
#   if( type == "variable" ) {
#     
#     variable_order = order(unlist( lapply(unique(attr(x, "variables")), function(i) {
#       sum(x[which(attr(x, "variables") == i)]) 
#     }) ))
#     
#     plot_data = data.frame(reduction = as.numeric(x), 
#                            blearner = factor(names(x)),
#                            varname = ordered(attr(x, "variables"), 
#       levels = unique(attr(x, "variables")[variable_order])) )
#     
#     barchart(varname ~ reduction, groups = blearner, stack = TRUE,
#              horizontal = TRUE, data = plot_data, auto.key = list())
#   }
  
  
  ### sort values
  xsorted <- sort(x)
  # if number of baselearners exceeds nbars, aggregate values of other blearners
  if( length(xsorted) > nbars ) {
    xsorted <- c( sum( head(xsorted, length(xsorted)-nbars+1) ), 
      tail(xsorted, nbars-1) ) 
    names(xsorted)[1] <- "other"
  }
  
  ### create labels for bars
  names(xsorted) <- sapply(names(xsorted), FUN = function(name) {
    paste(strtrim(name, maxchar), if( nchar(name) < maxchar ) "" else "..") })
  
  # add selection probabilities to bar labels
  selprob_char <- as.character(attr(x,"selprob")[order(x, decreasing = FALSE)])
  names(xsorted) = paste(names(xsorted), "\n" , paste0("sel. prob: ~", selprob_char))
  
  barchart(x = xsorted, horizontal = TRUE, xlab = xlab, ylab = "Baselearner", 
    xlim = xlim, ...)
}
