varimp_mboost <- function(object, percent = FALSE) {
  
  ### check arguments
  if( !(inherits(object, "gamboost") || inherits(object, "glmboost")) )
    stop(paste(deparse(substitute(object)), 
      "has to be a gamboost or glmboost object."))
  
  ### which baselearners were selected in boosting steps
  # differentiation between gamboost- and glmboost-object
  blearner_names <- if( inherits(object, "gamboost") ) 
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
  
  ### compute empirical risk (according to output in cvrisk)
  riskdiff <- riskdiff / length(object$response)
  
  ### explained Risk attributed to baselearners
  explained <- sapply(seq_along(blearner_names), FUN = function(i) {
    sum( riskdiff[which(blearner_selected == i)] ) 
  })
  names(explained) <- blearner_names
  if( percent ) { explained <- explained / sum(explained) }
  
  class(explained) <- "varimp_mboost"
  attr(explained, "percent") <- percent
  return(explained)
}


plot.varimp_mboost <- function(x, nbars = 20L, maxchar = 20L, xlim, ...) {
  args <- as.list(match.call())
  
  ### check arguments
  if( !(class(x) == "varimp_mboost") ) 
    stop(paste(deparse(substitute(x)), "is not of class varimp_mboost."))
  
  if( !(is.numeric(nbars) && nbars > 0 && length(nbars) == 1) )
    stop("Parameter nbars has to be a positive integer.")
  
  if( !(is.numeric(maxchar) && maxchar > 0 && length(maxchar) == 1) )
    stop("Parameter maxchar has to be a positive integer.")
  
  if( hasArg(xlab) || hasArg(ylab) )
    stop("xlab, ylab already defined by default.")  
  
  
  ### set x-axis label
  if( attr(x, "percent") ) xlab = "Rel. Empirical Risk Reduction" else 
    xlab = "Empirical Risk Reduction"
  
  ### set range for x-axis
  if( !("xlim" %in% names(args)) ) {
    if( attr(x, "percent") ) xlim = c(0,1)
    else xlim = c(0,max(x))
  }  
  
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
  
  # adjust left margin to length of horizontal y labels
  leftmargin <- max(strwidth(names(xsorted), "inch")+.4, na.rm = TRUE)
  opar <- par(mai=c(1.02, leftmargin, 0.82, 0.42))
  
  # plot risk reduction per variable 
  barplot(height = xsorted, horiz = TRUE, las = 1,
          xlab = "Empirical Risk Reduction", #ylab = "Baselearner",
          xlim = xlim, ...)
  box()
  par(opar)
}
