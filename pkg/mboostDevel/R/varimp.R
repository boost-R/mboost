varimp <- function(object, ...)
  UseMethod("varimp")


varimp.mboost <- function(object) {  

  ### which baselearners/variables were selected in boosting steps
  learner_names <- if( !(inherits(object, "glmboost")) ) 
    names(object$baselearner) else variable.names(object)  
  learner_selected <- object$xselect()
  
  
  ### --------------------------------------------------
  ### compute risks for each step
  # initial risk for the intercept model
  y <- object$response
  if(is.factor(y)) y <- 2*as.numeric(y) - 3
  risk0 <- object$family@risk( y = y, f = object$offset ) 
  # risk after each boosting-steps
  riskstep <- object$risk()
  # risk reduction per step
  riskdiff <- c(risk0, riskstep[-length(riskstep)]) - riskstep  
  
  ### compute empirical in-bag risk (according to output in cvrisk)
  riskdiff <- riskdiff / length(object$response)
  
  
  ### --------------------------------------------------
  ### explained risk attributed to baselearners
  explained <- sapply(seq_along(learner_names), FUN = function(i) {
    sum( riskdiff[which(learner_selected == i)] ) 
  })
  names(explained) <- learner_names
  
  
  ### --------------------------------------------------
  ### define new varimp-object
  class(explained) <- "varimp"
  
  # add selection probabilities to varimp-object
  attr(explained, "selprobs") <- sapply(seq_along(learner_names), function(i) {
    mean(learner_selected == i) 
  })
  
  # add variable names per baselearner
  attr(explained, "variable_names") <- unlist( 
    lapply(variable.names(object) , function(x) {
      do.call(function(...) paste(..., sep = ","), as.list(x)) 
    }) 
  )  
  return(explained)
}


plot.varimp <- function(x, percent = TRUE, type = "blearner", 
  nbars = 20L, maxchar = 20L, xlim, ...) {
  
  args <- as.list(match.call())
  
  
  ### --------------------------------------------------
  ### check arguments
  if( !(class(x) == "varimp") ) 
    stop(paste(deparse(substitute(x)), "is not of class varimp."))
  
  if( !(is.logical(percent)) )
     stop("Parameter percent has of type logical.")
  
  if( !(type %in% c("blearner","variable")) )
     stop("Parameter type has to be 'blearner' or 'variable'.")
  
  if( !(is.numeric(nbars) && nbars > 0 && length(nbars) == 1) )
    stop("Parameter nbars has to be a positive integer.")
  
  if( !(is.numeric(maxchar) && maxchar > 0 && length(maxchar) == 1) )
    stop("Parameter maxchar has to be a positive integer.")
  
  if( hasArg(xlab) || hasArg(ylab) )
    stop("xlab, ylab already defined by default.")  
  
  
  ### --------------------------------------------------
  ### set range for x-axis
  if( !("xlim" %in% names(args)) ) {
    if( percent ) xlim <- c(0,1)
    else xlim <- c(0,sum(x))
  }  
  
  ### set x-axis label and transform values to percentages of necessary
  if( percent ) {
    xlab <- "In-bag Risk Reduction (%)"
    x    <- x / sum(x)
  } else xlab <- "In-bag Risk Reduction"
  
  
  ### --------------------------------------------------
  ### create data.frame for all values shown in barchart
  plot_data <- data.frame(
    reduction = as.numeric(x), 
    blearner  = names(x),
    variable  = attr(x, "variable_names"),
    selprob   = attr(x, "selprob")
  )
  # sort rows by risk reduction (ascending)
  plot_data_sorted <- plot_data[order(plot_data$reduction),]
  
  
  ### --------------------------------------------------
  ### if number of baselearners/variables exceeds nbars, aggregate some values
  plot_data_rows <- nrow(plot_data_sorted)
  
  if( plot_data_rows > nbars ) {
    # add new level OTHER to factor variable blearner/variable
    plot_data_sorted[, type] <- factor(plot_data_sorted[, type],
      levels = c(levels(plot_data_sorted[, type]), "other"))  
    # set baselearner/variable names to OTHER depending on nbars
    plot_data_sorted[1:(plot_data_rows-nbars+1), type] <- "other"  
        
    # add up risk reduction and selprobs for OTHER
    plot_data_sorted[1:(plot_data_rows-nbars+1), "reduction"] = 
      sum( head(plot_data_sorted$reduction, plot_data_rows-nbars+1) )
    plot_data_sorted[1:(plot_data_rows-nbars+1), "selprob"] = 
      sum( head(plot_data_sorted$selprob, plot_data_rows-nbars+1) )
    
    # use only number of observations corresponding to nbars
    plot_data_sorted = tail(plot_data_sorted, nbars)
  }
  
  
  ### --------------------------------------------------
  ### create bar labels
  # (maybe) truncate labels for bars
  plot_data_sorted[, type] <- sapply(as.character(plot_data_sorted[, type]), 
    FUN = function(name) {
      paste(strtrim(name, maxchar), if( nchar(name) < maxchar ) "" else "..") 
    }
  )
  # add selprobs to bar labels
  plot_data_sorted[, type] = paste( plot_data_sorted[, type], "\n", 
    paste0("sel. prob: ~", as.character(round(plot_data_sorted[, "selprob"], 
    digits = 2))) )

  # use an ordered factor for correct order of bars (descending)
  plot_data_sorted[, type] = ordered(plot_data_sorted[, type],
    levels = unique(plot_data_sorted[, type]))

  
  ### --------------------------------------------------
  ### create final plot depending on type
  if( type == "variable" ) {
    barchart(variable ~ reduction, groups = blearner, data = plot_data_sorted,
      horizontal = TRUE, xlab = xlab, ylab = "Variables", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      stack = TRUE, auto.key = list(), ...)
  } else {
    barchart(blearner ~ reduction, data = plot_data_sorted,
      horizontal = TRUE, xlab = xlab, ylab = "Baselearner", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      ...)
  }
}  
