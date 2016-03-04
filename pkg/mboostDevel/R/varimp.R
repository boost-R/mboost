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
  # risk after each boosting-step
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
  
  # add variable names per baselearner to varimp-object
  var_names <- unname( variable.names(object) )
  var_order <- order( sapply(var_names, function(i) {
    sum(explained[var_names == i]) 
  }) )
  
  # add names as ordered factor (for order in graphical output)
  attr(explained, "variable_names") <- ordered(var_names, levels = 
      unique(var_names[var_order]))
  
  return(explained)
}


as.data.frame.varimp <- function(x, optional = FALSE, ...) {
  data.frame(
    reduction = as.numeric(x), 
    # blearner as ordered factor (corresponding to variable(_names))
    blearner  = ordered(names(x), levels = unique(names(x)[order(x)])),
    variable  = attr(x, "variable_names"),
    selprob   = attr(x, "selprob")
  )
}


plot.varimp <- function(x, percent = TRUE, type = "variable", 
  nbars = 10L, maxchar = 20L, xlim, auto.key, ...) {
  
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
  ### set plotting arguments
  # set range for x-axis
  if( !("xlim" %in% names(args)) ) {
    # special calculation to include the case that risk reduction is negative
    xlim <- c(sum(x*(x < 0)), sum(x*(x >= 0)))     
    if( percent ) xlim <- xlim / sum(abs(x))
  }  
  
  # set x-axis label and transform values to percentages if necessary
  if( percent ) {
    xlab <- "In-bag Risk Reduction (%)"
    x    <- x / sum(abs(x))
    # special handling if risk reduction is negative
    if( any(x < 0) ) 
      warning(paste("At least one risk reduction value is negative. Percental", 
        "reduction is thus calculated as 'reduction/sum(abs(reduction))'."))
  } else xlab <- "In-bag Risk Reduction"
  
  
  ### --------------------------------------------------
  ### create data.frame for all values shown in barchart
  plot_data <- data.frame(x)
  
  
  ### --------------------------------------------------
  ### if number of baselearners/variables exceeds nbars, aggregate some values
  if( nlevels(plot_data[, type]) > nbars ) {
    # logical vector indicating rows to be subsumed
    others <- !(plot_data[, type] %in% tail(levels(plot_data[, type]), nbars-1))    
    
    for( i in c("blearner","variable") ) {
      # add new level OTHER to factor variable blearner/variable
      plot_data[, i] <- 
        ordered(plot_data[, i], levels = c("other", levels(plot_data[, i])))
      # set baselearner/variable names to OTHER depending on nbars
      plot_data[others, i] <- "other"
    }    
    
    # add up risk reduction and selprobs for OTHER
    plot_data[others, "reduction"] <- sum( plot_data[others, "reduction"] )
    plot_data[others, "selprob"]   <- sum( plot_data[others, "selprob"] )
    
    # use only number of observations corresponding to nbars
    plot_data <- plot_data[c(which(!others), which(others)[1]),] 
    
    # update ordered factors
    for( i in c("blearner","variable") ) 
      plot_data[, i] <- factor(plot_data[, i])
  }  
  
  
  ### --------------------------------------------------
  ### create bar labels
  # (maybe) truncate bar labels for baselearner or variable names
  trunc_levels <- sapply(levels(plot_data[, type]), 
    function(name) {
      paste(strtrim(name, maxchar), if( nchar(name) < maxchar ) "" else "..") 
    } 
  )  
  
  # additional warning message if bar labels are not distinct (after truncation)
  if( length(unique(trunc_levels)) != nlevels(plot_data[, type]) ) {
    warning(paste("'maxchar' set too small to distinguish", type, "labels.",
      "Different variables are probably summed up.",
      "Value for argument 'maxchar' should be incremented."))
  }
  levels(plot_data[, type]) <- trunc_levels
  
  # convert rounded selprobs to character
  plot_data[, "selprob"] <- 
    as.character(round(plot_data[, "selprob"], digits = 2))  
  
  # for type = "variable" additionally accumulate selprobs per variable 
  # in order of risk reduction of involved baselearners
  if( type == "variable" ) {
    # reverse order of baselearners (larger stacks first)
    plot_data[, "blearner"] <- 
      ordered(plot_data[, "blearner"], rev(levels(plot_data[, "blearner"])))
    # sum up selprobs
    plot_data[, "selprob"] <- sapply(plot_data[, "variable"], function(i) {
      do.call( function(...) paste(..., sep = " + "), as.list( 
        plot_data[plot_data$variable == i, "selprob"]
        [order(plot_data[plot_data$variable == i, "blearner"])] 
      ) )
    })
  }
  
  # add selprobs to bar labels
  selprob_labels = unlist( unique(plot_data[, c(type, "selprob")])[2] )
  
  levels(plot_data[, type]) = paste0( levels(plot_data[, type]), "\n ",
    "sel. prob: ~", selprob_labels[order(unique(plot_data[, type]))] )  
  

  ### --------------------------------------------------
  ### activate auto.key if bars are stacked (only other bars than 'others')
  if( !("auto.key" %in% names(args)) ) {
    auto.key <- nlevels(plot_data[, "variable"]) < nrow(plot_data)
  }
  
  
  ### --------------------------------------------------
  ### create final plot depending on type
  if( type == "variable" ) {
    barchart(variable ~ reduction, groups = blearner, data = plot_data,
      horizontal = TRUE, xlab = xlab, ylab = "Variables", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))), 
      stack = TRUE, auto.key = auto.key, ...)
  } else {
    barchart(blearner ~ reduction, data = plot_data,
      horizontal = TRUE, xlab = xlab, ylab = "Baselearner", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      ...)
  }
}  
