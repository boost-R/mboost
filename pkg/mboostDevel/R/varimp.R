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
  
  var_names <-  variable.names(object)   
  var_order <- order( sapply(var_names, function(i) sum(explained[var_names==i]) ) )
  var_names <- ordered(var_names, levels = unique(var_names[var_order]))
  
  # add variable names per baselearner
  attr(explained, "variable_names") <- var_names
  
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
    xlim <- c(sum(x*(x<0)),sum(x*(x>=0)))     # Special calculation to include the case that risk reduction might be negative
    if( percent ) xlim <- xlim / sum(abs(x))
  }  
  
  ### set x-axis label and transform values to percentages of necessary
  if( percent ) {
    xlab <- "In-bag Risk Reduction (%)"
    x    <- x / sum(abs(x))
    if( any(x<0) ) warning("At least one risk reduction value is negative. Percental reduction is thus calculated as 'reduction/sum(abs(reduction))'.")
  } else xlab <- "In-bag Risk Reduction"
  
  
  ### --------------------------------------------------
  ### create data.frame for all values shown in barchart
  plot_data <- data.frame(
    reduction = as.numeric(x), 
    blearner  = names(x),
    variable  = attr(x, "variable_names"),
    selprob   = attr(x, "selprob")
  )
  
  ### --------------------------------------------------
  ### equip variable- or blearner - names with order depending on type and add "other" as lowest level
  if( type == "blearner" ) 
  {
    plot_data$blearner <- factor(plot_data$blearner, levels = unique(plot_data$blearner[order(plot_data$reduction)]))
  }
  
  ### --------------------------------------------------
  ### if number of baselearners/variables exceeds nbars, aggregate some values
  
  if( nlevels(plot_data[, type]) > nbars ) {
    # logical vector indicating rows to be subsumed
    others <- !(plot_data[, type] %in% tail(levels(plot_data[,type]), nbars-1))
    # add new level OTHER to factor variable blearner/variable
    for(i in c("blearner","variable")) plot_data[, i] <- factor(plot_data[, i], levels = c("other", levels(plot_data[, i])))
    # set baselearner/variable names to OTHER depending on nbars
    plot_data[others, c("blearner","variable")] <- "other"  
    # add up risk reduction and selprobs for OTHER
    plot_data[others, "reduction"] <- 
      sum( plot_data[others, "reduction"] )
    plot_data[others, "selprob"] <- 
      sum( plot_data[others, "selprob"] )
    
    # use only number of observations corresponding to nbars
    plot_data <- rbind( plot_data[ c(which(!others), which(others)[1]), ] )
    for(i in c("blearner","variable")) plot_data[, i] <- factor(plot_data[, i])
  }
  
  
  ### --------------------------------------------------
  ### create bar labels
  # (maybe) truncate labels for bars
  trunc_levels <- sapply(levels(plot_data[, type]), 
    FUN = function(name) {
      paste(strtrim(name, maxchar), if( nchar(name) < maxchar ) "" else "..") 
    } 
  )
  if(length(unique(trunc_levels))==nlevels(plot_data[, type])) levels(plot_data[, type]) <- trunc_levels
  else {
    warning(paste("'maxchar' to small to distiguish", type, "labels. Instead they are labeled with numberated according to their (first) apprearences in the model and 'o' for 'other'."))
    levels(plot_data[,type]) <- if(type=="blearner") order(as.numeric(x)) else sapply(levels(plot_data$variable), function(i) if(i=="other") "o" else which(unique(attr(x,"variable"))==i)[1])
  }
  
  # convert rounded selprobs to character
  plot_data[,"selprob"] <- as.character(round(plot_data[, "selprob"], digits = 2))
                                               
  # for type = "variable" accumulate selprobs per variable in alphabetical order of particular baselearners
  
  if(type == "variable") plot_data[,"selprob"] <- sapply(plot_data[, "variable"], 
                                                               function(i) do.call(function(...) paste(..., sep = " + "),
                                                                                   as.list( plot_data[plot_data$variable == i,"selprob"]
                                                                                            [order(plot_data[plot_data$variable == i,"blearner"])])))
  # add selprobs to bar labels
  levels(plot_data[, type]) <- paste( levels(plot_data[ ,type]), "\n", 
    paste0("sel. prob: ~", unique( plot_data[order(plot_data[, type]), "selprob"])) )

  # use an ordered factor for correct order of bars (descending)
  plot_data[, type] <- ordered(plot_data[, type])

  
  ### --------------------------------------------------
  ### create final plot depending on type
  if( type == "variable" ) {
    barchart(variable ~ reduction, groups = blearner, data = plot_data,
      horizontal = TRUE, xlab = xlab, ylab = "Variables", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      stack = TRUE, auto.key = list(), ...)
  } else {
    barchart(blearner ~ reduction, data = plot_data,
      horizontal = TRUE, xlab = xlab, ylab = "Baselearner", xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      ...)
  }
}  
