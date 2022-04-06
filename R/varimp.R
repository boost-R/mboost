varimp <- function(object, ...)
  UseMethod("varimp")


varimp.mboost <- function(object, ...) {  

  ### which baselearners/variables were selected in boosting steps
  learner_names <- if( !(inherits(object, "glmboost")) ) 
    names(object$baselearner) else variable.names(object)  
  learner_selected <- object$xselect()
  
  
  ### --------------------------------------------------
  ### compute inbag risk reduction for each step
  riskdiff <- -1 * diff(risk(object))
  
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
  
  # add selection frequencies to varimp-object
  attr(explained, "selfreqs") <- sapply(seq_along(learner_names), function(i) {
    mean(learner_selected == i)
  })
  
  # add variable names per baselearner to varimp-object
  var_names <- unname( variable.names(object) )
  # for identifiability sort variable names in interactions 
  # (not required for glmboost interactions)
  var_names <- sapply(strsplit(var_names, ", "), function(x) {
    do.call( function(...) paste(... , sep = ", " ), 
             as.list(x[order(x)]) ) })
  var_order <- order( sapply(var_names, function(i) {
    sum(explained[var_names == i]) 
  }) )
  
  # add names as ordered factor (for order in graphical output)
  attr(explained, "variable_names") <- ordered(var_names, levels = 
      unique(var_names[var_order]))
  
  return(explained)
}


as.data.frame.varimp <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    reduction = as.numeric(x), 
    # blearner as ordered factor (corresponding to variable(_names))
    blearner  = ordered(names(x), levels = unique(names(x)[order(x)])),
    variable  = attr(x, "variable_names"),
    selfreq   = attr(x, "selfreqs")
  )
}


plot.varimp <- function(x, percent = TRUE, type = c("variable", "blearner"), 
  blorder = c("importance", "alphabetical", "rev_alphabetical", "formula"), 
  nbars = 10L, maxchar = 20L, xlab = NULL, ylab = NULL, xlim, auto.key, ...) {
  
  ### --------------------------------------------------
  ### get arguments
  args <- as.list(match.call())  
    
  type <- match.arg(type) 
  blorder <- match.arg(blorder)   
  
  
  ### --------------------------------------------------
  ### check arguments
  if( !inherits(x, "varimp") ) 
    stop(paste(deparse(substitute(x)), "is not of class varimp."))
  
  if( !(is.logical(percent)) )
     stop("Parameter percent has of type logical.")
  
  if( !(type %in% c("blearner","variable")) )
     stop("Parameter type has to be 'blearner' or 'variable'.")
  
  if( !(is.numeric(nbars) && nbars > 0 && length(nbars) == 1) )
    stop("Parameter nbars has to be a positive integer.")
  
  if( !(is.numeric(maxchar) && maxchar > 0 && length(maxchar) == 1) )
    stop("Parameter maxchar has to be a positive integer.")
  
  if( !(blorder %in% c("importance", "alphabetical", "rev_alphabetical",
                       "formula")) )
    stop(paste("Parameter blorder has to be one of 'importance',",
      "'alphabetical', 'rev_alphabetical' or 'formula'."))
  
  
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
    x <- x / sum(abs(x))
    # special handling if risk reduction is negative
    if( any(x < 0) ) 
      warning(paste("At least one risk reduction value is negative. Percental", 
        "reduction is thus calculated as 'reduction/sum(abs(reduction))'."))
  }
  
  # set axis labels
  xlab <- ifelse(is.null(xlab), ifelse(percent, "In-bag Risk Reduction (%)", 
    "In-bag Risk Reduction"), xlab) 
  ylab <- ifelse(is.null(ylab), ifelse(type == "variable", "Variables", 
    "Baselearner"), ylab) 
  
  
  
  ### --------------------------------------------------
  ### create data.frame for all values shown in barchart
  plot_data <- data.frame(x)

  
  ### --------------------------------------------------
  # specify baselearner order (if blorder != "importance", the default order set
  # via as.data.frame.varimp)
  plot_data[, "blearner"] <- switch(blorder,
    # importance is default, in this case no change
    importance       = plot_data[, "blearner"],
    # as blearner order is reverted again later, they are specified in 
    # reverted order here
    alphabetical     = ordered(names(x), rev(levels(ordered(names(x))))),
    rev_alphabetical = ordered(names(x)),
    formula          = ordered(names(x), rev(names(x)))
  )
  
  
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
    
    # add up risk reduction and selfreqs for OTHER
    plot_data[others, "reduction"] <- sum( plot_data[others, "reduction"] )
    plot_data[others, "selfreq"]   <- sum( plot_data[others, "selfreq"] )
    
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
  
  # convert rounded selfreqs to character
  plot_data[, "selfreq"] <- 
    as.character(round(plot_data[, "selfreq"], digits = 2))  
  
  # for type = "variable" additionally accumulate selfreqs per variable 
  # in order of risk reduction of involved baselearners
  if( type == "variable" ) {    
    # reverse order of baselearners (larger stacks first)
    # (also compare line 128ff - setting of baselearner order)
    plot_data[, "blearner"] <- 
      ordered(plot_data[, "blearner"], rev(levels(plot_data[, "blearner"])))
    # sum up selfreqs
    plot_data[, "selfreq"] <- sapply(plot_data[, "variable"], function(i) {
      do.call( function(...) paste(..., sep = " + "), as.list( 
        plot_data[plot_data$variable == i, "selfreq"]
        [order(plot_data[plot_data$variable == i, "blearner"])] 
      ) )
    })
  }
  
  # add selfreqs to bar labels
  selfreq_labels = unlist( unique(plot_data[, c(type, "selfreq")])[2] )
  
  levels(plot_data[, type]) = paste0( levels(plot_data[, type]), "\n ",
    "sel. freq: ~", selfreq_labels[order(unique(plot_data[, type]))] )  
  

  ### --------------------------------------------------
  ### activate auto.key if bars are stacked (only other bars than 'others')
  if( !("auto.key" %in% names(args)) ) {
    auto.key <- nlevels(plot_data[, "variable"]) < nrow(plot_data)
  }
  
  
  ### --------------------------------------------------
  ### create final plot depending on type
  if( type == "variable" ) {
    barchart(variable ~ reduction, groups = plot_data[, "blearner"], 
      data = plot_data, horizontal = TRUE, xlab = xlab, ylab = ylab, xlim =xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))), 
      stack = TRUE, auto.key = auto.key, ...)
  } else {
    barchart(blearner ~ reduction, data = plot_data,
      horizontal = TRUE, xlab = xlab, ylab = ylab, xlim = xlim,
      scales = list(x = list(tck = c(1,0), at = seq(0,sum(x), length.out = 5))),
      ...)
  }
}  

print.varimp <- function(x, ...) {
    cat("Variable importance (fraction of in-bag risk reduction):\n\n")
    attr(x, "selfreqs") <- NULL
    attr(x, "variable_names") <- NULL
    class(x) <- "numeric"
    NextMethod()
}
