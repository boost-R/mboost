\name{varimp}
\alias{varimp}
\alias{varimp.mboost}
\alias{plot.varimp}
\alias{as.data.frame.varimp}

\title{ Variable Importance }
\description{
In-bag risk reduction per base-learner as variable importance for boosting.
}
\usage{
\method{varimp}{mboost}(object, ...)

\method{plot}{varimp}(x, percent = TRUE, type = c("variable", "blearner"), 
  blorder = c("importance", "alphabetical", "rev_alphabetical", "formula"),
  nbars = 10L, maxchar = 20L, xlab = NULL, ylab = NULL, xlim, auto.key, ...)
\method{as.data.frame}{varimp}(x, row.names = NULL, optional = FALSE, ...) 
}


\arguments{
  
  \item{object}{ an object of class \code{mboost}.}
  
  \item{x}{ an object of class \code{varimp}.}
  
  \item{percent}{ logical, indicating whether variable importance should be 
  specified in percent.}
  
  \item{type}{ a character string specifying whether to draw bars for variables
  (\code{"variable"}, default) or base-learners (\code{"blearner"}) in the 
  model (no effect for a glmboost object).}
  
  \item{blorder}{ a character string specifying the order of the base-learners
  in the plot. The default \code{"importance"} corresponds to the order of the
  base-learner importance, \code{"alphabetical"} and \code{"rev_alphabetical"} 
  to alphabetical order, respectively its reverse, and \code{"formula"} to 
  their order in the model formula.}
  
  \item{nbars}{ integer, maximum number of bars to be plotted. If \code{nbars} 
  is exceeded, least important variables / base-learners are summarized as
  "other".}
  
  \item{maxchar}{ integer, maximum number of characters in bar labels.}
  
  \item{xlab}{text for the x-axis label. If not set (default is \code{NULL}) 
    x-axis label is generated automatically depending on argument 
    \code{percent}.}  
  
  \item{ylab}{text for the y-axis label. If not set (default is \code{NULL}) 
    y-axis label is generated automatically depending on argument \code{type}.}
  
  \item{xlim}{ the x limits of the plot. Defaults are from \code{0} to total 
  reduction, or from \code{0} to \code{1} for \code{percent = TRUE}. (In case 
  of negative risk reductions, default limits are from total negative to total
  positve reduction, or the latter normalized by the total absolute reduction 
  for \code{percent = TRUE}.)}
  
  \item{auto.key}{ logical, or a list passed to \code{lattice::barchart}. By 
  default \code{auto.key=TRUE} provides automatically generated legends showing 
  the underlying base-learners in the stacked barchart 
  (\code{type = "variable"}). If there is an unique base-learner for each 
  variable(-interaction), \code{auto.key = FALSE} is default setting. 
  For \code{type = "blearner"} the argument has no effect at all.} 
  
  \item{...}{ additional arguments passed to \code{lattice::barchart}.}
  
  \item{row.names}{ NULL or a character vector giving the row names for the 
  data frame. Missing values are not allowed.}

  \item{optional}{logical. If TRUE, setting row names and converting column 
  names (to syntactic names: see make.names) is optional.}
  
}
\details{
This function extracts the in-bag risk reductions per boosting step of a 
fitted \code{mboost} model and accumulates it individually for each base-learner
contained in the model. This quantifies the individual contribution to risk 
reduction of each base-learner and can thus be used to compare the importance 
of different base-learners or variables in the model. Starting from offset only, 
in each boosting step risk reduction is computed as the difference between 
in-bag risk of the current and the previous model and is accounted for the 
base-learner selected in the particular step.

The results can be plotted in a bar plot either for the base-learners, or the
variables contained in the model. The bars are ordered according to variable 
importance. If their number  exceeds \code{nbars} the least important are 
summarized as "other". If bars are plotted per  variable, all base-learners 
containing the same variable will be accumulated in a stacked bar. This is of 
use for models including for example seperate base-learners for the linear and 
non-linear part of a covariate effect (see \code{?bbs} option 
\code{center=TRUE}). However, variable interactions are treated as individual 
variables, as their desired handling might depend on context.

As a comparison the selection frequencies are added to the respective 
base-learner labels in the plot (rounded to three digits). For stacked bars 
they are ordered accordingly.
}

\value{ An object of class \code{varimp} with available \code{plot} and
\code{as.data.frame} methods. 

Converting a \code{varimp} object results in a \code{data.frame} containing the 
risk reductions, selection frequencies and the corresponding base-learner and 
variable names as ordered \code{factors} (ordered according to their particular 
importance).}

\author{
Tobias Kuehn (\email{tobi.kuehn@gmx.de}), 
Almond Stoecker (\email{almond.stoecker@gmail.com})
}

\examples{

data(iris)
### glmboost with multiple variables and intercept
iris$setosa <- factor(iris$Species == "setosa")
iris_glm <- glmboost(setosa ~ 1 + Sepal.Width + Sepal.Length + Petal.Width +
                         Petal.Length,
                     data = iris, control = boost_control(mstop = 50), 
                     family = Binomial(link = c("logit")))
varimp(iris_glm)
### importance plot with four bars only
plot(varimp(iris_glm), nbars = 4)

### gamboost with multiple variables
iris_gam <- gamboost(Sepal.Width ~ 
                         bols(Sepal.Length, by = setosa) +
                         bbs(Sepal.Length, by = setosa, center = TRUE) + 
                         bols(Petal.Width) +
                         bbs(Petal.Width, center = TRUE) + 
                         bols(Petal.Length) +
                         bbs(Petal.Length, center = TRUE),
                     data = iris)
varimp(iris_gam)
### stacked importance plot with base-learners in rev. alphabetical order
plot(varimp(iris_gam), blorder = "rev_alphabetical")

### similar ggplot
\dontrun{ 
library(ggplot2)
ggplot(data.frame(varimp(iris_gam)), aes(variable, reduction, fill = blearner)) + 
    geom_bar(stat = "identity") + coord_flip() }
}
