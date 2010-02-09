\name{Family}
\alias{Family}
\alias{AdaExp}
\alias{Binomial}
\alias{GaussClass}
\alias{GaussReg}
\alias{Gaussian}
\alias{Huber}
\alias{Laplace}
\alias{Poisson}
\alias{CoxPH}
\alias{QuantReg}
\alias{ExpectReg}
\alias{NBinomial}
\alias{PropOdds}
\alias{Weibull}
\alias{Loglog}
\alias{Lognormal}
\alias{AUC}
\title{ Gradient Boosting Families }
\description{
    \code{boost_family} objects provide a convenient way to specify loss functions
    and corresponding risk functions to be optimized by one of the boosting
    algorithms implemented in this package.
}
\usage{
Family(ngradient, loss = NULL, risk = NULL,
       offset = function(y, w)
           optimize(risk, interval = range(y), 
                    y = y, w = w)$minimum,
       check_y = function(y) y,
       weights = c("any", "none", "zeroone", "case"),
       nuisance = function() return(NULL),
       name = "user-specified", fW = NULL, 
       response = function(f) NA)
AdaExp()
AUC()
Binomial()
GaussClass()
GaussReg()
Gaussian()
Huber(d = NULL)
Laplace()
Poisson()
CoxPH()
QuantReg(tau = 0.5, qoffset = 0.5)
ExpectReg(tau = 0.5)
NBinomial(nuirange = c(0, 100))
PropOdds(nuirange = c(-0.5, -1), offrange = c(-5, 5))
Weibull(nuirange = c(0, 100))
Loglog(nuirange = c(0, 100))
Lognormal(nuirange = c(0, 100))
}
\arguments{
  \item{ngradient}{ a function with arguments \code{y}, \code{f} and \code{w} implementing the
                    \emph{negative} gradient of the \code{loss} function (which is to be minimized). }
  \item{loss}{ an optional loss function with arguments \code{y} and \code{f}. }
  \item{risk}{ an optional risk function with arguments \code{y}, \code{f} and \code{w} to be minimized (!),
               the weighted mean of the loss function by default. }
  \item{offset}{ a function with argument \code{y} and \code{w} (weights)
                 for computing a \emph{scalar} offset. }
  \item{fW}{ transformation of the fit for the diagonal weights matrix for an
             approximation of the boosting hat matrix for loss functions other than
             squared error.}
  \item{response}{ inverse link function of a GLM or any other transformation
                   on the scale of the response.}
  \item{check_y}{ a function for checking and transforming 
                  the class / mode of a response variable.}
  \item{nuisance}{ a function for extracting nuisance parameters from the family.}
  \item{weights}{ a character indicating if weights are allowed. }
  \item{name}{ a character giving the name of the loss function for pretty printing. }
  \item{d}{ delta parameter for Huber loss function. If omitted, it is chosen adaptively.}
  \item{tau}{ the quantile or expectile to be estimated, a number strictly between 0 and 1.}
  \item{qoffset}{ quantile of response distribution to be used as offset.}
  \item{nuirange}{ range of starting values for the nuisance parameter.}
  \item{offrange}{ interval to search for offset in.}
}
\details{
  The boosting algorithm implemented in \code{\link{mboost}} minimizes the 
  (weighted) empirical risk function \code{risk(y, f, w)} with respect to \code{f}. 
  By default, the risk function is the weighted sum of the loss function \code{loss(y, f)} 
  but can be chosen arbitrarily. The \code{ngradient(y, f)} function is the negative 
  gradient of \code{loss(y, f)} with respect to \code{f}.

  Pre-fabricated functions for the most commonly used loss functions are
  available as well. Buehlmann and Hothorn (2007) give a detailed
  overview of the available loss functions. The \code{offset} function
  returns the population minimisers evaluated at the response, i.e.,
  \eqn{1/2 \log(p / (1 - p))} for \code{Binomial()} or \code{AdaExp()}
  and \eqn{(\sum w_i)^{-1} \sum w_i y_i} for \code{Gaussian()} and the
  median for \code{Huber()} and \code{Laplace()}. A short summary of the
  available families is given in the following paragraphs:
  
\code{AdaExp()}, \code{Binomial()} and \code{AUC()} implement
families for binary classification. \code{AdaExp()} uses the
exponential loss, which essentially leads to the AdaBoost algorithm
of Freund and Schapire (1996). \code{Binomial()} implements the
negative binomial log-likelihood of a logistic regression model 
as loss function. Thus, using \code{Binomial} family closely corresponds 
to fitting a logistic model.
However, the coefficients resulting from boosting with family
\code{Binomial} are \eqn{1/2} of the coefficients of a logit model
obtained via \code{\link{glm}}. This is due to the internal recoding
of the response to \eqn{-1} and \eqn{+1} (see below).
Buehlmann and Hothorn (2007) argue that the
family \code{Binomial} is the preferred choice for binary
classification. For binary classification problems the response
\code{y} has to be a \code{factor}. Internally \code{y} is re-coded
to \eqn{-1} and \eqn{+1} (Buehlmann and Hothorn 2007).
\code{AUC()} uses \eqn{1-AUC(y, f)} as the loss function. 
The area under the ROC curve (AUC) is defined as 
\eqn{AUC = (n_{-1} n_1)^{-1} \sum_{i: y_i = 1} \sum_{j: y_j = -1} I(f_i > f_j)}.
Since this is not differentiable in \code{f}, we approximate the jump function
\eqn{I((f_i - f_j) > 0)} by the distribution function of the triangular 
distribution on \eqn{[-1, 1]} with mean \eqn{0}, similar to the logistic 
distribution approximation used in Ma and Huang (2005). 

  \code{Gaussian()} is the default family in \code{\link{mboost}}. It
  implements \eqn{L_2}Boosting for continuous response. Note
  that families \code{GaussReg()} and \code{GaussClass()} (for regression
  and classification) are deprecated now.
  \code{Huber()} implements a robust version for boosting with
  continuous response, where the Huber-loss is used. \code{Laplace()}
  implements another strategy for continuous outcomes and uses the
  \eqn{L_1}-loss instead of the \eqn{L_2}-loss as used by
  \code{Gaussian()}.

  \code{Poisson()} implements a family for fitting count data with
  boosting methods. The implemented loss function is the negative
  Poisson log-likelihood. Note that the natural link function
  \eqn{\log(\mu) = \eta} is assumed.

  \code{CoxPH()} implements the negative partial log-likelihood for Cox
  models. Hence, survival models can be boosted using this family.

  \code{QuantReg()} implements boosting for quantile regression, which is
  introduced in Fenske et al. (2009). \code{ExpectReg} works in analogy, 
  only for expectiles, which were introduced to regression by Newey and Powell (1987).

  Families with an additional scale parameter can be used for fitting
  models as well: \code{PropOdds()} leads to proportional odds models
  for ordinal outcome variables. When using this family, an ordered set of
  threshold parameters is re-estimated in each boosting iteration.
  \code{NBinomial()} leads to regression models with a negative binomial
  conditional distribution of the response. \code{Weibull()}, \code{Loglog()},
  and \code{Lognormal()} implement the negative log-likelihood functions
  of accelerated failure time models with Weibull, log-logistic, and
  lognormal distributed outcomes, respectively. Hence, parametric survival
  models can be boosted using these families. For details see Schmid and
  Hothorn (2008) and Schmid et al. (2010).

}
\section{Warning}{
  The coefficients resulting from boosting with family
  \code{Binomial} are \eqn{1/2} of the coefficients of a logit model
  obtained via \code{\link{glm}}. This is due to the internal recoding
  of the response to \eqn{-1} and \eqn{+1} (see above). 
  
  For \code{AUC()}, variables should be centered and scaled and observations with weight > 0 must not contain missing values. 
  The estimated coefficients for \code{AUC()} have no probabilistic interpretation.      
}
\value{
  An object of class \code{boost_family}.
}
\references{
    Peter Buehlmann and Torsten Hothorn (2007),
    Boosting algorithms: regularization, prediction and model fitting.
    \emph{Statistical Science}, \bold{22}(4), 477--505.

    Nora Fenske, Thomas Kneib, and Torsten Hothorn (2009).
    Identifying risk factors for severe childhood malnutrition by
    boosting additive quantile regression. Technical Report Nr. 52,
    Institut fuer Statistik, LMU Muenchen.
    \url{http://epub.ub.uni-muenchen.de/10510/}

    Yoav Freund and Robert E. Schapire (1996),
    Experiments with a new boosting algorithm.
    In \emph{Machine Learning: Proc. Thirteenth International Conference},
    148--156.

	Shuangge Ma and Jian Huang (2005),
	Regularized ROC method for disease classification and biomarker
	selection with microarray data. 
	\emph{Bioinformatics}, \bold{21}(24), 4356--4362.

    Whitney K. Newey and James L. Powell (1987),
    Asymmetric least squares estimation and testing.
    \emph{Econometrika}, \bold{55}, 819--847.

    Matthias Schmid and Torsten Hothorn (2008), 
    Flexible boosting of accelerated failure time models. 
    \emph{BMC Bioinformatics}, \bold{9}(269).

    Matthias Schmid, Sergej Potapov, Annette Pfahlberg, 
    and Torsten Hothorn (2010). Estimation and regularization techniques for 
    regression models with multidimensional prediction functions. 
    \emph{Statistics and Computing}, in press.

  }
\seealso{\code{\link{mboost}} for the usage of \code{Family}s. See
    \code{\link{boost_family-class}} for objects resulting from a call to \code{Family}. }
\author{
    \code{ExpectReg()} was donated by Fabian Sobotka.
    \code{AUC()} was donated by Fabian Scheipl.
}
\examples{

    Laplace()

    Family(ngradient = function(y, f) y - f,
           loss = function(y, f) (y - f)^2,
           name = "My Gauss Variant")

}
\keyword{models}
