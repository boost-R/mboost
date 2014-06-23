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
\alias{GammaReg}
\alias{CoxPH}
\alias{QuantReg}
\alias{ExpectReg}
\alias{NBinomial}
\alias{PropOdds}
\alias{Weibull}
\alias{Loglog}
\alias{Lognormal}
\alias{AUC}
\alias{Gehan}
\alias{Hurdle}
\alias{Multinomial}
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
       nuisance = function() return(NA),
       name = "user-specified", fW = NULL,
       response = function(f) NA,
       rclass = function(f) NA)
AdaExp()
AUC()
Binomial(link = c("logit", "probit"), ...)
GaussClass()
GaussReg()
Gaussian()
Huber(d = NULL)
Laplace()
Poisson()
GammaReg(nuirange = c(0, 100))
CoxPH()
QuantReg(tau = 0.5, qoffset = 0.5)
ExpectReg(tau = 0.5)
NBinomial(nuirange = c(0, 100))
PropOdds(nuirange = c(-0.5, -1), offrange = c(-5, 5))
Weibull(nuirange = c(0, 100))
Loglog(nuirange = c(0, 100))
Lognormal(nuirange = c(0, 100))
Gehan()
Hurdle(nuirange = c(0, 100))
Multinomial()
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
  \item{rclass}{ function to derive class predictions from conditional class
                   probabilities (for models with factor response variable).}
  \item{check_y}{ a function for checking and transforming
                  the class / mode of a response variable.}
  \item{nuisance}{ a function for extracting nuisance parameters from the family.}
  \item{weights}{ a character indicating what type of weights are
    allowed. These can be either arbitrary (non-negative) weights
    code{"any"}, only zero and one weights \code{"zeroone"},
    (non-negative) interger weights \code{"case"}, or no weights are
    allowed \code{"none"}.}
  \item{name}{ a character giving the name of the loss function for pretty printing. }
  \item{link}{ link function for binomial family. Alternatively,
               one may supply the name of a distribution (for example
               \code{link = "norm"}), parameters of which may
               be specified via the \dots argument.}
  \item{d}{ delta parameter for Huber loss function. If omitted, it is chosen adaptively.}
  \item{tau}{ the quantile or expectile to be estimated, a number strictly between 0 and 1.}
  \item{qoffset}{ quantile of response distribution to be used as offset.}
  \item{nuirange}{ a vector containing the end-points of the interval to be
          searched for the minimum risk w.r.t. the nuisance parameter.
	  In case of \code{PropOdds}, the starting values for
	  the nuisance parameters. }
  \item{offrange}{ interval to search for offset in.}
  \item{\dots}{ additional arguments to link functions.}
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
  returns the population minimizers evaluated at the response, i.e.,
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
to fitting a logistic model. Alternative link functions
can be specified via the name of the corresponding distribution, for
example \code{link = "cauchy"} lead to \code{\link{pcauchy}}
used as link function. This feature is still experimental and
not well tested.

However, the coefficients resulting from boosting with family
\code{Binomial(link = "logit")} are \eqn{1/2} of the coefficients of a logit model
obtained via \code{\link{glm}}. This is due to the internal recoding
of the response to \eqn{-1} and \eqn{+1} (see below).
  However, Buehlmann and Hothorn (2007) argue that the
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
  \eqn{\log(\mu) = \eta} is assumed. The default step-site \code{nu = 0.1}
  is probably too large for this family (leading to
  infinite residuals) and smaller values are more appropriate.

  \code{GammaReg()} implements a family for fitting nonnegative response
  variables. The implemented loss function is the negative Gamma
  log-likelihood with logarithmic link function (instead of the natural
  link).

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

  \code{Gehan()} implements rank-based estimation of survival data in an
  accelerated failure time model. The loss function is defined as the sum
  of the pairwise absolute differences of residuals. The response needs to
  be defined as \code{Surv(y, delta)}, where \code{y} is the observed survial
  time (subject to censoring) and \code{delta} is the non-censoring indicator
  (see \code{\link{Surv}} for details). For details on \code{Gehan()} see
  Johnson and Long (2011).

  Hurdle models for zero-inflated count data can be fitted by using a combination
  of the \code{Binomial()} and \code{Hurdle()} families. While the \code{Binomial()}
  family allows for fitting the zero-generating process of the Hurdle model,
  \code{Hurdle()} fits a negative binomial regression model to the non-zero
  counts. Note that the specification of the Hurdle model allows for using
  \code{Binomial()} and \code{Hurdle()} independently of each other.

  Linear or additive multinomial logit models can be fitted using
  \code{Multinomial()}; although is family requires some extra effort for
  model specification (see example).  More specifically, the predictor must
  be in the form of a linear array model (see \code{\link{\%O\%}}).  Note
  that this family does not work with tree-based base-learners at the
  moment.  The class corresponding to the last level of the factor coding
  the response is uses as reference class.

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

    Nora Fenske, Thomas Kneib, and Torsten Hothorn (2011),
    Identifying risk factors for severe childhood malnutrition by
    boosting additive quantile regression.
	\emph{Journal of the American Statistical Association}, \bold{106}:494-510.

    Yoav Freund and Robert E. Schapire (1996),
    Experiments with a new boosting algorithm.
    In \emph{Machine Learning: Proc. Thirteenth International Conference},
    148--156.

    Shuangge Ma and Jian Huang (2005), Regularized ROC method for
    disease classification and biomarker selection with microarray
    data. \emph{Bioinformatics}, \bold{21}(24), 4356--4362.

    Whitney K. Newey and James L. Powell (1987),
    Asymmetric least squares estimation and testing.
    \emph{Econometrika}, \bold{55}, 819--847.

    Matthias Schmid and Torsten Hothorn (2008),
    Flexible boosting of accelerated failure time models.
    \emph{BMC Bioinformatics}, \bold{9}(269).

    Matthias Schmid, Sergej Potapov, Annette Pfahlberg,
    and Torsten Hothorn (2010). Estimation and regularization techniques for
    regression models with multidimensional prediction functions.
    \emph{Statistics and Computing}, \bold{20}, 139-150.

    Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
    (2014). Model-based Boosting in R: A Hands-on Tutorial Using the R
    Package mboost. \emph{Computational Statistics}, \bold{29}, 3--35.\cr
    \url{http://dx.doi.org/10.1007/s00180-012-0382-5}

    Available as vignette via: \code{vignette(package = "mboostDevel", "mboost_tutorial")}

    Brent A. Johnson and Qi Long (2011) Survival ensembles by the sum of pairwise
    differences with application to lung cancer microarray studies.
    \emph{Annals of Applied Statistics}, \bold{5}, 1081--1101.

  }
\seealso{\code{\link{mboost}} for the usage of \code{Family}s. See
    \code{\link{boost_family-class}} for objects resulting from a call to \code{Family}. }
\author{
    \code{ExpectReg()} was donated by Fabian Sobotka.
    \code{AUC()} was donated by Fabian Scheipl.
}
\examples{

    Laplace()

    MyGaussian <- function(){
           Family(ngradient = function(y, f, w = 1) y - f,
           loss = function(y, f) (y - f)^2,
           name = "My Gauss Variant")
    }

    ### fitting multinomial logit model via a linear array model
    X0 <- K0 <- diag(nlevels(iris$Species) - 1)
    colnames(X0) <- levels(iris$Species)[-nlevels(iris$Species)]
    mlm <- mboost(Species ~ bols(Sepal.Length, df = 2) \%O\%
                            buser(X0, K0, df = 2), data = iris,
                  family = Multinomial())
    head(round(predict(mlm, type = "response"), 2))

    \donttest{
    ### compare results with nnet::multinom
    mlmn <- multinom(Species ~ Sepal.Length, data = iris)
    max(abs(fitted(mlm[1000], type = "response") -
            fitted(mlmn, type = "prob")))
    }

}
\keyword{models}
