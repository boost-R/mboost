\name{bodyfat}
\alias{bodyfat}
\docType{data}
\title{ Prediction of Body Fat by Skinfold Thickness, Circumferences, and
        Bone Breadths}
\description{
    For 71 healthy female subjects, body fat measurements and several
    anthropometric measurements are available for predictive modelling
    of body fat.
}
\usage{data("bodyfat")}
\format{
  A data frame with 71 observations on the following 10 variables.
  \describe{
    \item{\code{age}}{age in years.}
    \item{\code{DEXfat}}{body fat measured by DXA, response variable.}
    \item{\code{waistcirc}}{waist circumference.}
    \item{\code{hipcirc}}{hip circumference.}
    \item{\code{elbowbreadth}}{breadth of the elbow.}
    \item{\code{kneebreadth}}{breadth of the knee.}
    \item{\code{anthro3a}}{sum of logarithm of three anthropometric measurements.}
    \item{\code{anthro3b}}{sum of logarithm of three anthropometric measurements.}
    \item{\code{anthro3c}}{sum of logarithm of three anthropometric measurements.}
    \item{\code{anthro4}}{sum of logarithm of three anthropometric measurements.}
  }
}
\details{

  Garcia et al. (2005) report on the development of predictive regression equations
  for body fat content by means of common anthropometric
  measurements which were obtained for 71 healthy German women.
  In addition, the women's body composition was measured by
  Dual Energy X-Ray Absorptiometry (DXA). This reference method
  is very accurate in measuring body fat but finds little applicability
  in practical environments, mainly because of high costs and the
  methodological efforts needed. Therefore, a simple regression equation
  for predicting DXA measurements of body fat is of special interest for the practitioner.
  Backward-elimination was applied to select
  important variables from the available anthropometrical measurements, and
  Garcia (2005) report a final linear model utilizing
  hip circumference, knee breadth and a compound covariate which is defined as
  the sum of log chin skinfold, log triceps skinfold and log subscapular skinfold.

}
\source{

    Ada L. Garcia, Karen Wagner, Torsten Hothorn, Corinna Koebnick,
    Hans-Joachim F. Zunft and Ulrike Trippo (2005),
    Improved prediction of body fat by measuring skinfold
    thickness, circumferences, and bone breadths. \emph{Obesity Research},
    \bold{13}(3), 626--634.

    Peter Buehlmann and Torsten Hothorn (2007),
    Boosting algorithms: regularization, prediction and model fitting.
    \emph{Statistical Science}, \bold{22}(4), 477--505.

  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid
  (2012). Model-based Boosting in R: A Hands-on Tutorial Using the R
  Package mboost. \emph{Department of Statistics, Technical Report No. 120}.\cr
  \url{http://epub.ub.uni-muenchen.de/12754/}

  Available as vignette via: vignette(package = "mboost", "mboost_tutorial")
}
\examples{

    data("bodyfat", package = "mboost")

    ### final model proposed by Garcia et al. (2005)
    fmod <- lm(DEXfat ~ hipcirc + anthro3a + kneebreadth, data = bodyfat)
    coef(fmod)

    ### plot additive model for same variables
    amod <- gamboost(DEXfat ~ hipcirc + anthro3a + kneebreadth,
                     data = bodyfat, baselearner = "bbs")
    layout(matrix(1:3, ncol = 3))
    plot(amod[mstop(AIC(amod, "corrected"))], ask = FALSE)

}
\keyword{datasets}
