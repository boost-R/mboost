\name{wpbc}
\alias{wpbc}
\docType{data}
\title{ Wisconsin Prognostic Breast Cancer Data}
\description{

  Each record represents follow-up data for one breast cancer
  case.  These are consecutive patients seen by Dr. Wolberg
  since 1984, and include only those cases exhibiting invasive
  breast cancer and no evidence of distant metastases at the
  time of diagnosis. 

}
\usage{data("wpbc")}
\format{
  A data frame with 198 observations on the following 34 variables.
  \describe{
    \item{\code{status}}{a factor with levels \code{N} (nonrecur) and 
                         \code{R} (recur)}
    \item{\code{time}}{recurrence time (for \code{status == "R"}) or
                       disease-free time (for \code{status == "N"}). }
    \item{\code{mean_radius}}{radius (mean of distances from center to points on the perimeter) (mean).}
    \item{\code{mean_texture}}{texture (standard deviation of gray-scale values) (mean).}
    \item{\code{mean_perimeter}}{perimeter (mean).}
    \item{\code{mean_area}}{area (mean).}
    \item{\code{mean_smoothness}}{smoothness (local variation in radius lengths) (mean).}
    \item{\code{mean_compactness}}{compactness (mean).}
    \item{\code{mean_concavity}}{concavity (severity of concave portions of the contour) (mean).}
    \item{\code{mean_concavepoints}}{concave points (number of concave portions of the contour) (mean).}
    \item{\code{mean_symmetry}}{symmetry (mean).}
    \item{\code{mean_fractaldim}}{fractal dimension (mean).}
    \item{\code{SE_radius}}{radius (mean of distances from center to points on the perimeter) (SE).}
    \item{\code{SE_texture}}{texture (standard deviation of gray-scale values) (SE).}
    \item{\code{SE_perimeter}}{perimeter (SE).}
    \item{\code{SE_area}}{area (SE).}
    \item{\code{SE_smoothness}}{smoothness (local variation in radius lengths) (SE).}
    \item{\code{SE_compactness}}{compactness (SE).}
    \item{\code{SE_concavity}}{concavity (severity of concave portions of the contour) (SE).}
    \item{\code{SE_concavepoints}}{concave points (number of concave portions of the contour) (SE).}
    \item{\code{SE_symmetry}}{symmetry (SE).}
    \item{\code{SE_fractaldim}}{fractal dimension (SE).}
    \item{\code{worst_radius}}{radius (mean of distances from center to points on the perimeter) (worst).}
    \item{\code{worst_texture}}{texture (standard deviation of gray-scale values) (worst).}
    \item{\code{worst_perimeter}}{perimeter (worst).}
    \item{\code{worst_area}}{area (worst).}
    \item{\code{worst_smoothness}}{smoothness (local variation in radius lengths) (worst).}
    \item{\code{worst_compactness}}{compactness  (worst).}
    \item{\code{worst_concavity}}{concavity (severity of concave portions of the contour) (worst).}
    \item{\code{worst_concavepoints}}{concave points (number of concave portions of the contour) (worst).}
    \item{\code{worst_symmetry}}{symmetry (worst).}
    \item{\code{worst_fractaldim}}{fractal dimension (worst).}
    \item{\code{tsize}}{diameter of the excised tumor in centimeters.}
    \item{\code{pnodes}}{number of positive axillary lymph nodes observed at time of surgery.}
  }
}
\details{

    The first 30 features are computed from a digitized image of a
    fine needle aspirate (FNA) of a breast mass.  They describe
    characteristics of the cell nuclei present in the image.

    There are two possible learning problems: predicting \code{status} or predicting
    the time to recur.

	1) Predicting field 2, outcome: R = recurrent, N = non-recurrent
	- Dataset should first be filtered to reflect a particular
	endpoint; e.g., recurrences before 24 months = positive,
	non-recurrence beyond 24 months = negative.
	- 86.3% accuracy estimated accuracy on 2-year recurrence using
	previous version of this data.  

	2) Predicting Time To Recur (field 3 in recurrent records)
	- Estimated mean error 13.9 months using Recurrence Surface
	Approximation.

   The data are originally available from the UCI machine learning repository, see
   \url{http://www.ics.uci.edu/~mlearn/databases/breast-cancer-wisconsin/}.

}
\source{

  W. Nick Street, Olvi L. Mangasarian and William H. Wolberg (1995). 
  An inductive learning approach to prognostic prediction. 
  In A. Prieditis and S. Russell, editors, \emph{Proceedings of the
  Twelfth International Conference on Machine Learning}, pages
  522--530, San Francisco, Morgan Kaufmann.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, accepted.
  \url{http://www.imstat.org/sts/future_papers.html}

}
\examples{

    data("wpbc", package = "mboost")

    ### fit logistic regression model with 100 boosting iterations
    coef(glmboost(status ~ ., data = wpbc[,colnames(wpbc) != "time"], 
                  family = Binomial()))

}
\keyword{datasets}
