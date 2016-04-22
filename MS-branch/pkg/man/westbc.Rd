\name{westbc}
\alias{westbc}
\docType{data}
\title{ Breast Cancer Gene Expression }
\description{
    Gene expressions for 7129 genes in 49 breast cancer samples and the 
    status of lymph node involvement.
}
\usage{data("westbc")}
\format{
  An object of class \code{exprSet} (see package \code{Biobase}).
}
\details{

  A full description of the data 
  can be found in West et al. (2001) and an application
  of boosted linear models is given by Buhlmann (2006).

}
\references{

  Peter Buhlmann (2006), Boosting for high-dimensional linear models.
  \emph{The Annals of Statistics}, \bold{34}(2), 559--583.

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, accepted.
  \url{http://www.imstat.org/sts/future_papers.html}

}
\source{
  Mike West, Carrie Blanchette, Holly Dressman, Erich Huang,
  Seiichi Ishida, Rainer Spang, Harry Zuzan, John A. Olson Jr.,
  Jeffrey R. Marks and Joseph R. Nevins (2001), 
  Predicting the clinical status of human breast cancer by using
  gene expression profiles, \emph{Proceedings of the National Academy of Sciences},
  \bold{98}, 11462-11467.
  \url{http://data.cgt.duke.edu/west.php}
}
\keyword{datasets}
