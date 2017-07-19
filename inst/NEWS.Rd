\name{NEWS}
\title{News for Package 'mboost'}

\section{Changes in mboost version 2.8-1 (2017-07-19)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Added all possible options to the specific boosting functions
        instead of passing the options via \code{...} to \code{mboost_fit}.
        Closes \href{https://github.com/boost-R/mboost/issues/81}{#81}.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Minor speed ups in \code{df2lambda} (i.e., when computing penalty parameter 
        for the defined degrees of freedom). Changes proposed by Benjamin 
        Christoffersen.
      \item Updated kernel boosting reference. 
        Closes \href{https://github.com/boost-R/mboost/issues/84}{#84}.  
      \item Rebuilt package with LF instead of CRLF to fix \file{cleanup} script
        as requested by CRAN. Fixes \href{https://github.com/boost-R/mboost/issues/81}{#82}
      \item Use "old" definition of degrees of freedom in 
        \code{vignette("mboost", package = "mboost")} to make results reproducible.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Fix handling of missing values in \code{mboost} and \code{gamboost}
        when weights are specified. Fixes \href{https://github.com/boost-R/mboost/issues/80}{#80}.
    }
  }
}
  
\section{Changes in mboost version 2.8-0 (2017-05-04)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Models with zero steps (i.e., models containing only the offset) can
        now be fitted. Furthermore, cross-validation can now also select a model
        without base-learners.
        Fixes \href{https://github.com/boost-R/mboost/issues/64}{#64}, 
        \href{https://github.com/boost-R/mboost/issues/66}{#66}, and 
        \href{https://github.com/boost-R/mboost/issues/69}{#69}.        
      \item \code{Binomial} now uses link functions by making use of \code{make.link}. 
        Furthermore, an alternative implementation of \code{Binomial} models 
        along the lines of the \code{glm} implementation can be used via 
        \code{Binomial(type = "glm")}. Additionally, it works not only with a 
        two-level factor but also with a two-column matrix containing the number 
        of successes and number of failures. 
        Fixes \href{https://github.com/boost-R/mboost/issues/34}{#34}, 
        \href{https://github.com/boost-R/mboost/issues/63}{#63} and
        \href{https://github.com/boost-R/mboost/issues/65}{#65}.
      \item Added new base-learner \code{bkernel} for kernel boosting as described in
      
        S. Friedrichs, J. Manitz, P. Burger, C.I. Amos, A. Risch, J.C. Chang-Claude, 
        H.E. Wichmann, T. Kneib, H. Bickeboeller, and B. Hofner (2017), 
        Pathway-Based Kernel Boosting for the Analysis of Genome-Wide Association Studies.
        \emph{Computational and Mathematical Methods in Medicine}. 2017(6742763), 1-17. 
        \doi{10.1155/2017/6742763}.
        
      \item Removed check if \code{df2lambda} is stable. Hence,  
        \code{options(mboost_check_df2lambda)} (introduced in mboost 2.5-0) 
        is no longer used. 
        Closes \href{https://github.com/boost-R/mboost/issues/26}{#26}.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Added Andreas Mayr as contributor.
      \item Updated references and added reference to \code{citation("mboost")}.
      \item Fixed code of India example, which can be used to reproduce the data analysis 
        presented in
        
        N. Fenske, T. Kneib, and T. Hothorn (2011), 
        Identifying risk factors for severe childhood malnutrition by boosting 
        additive quantile regression.
        \emph{Journal of the American Statistical Association}, \bold{106}:494-510.
	      
	    (see \code{system.file("India_quantiles.R", package = "mboost")})
	  \item Fixed package citation.
	  \item Register C routines to make CRAN happy (again).
        Fixes \href{https://github.com/boost-R/mboost/issues/77}{#77}. 
    }
  }
  \subsection{Bug-fixes}{
     \itemize{
       \item Make sure that \code{family = Multinomial} is only used with Kronecker 
         product base-learners. Fixes \href{https://github.com/boost-R/mboost/issues/46}{#46}.
      \item Use argument \code{PACKAGE} in \code{.Call}. 
         Fixes \href{https://github.com/boost-R/mboost/issues/72}{#72}.
      \item If \code{center} is specified as boolean value in \code{bols}, we now throw an error. 
         Fixes \href{https://github.com/boost-R/mboost/issues/70}{#70}.
      \item  Fixed \code{AUC} family which expected fit to be equal to a constant in the 
        first iteration.
      \item Check for new data, e.g., in \code{predict}, was broken.
        Fixes \href{https://github.com/boost-R/mboost/issues/68}{#68}. 
     \item Make sure that \code{newdata} is discarded in \code{fitted}. 
        Fixes \href{https://github.com/boost-R/mboost/issues/76}{#76}.
     }
  }
}

\section{Changes in mboost version 2.7-0 (2016-11-23)}{
  \subsection{User-visible changes}{
    \itemize{
      \item New \code{Cindex} family to optimize survival models w.r.t. the concordance index.
      Fixes \href{https://github.com/boost-R/mboost/issues/53}{#53}.
      \item Added function \code{varimp} to extract variable importance. 
      A dedicated plot function exists (\code{plot(varimp())}). Code was provided 
      by Tobias Kuehn and Almond Stoecker. 
      See pull request \href{https://github.com/boost-R/mboost/pull/29}{#29}.
      \item Improved \code{plot} function for boosting models:
      \itemize{
        \item \code{plot} fails earlier in case of multiple \code{levelplots},
          i.e., maps (thanks to Mikko Korpela).
          See pull request \href{https://github.com/boost-R/mboost/pull/39}{#39}.
        \item Provide sensible defaults for \code{xlab} and \code{ylab} and allow 
          user-specified axis labels for bi- and multivariate plots.
          Fixes \href{https://github.com/boost-R/mboost/issues/51}{#51}.
        \item Export \code{plot} functions (\code{plot.glmboost}, \code{plot.mboost}, 
          \code{lines.mboost}, \code{plot.varimp} and \code{plot.cvrisk}) for 
          better usability and visibility.  
      }
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Updated manual regarding the usage of families and clarified the usage 
        of argument \code{qoffset}.
      \item Updated manual for base-learners:
      \itemize{
        \item Highlight that \code{x} should be centered if \code{bols(x, intercept = FALSE)} 
          is used.
        \item Discourage using \code{bbs(, constraint != "none")}; Preferably use 
          \code{bmono} for constrained effect estimates.
          Fixes \href{https://github.com/boost-R/mboost/issues/36}{#36}.
      }
      \item Improved vignettes (thanks to Mikko Korpela).
        See pull request \href{https://github.com/boost-R/mboost/pull/38}{#38}.
    }
  }
  \subsection{Bug-fixes}{
     \itemize{
      \item Solve potential problem with \code{IPCweights()}. Fixes
      \href{https://github.com/boost-R/mboost/issues/54}{#54}.
       \item Drop unobserved factor levels from \code{bols()}. Fixes
      \href{https://github.com/boost-R/mboost/issues/47}{#47}.
      \item Adapt \code{btree} to changes introduced in package \pkg{party}.
        Fixes \href{https://github.com/boost-R/mboost/issues/58}{#58}.
      \item Improved \code{cvrisk} to be more robust in various use cases 
        (thanks to Mikko Korpela). 
        See pull request \href{https://github.com/boost-R/mboost/pull/42}{#42}.
      \item Be more careful regarding namespace scoping rules.
        Fixes \href{https://github.com/boost-R/mboost/issues/45}{#45}.
      }
   }
}      

\section{Changes in mboost version 2.6-0 (2016-03-11)}{
  \subsection{User-visible changes}{
    \itemize{
      \item New maintainer: Benjamin Hofner follows Torsten Hothorn as
      maintainer.
      \item Package \pkg{party} is now imported. \pkg{mboost} no longer
      directly relies on unexported functions.
      \item  Allow extrapolation for predictions if kronecker products,
      tensor  products or sums are used. Fixes
      \href{https://github.com/boost-R/mboost/issues/23}{#23}.
      \item Development now hosted entirely on github as
      \href{https://github.com/boost-R/mboost/}{boost-R/mboost}.
      \item Started using \pkg{testthat}.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Improved checks for \code{newdata}: Warnings are no longer
      issued if data has just different types of numeric values (i.e.,
      integer vs. double). Resolves issue
      \href{https://github.com/boost-R/mboost/issues/17}{#17}.
      \item Fixed \file{CITATION} by removing duplicated string 'R package
      version' (spotted by Heidi Seibold).
      \item \code{predict}: Improve warning when \code{length(offset) >
	1}. Closes issue
      \href{https://github.com/boost-R/mboost/issues/20}{#20}.
      \item Added test coverage using package \pkg{covr}.
      \item Better error handling in \code{cvrisk} also for parallel
      processes.
      \item Suppress warning of \code{rankMatrix}. (Resolves issue
      \href{https://github.com/boost-R/mboost/issues/24}{#24}).
      \item Stop exporting internal functions for \pkg{FDboost}. Use
      \code{mboost_intern()} instead. Caution: Do not use this function.
    }
  }
  \subsection{Bug-fixes}{
     \itemize{
       \item Handling of missing values has been improved. Resolves issue
       \href{https://github.com/boost-R/mboost/issues/12}{#12}.
       \item Minor bug fixed in vignette
       \file{mboost_illustrations.Rnw}.
       \item Throw an error if model cannot be fitted. Fixes issue
       \href{https://github.com/boost-R/mboost/issues/18}{#18}.
       \item Fixed bug in bkronecker with dense matrices. Resolves issue
       \href{https://github.com/boost-R/mboost/issues/30}{#30}.
     }
   }
}

\section{Changes in mboost version 2.5-0 (2015-08-13)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Added documentation for \code{plot.mboost} function and moved
      documentation of \code{plot.glmboost} to the same help page.
      Resolves issue \href{https://github.com/boost-R/mboost/issues/14}{#14}.
      \item \code{bbs} and \code{bmono} no longer allow data outside of
      the \code{boundary.knots} during model fitting.
      \item Predictions for \code{bbs} and \code{bmono} now use linear
      extrapolation (user request inspired by
      \code{mgcv::Predict.matrix.pspline.smooth}).
      \item Better handling of errors in (single) folds of \code{cvrisk}:
      results of folds without errors are used and a \code{warning} is
      issued.
      \item Parallel computing via \code{mclapply}: Set
      \code{mc.preschedule = FALSE} per default.
      \item Added new option \code{options(mboost_check_df2lambda =
	TRUE)}, which controls if a stability check in \code{df2lambda}
      is performed. If set to \code{FALSE} this might speed up the
      computation of \code{df2lambda} especially with large design
      matrices.
      \item Prediction now also possible with \code{newdata = list()}.
      Resolves issue
      \href{https://github.com/boost-R/mboost/issues/15}{#15}.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item \code{PropOdds()}: Updated manual for proportional odds model.
      \item \code{Multinomial()}: Updated manual for multinomial logit
      model. Predictions for new data are now
      possible (resolves issue
      \href{https://github.com/boost-R/mboost/issues/13}{#13}, thanks to
      Sarah Brockhaus).
      \item \file{inst/CITATION}: Added subheadings and
      tutorial paper.
      \item Stopped computing the singular vectors in \code{df2lambda}
      as the singular values are sufficient and as
      \dQuote{computing the singular vectors is the slow part for large
	matrices} (proposed by Fabian Scheipl).
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Fixed bug in \code{PropOdds()} which occurred if
      \code{offset} was specified: nuisance parameters \code{delta}
      and \code{sigma} were not properly initialized (spotted by Madlene Nussbaum).
      \item Bug in \code{plot.mboost()} fixed which occurred if a factor
      with equal effect estimates for different categories was plotted.
      \item Bug in \code{df2lambda} fixed: Make sure that \code{A} is
      symmetric if it is \code{Matrix}-object (spotted by Souhaib Ben
      Taieb).
      \item Bug in \code{df2lambda} fixed. Design matrices were always
      assumed to be of full rank.
      \item Truncate output of complete data structure when model is
      printed. Resolves issue
      \href{https://github.com/boost-R/mboost/issues/11}{#11}.
      \item  Adhere to CRAN policies regarding import of base packages
      (closes \href{https://github.com/boost-R/mboost/issues/9}{#9}).
    }
  }
}

\section{Changes in mboost version 2.4-2 (2015-02-12)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Export \code{df2lambda}, \code{hyper_bbs} and \code{bl_lin}
        to make package \pkg{FDboost} happy. Note: These functions
	usually should not be called directly by users.
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Added Hothorn et al (2010) to \file{inst/CITATION}
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Changes in \file{inst/CITATION} to make CRAN happy:
        Citations can now be extracted without the need to install the
        package.
      \item Removed \code{EISPACK = FALSE} from \code{eigen()} as the
        argument is defunct and ignored.
      \item Changed \code{require} to \code{requireNamespace}
    }
  }
}

\section{Changes in mboost version 2.4-1 (2014-12-15)}{
  \subsection{Miscellaneous}{
    \itemize{
      \item Moved generic definition of \code{selected} to \pkg{stabs}
      which is required anyway (thus, \pkg{stabs} >= 0.5-0 is now
      required)
      \item load AML dataset (\file{AML_Bullinger.rda}) from package
      \pkg{TH.data}
      \item Updated references (for stability selection, confidence
      intervals and constrained regression)
      \item fixed \file{inst/CITATION}
      \item Refer to \code{news(package = "mboost")} instead of to the
      \file{NEWS} file.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item Cross-validation was potentially wrong for \code{CoxPH()}
      models. Users can now choose if they want the naive
      cross-validation or the improved version by Verweij and van
      Houwelingen (1993); (spotted by Holger Reulen <hreulen _ at _
      uni-goettingen.de>)
      \item Examples in \code{\\dontrun} are now executable and all
      dependencies are properly stated in \file{DESCRIPTION}
    }
  }
}

\section{Changes in mboost version 2.4-0 (2014-10-02)}{
  \subsection{User-visible changes}{
    \itemize{
      \item Added \code{confint} function to compute (bootstrap)
      confidence intervals together with plot and print methods

      \item \code{stabsel()} now depends on the new package \pkg{stabs}
      where the back end and methods such as \code{plot} and
      \code{print} are implemented

      \item Improved \code{plot} method for varying coefficients
      (\code{ylim} now suitable) and base-learners of factor variables.

      \item Tweaked \code{update} function: we now can turn the
      \code{trace} on and off, and specify the type of \code{risk} as
      well as the \code{oobweight} to \code{update()}

    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item Updated vignette \file{mboost_tutorial} to reflect latest
      changes in \pkg{mboost}.
      \item Changed plain text \file{NEWS} to \file{inst/NEWS.Rd}
      \item Removed links to archived package \pkg{mfp}.
      \item Explicitly specify the packages for functions that are
      implemented in packages that are listed as \code{Suggests:}, e.g
      we now use \code{party::ctree_control} etc.
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item \code{glmboost()$model.frame()} was broken
      \item \code{glmboost()$update()} was broken
      \item \code{predict()} for models with non-scalar offsets was
      broken
    }
  }
}

\section{Changes in mboost version 2.3-0 (2014-06-26)}{
  \subsection{User-visible changes}{
    \itemize{
      \item \code{stabsel} was recoded and now uses different
      terminology, much more options and a better tested code base

      \item new replacement function \code{mstop<-} as an alternative to
      \code{<mboost>[i]} (suggested by Achim Zeileis).

      \item \code{bmono}
      \itemize{
	\item new and faster algorithm to compute monotonic P-splines
	(\code{type = "quad.prog"})
	\item new constraints added for positive and negative spline
	estimates
      }
      \item \code{bbs}
      \itemize{
	\item allows monotone T-splines (experimental)
	\item new argument \code{deriv} to bbs for computing derivatives of
	B-splines
      }
      \item \code{bmrf} can now also handle neighborhood matrices as an
      argument to \code{bnd}
      \item added new families \code{Hurdle} and \code{Multinomial}
      \item \code{boost_control}: added new argument \code{stopintern} for
      internal stopping (based on oobag data) during fitting
      \item All data sets have been moved to the new package set \pkg{TH.data}
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item added new argument \code{which} to \code{variable.names()}
      \item added new method \code{risk} to extract risks
      \item \code{brandom} now checks that a factor is given
      \item speed improvements when updating a model via \code{mod[mstop]}
      \item changed \code{\\dontrun} to \code{\\donttest}
      \item updated references
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item fixed a problem with \code{extract()} of single base-learners
      \item fixed bug in \code{AIC.mboost}: \code{df = "actset"} can
      only be used with glmboost models
      \item fixed package start up messages
      \item fixed a problem in \code{mboost_fit} (when names of
      base-learners were missing)
    }
  }
}

\section{Changes in mboost version 2.2-3 (2013-09-09)}{
  \itemize{
    \item fixed bugs in survival families:
    \itemize{
      \item  \code{offset} in all survival families was based on
      \code{max(survtime)} instead of \code{max(log(survtime))};
      \item  \code{offset} in \code{CoxPH} can't be computed from Cox
      Partial LH as constants are canceled out; Use fixed \code{offset}
      instead;
    }
    \item speed up checking of manual by changing some computations (e.g. reduce
    \code{mstop}) or exclude code from checking via \code{\\dontrun\{\}}
    \item removed dependency on \pkg{ipred} (replaced with \pkg{TH.data})
    \item small improvements in manual
  }
}
\section{Changes in mboost version 2.2-2 (2013-02-08)}{
  \itemize{

    \item \code{bbs(..., center = "spectralDecomp")} computes the spectral
    decomposition of the penalty matrix and the penalized part of the
    design matrix is defined by this decomposition.

    Experiments show that \code{bols(x) + bbs(x, center = "spectralDecomp")}
    is a little better in recovering the true underlying functions than
    the default \code{bols(x) + bbs(x, center = TRUE)} or, equivalently,
    \code{bols(x) + bbs(x, center = "differenceMatrix")}.

    For \code{bbs(x, y, center = TRUE)} or \code{bmrf(x, center = TRUE)},
    the spectral decomposition is (and was) always used.

    \item fixed bug in \code{stabsel}: \code{'...'} was not passed to
         \code{cvrisk} and thus one could not specify options for \code{mclapply}
    \item fixed bug in \code{brandom}: now really use
    \code{contrasts.arg = "contr.dummy"} per default.
    \item removed \code{tests/} folder and \code{.Rout.save} files for
    vignettes from the CRAN release
    \item small improvements in manual
  }
}
\section{Changes in mboost version 2.2-1 (2013-01-14)}{
  \itemize{
    \item included warnings in \code{stabsel()} for better
         guidance of the user:
    \itemize{
         \item A warning is issued if the upper bound for the
         \code{FWER} in stability selection is greater (by a certain
         margin) than the specified bound.
         \item A warning is also issued if \code{mstop} is too small to select
	 \code{q} variables.
    }
    \item improved output of errors and warnings in \code{stabsel}.
    \item suppress the notes from package \pkg{Matrix} about method ambiguity
         (\code{"Note: method with signature ... chosen, ... would also be valid"})
    \item updated manual on base-learners to reflect the change in the default for
         degrees of freedom (additionally, all options are now discussed in a
         separate section of the base-learner manual)
    \item updated vignette \file{mboost_tutorial}
    \item updated \file{mboost_package.Rd}:
         now all important changes since \pkg{mboost} 2.0 are documented there
    \item changed roles of contributors to \code{ctb}
    \item suggested packages are now only used inside \code{if(require(pkg))} statements
    \item changed start up message
  }
}
\section{Changes in mboost version 2.2-0 (2012-11-21)}{
  \itemize{
    \item switch from packages \pkg{multicore} and \pkg{snow} to \pkg{parallel}
    \item changed behavior of \code{bols(x, intercept = FALSE)} when
    \code{x} is a factor:
    \itemize{
      \item now the intercept is simply dropped from the design matrix
      \item coding can be specified as usually for factors.
    }
    \item changed default for \code{options("mboost_dftraceS")} to \code{FALSE}, i.e.,
         degrees of freedom are now computed from smoothing parameter
         as described in B. Hofner, T. Hothorn, T. Kneib, M. Schmid (2011).
    \item changed computation of B-spline basis at the boundaries:
         now also use equidistant knots in the boundaries (per default)
    \item improved \code{plot} function when dealing with spatial plots
         (now builds suitable grid based on the observations if no
         \code{newdata} is given)
    \item increased default number of subsampling replicates in \code{stabsel} to 100
    \item [experimental] \code{bmono()} now implements constraints at the boundaries of
    (monotonic) P-splines
    \item [experimental] added family \code{Gehan()} for rank-based
    estimation of survival models in an accelerated failure time
    framework (contributed by Brent Johnson \email{bajohn3@emory.edu})
  }
}
\section{Changes in mboost version 2.1-3 (2012-09-27)}{
  \itemize{
    \item matrices with one column are now handled as vectors in base-learners
    \item improved manual
    \item fixed error that occurs with R (>= 2.16) due to internal changes in R
  }
}
\section{Changes in mboost version 2.1-2 (2012-02-29)}{
  \itemize{
    \item improved handling of missing values
         (throws warnings and fixed a bug that occurred for missings in the response)
    \item improved manual for the handling of contrasts in bols
    \item added tutorial vignette
    \item updated references
  }
}
\section{Changes in mboost version 2.1-1 (2011-11-28)}{
  \itemize{
    \item new option "mboost_eps" for factor in Demmler-Reinsch orthogonalization
  }
}
\section{Changes in mboost version 2.1-0 (2011-11-15)}{
  \itemize{
    \subsection{Base-learners}{
      \itemize{
	\item added base-learners for smooth monotonic (or convex/concave)
	functions of one or two variables (bmono())
	\item added base-learners for radial basis functions (brad())
	\item added base-learners for Markov random fields (bmrf())
	\item bbs(x, cyclic = TRUE) for cyclic covariates ensures that
	predictions at the boundaries coincide and that the resulting
	function estimate is smoothly joined
	\item bols(x, intercept = FALSE) only reasonable if x is centered.
	A warning is now issued if x is not centered.
	\item changed default for degrees of freedom in bspatial()
	to df = 6
	\item added checks in bbs (and brandom) to ensure that the specified
	degrees of freedom are greater than the range of the (unpenalized)
	null space
	\item bolscw can be mixed with other base-learners (although
	not yet exported and not via the formula interface)
	\item new experimental base-learner \%O\% for smoothing
	matrix-values responses
      }
    }
    \subsection{Families}{
      \itemize{
	\item add Binomial(link = "probit") and general cdf's as link
	functions (experimental)
	\item added new families:
	\itemize{
	  \item AUC() for AUC loss function
	  \item GammaReg() for gamma regression models
	}
      }
    }
    \subsection{Methods}{
      \itemize{
	\item added extract() methods for base-learners and fitted models
	\item added residuals() function to extract residuals from the
        model
	\item improved predict.mboost(): added names where missing and the
        offset as attribute where applicable.
	\item fixed bug in predict() with glmboost.matrix(..., center = TRUE)
	\item coef now also works with tree base-learners
        (returns NULL in this case)
	\item changed coef.gamboost to coef.mboost
	\item various improvements in plot.mboost function
      }
    }

    \subsection{Miscellaneous}{
      \itemize{
	\item changed default in glmboost() to center = TRUE
	\item speed up glmboost() a little bit
	\item changed behavior of cvrisk() if weights are used:
        out-of-bag-risk now weighted according to "weights" as
        specified in call to mboost
	\item added warning if df2lambda is likely to become numerically
        unstable (i.e. in the case of large entries in the design matrix)
	\item improved storage, speed and stability using Matrix technology for
	bols() for factors with many levels and brandom();
	further improvements in base-learners that are combined via \%+\%.
	\item various improvements and fixes in manuals
      }
    }
  }
}
\section{Changes in mboost version 2.0-12 (2011-08-22)}{
  \itemize{
    \item minor bug-fixes to make mboost work with gamboostLSS
    \item replaced writeLines with packageStartupMessage in .onAttach()
    \item replaced partially matched function arguments by full arguments
    \item minor fixes in manuals
  }
}
\section{Changes in mboost version 2.0-11 (2011-03-17)}{
  \itemize{
    \item fix problem in bl_lin when using dense matrices from package "Matrix"
  }
}
\section{Changes in mboost version 2.0-10 (2011-02-20)}{
  \itemize{
    \item add rqss results for India childhood malnutrition data
  }
}
\section{Changes in mboost version 2.0-9 (2010-11-19)}{
  \itemize{
    \item add gbm to Suggests
  }
}
\section{Changes in mboost version 2.0-8 (2010-11-11)}{
  \itemize{
    \item make survival package happy again
  }
}
\section{Changes in mboost version 2.0-7 (2010-09-28)}{
  \itemize{
    \item vignette "mboost" updated
    \item remove problem with R CMD check that occurred on some 64bit systems
  }
}
\section{Changes in mboost version 2.0-6 (2010-05-22)}{
  \itemize{
    \item no not use multicore functionality in R CMD check, really.
  }
}
\section{Changes in mboost version 2.0-5 (2010-05-21)}{
  \itemize{
    \item no not use multicore functionality in R CMD check
  }
}
\section{Changes in mboost version 2.0-4 (2010-04-15)}{
  \itemize{
    \item new vignette "mboost" describing 2.0-x series features
    \item fixed bug in bols(): contrast.arg was ignored if not a named list
         (which is wasn't per default)
    \item added (missing) response functions to families Weibull(),
         Loglog(), Lognormal() and NBinomial()
    \item fixed bug in family CoxPH which occurred with NAs
    \item improvements and corrections in documentation
  }
}
\section{Changes in mboost version 2.0-3 (2010-03-10)}{
  \itemize{
    \item glmboost(..., center = TRUE) now also centers columns of the
         design matrix corresponding to contrasts of factors
         when an intercept term is present leading to faster risk
         minimization in these cases.
    \item coef.glmboost: New argument \code{off2int = TRUE} adds the
         offset to the intercept. In addition, the intercept
         term is now adjusted for centered covariates.
    \item check for infinite residuals in mboost_fit(). Especially
         for family = Poisson(), something like boost_control(nu = 0.01)
         fixes this problem.
    \item "by" (in bols() and bbs()) can now handle factors with more than
         two levels
    \item improved plot.mboost() for varying coefficients
    \item minor improvements in documentation
  }
}
\section{Changes in mboost version 2.0-2 (2010-03-04)}{
  \itemize{
    \item fixed bug in helper function get_index, which caused (in some circumstances)
         wrong handling of factors in gamboost() (spotted by
         Juliane Schaefer <JSchaefer _at_ uhbs.ch>)
    \item reduce memory footprint in blackboost (requires party 0.9-9993)
  }
}
\section{Changes in mboost version 2.0-1 (2010-03-01)}{
  \itemize{
    \item fixed bug in coef( , aggregate = "cumsum"): fraction "nu" was missing
  }
}
\section{Changes in mboost version 2.0-0 (2010-02-01)}{
  \itemize{
    \item generic implementation of component-wise functional gradient
         boosting in \code{mboost_fit}, specialized code for linear,
         additive and interaction models removed
    \item new families available for ordinal, expectile and censored regression
    \item computations potentially based on package Matrix
         (reduces memory usage)
    \item various speed improvements
    \item added interface to extract selected base-learners (selected())
    \item added interface for parallel computations in cvrisk with
         arbitrary packages (e.g. multicore, snow)
    \item added "which" argument in predict and coef functions and improved
         usability of "which" in plot-function. Users can specify "which" as
         numeric value or as a character string
    \item added function cv() to generate matrices for k-fold cross-validation,
         subsampling and bootstrap
    \item new function stabsel() for stability selection with error control
    \item added function model.weights() to extract the weights
    \item added interface to expand model by increasing mstop in
         model[mstop]
    \item alternative definition of degrees of freedom available
    \item Interface changes:
    \itemize{
      \item class definition / Family() arguments changed
      \item changed behavior of subset method (model[mstop]). Object
      is directly altered and not duplicated
      \item argument "center" in bols replaced with "intercept"
      \item argument "z" in base-learners replaced with "by"
      \item bns and bss deprecated;
    }
  }
}
\section{Changes in mboost version 1.1-4 (2009-11-18)}{
  \itemize{
    \item fixed bug in prediction with varying coefficients for binary
         effect modifiers
  }
}
\section{Changes in mboost version 1.1-3 (2009-09-21)}{
  \itemize{
    \item better x-axes in plot.cvrisk and possibility to change xlab
    \item parallel cvrisk on Unix systems only (\pkg{multicore} isn't safe on windows)
    \item included new penalty for ordinal predictors (in bols())
    \item corrected bug in bspatial (centering was not used for Xna)
    \item removed output of dfbase (which is seldom used) in gamboost
    \item changed manual for coef.gamboost
    \item make sure NAs are handled correctly when center = TRUE
         in glmboost
  }
}
\section{Changes in mboost version 1.1-2 (2009-07-21)}{
  \itemize{
    \item better weights and boundary knots handling in bspatial
    \item cvrisk runs in parallel if package \pkg{multicore} is available
    \item errors removed and minor improvements in the manuals
    \item center = TRUE in glmboost did only apply to numeric (not integer)
         predictors
    \item for safety reasons: na.action = na.omit again
         (causes slight changes in wpbc3 example)
  }
}
\section{Changes in mboost version 1.1-1 (2009-04-21)}{
  \itemize{
    \item new quantile regression facilities.
    \item fix problem with bbs base-learner and cvrisk
  }
}
\section{Changes in mboost version 1.1-0 (2009-03-27)}{
  \itemize{
    \item bbs instead of bss is the default base-learner in gamboost
    \item make sure bbs with weights and expanded observations
         returns numerically the very same results
    \item btree can now deal with multiple variables
    \item new gMDL criterion (contributed by Zhu Wang <zhu.wang@yale.edu>)
    \item make survival package happy again
  }
}
\section{Changes in mboost version 1.0-6 (2009-01-09)}{
  \itemize{
    \item bols allows to specify non-default contrasts.
  }
}
\section{Changes in mboost version 1.0-5 (2008-12-02)}{
  \itemize{
    \item remove experimental memory optimization steps
  }
}
\section{Changes in mboost version 1.0-4 (2008-11-12)}{
  \itemize{
    \item negative gradient of GaussClass() was wrong, spotted by
         Kao Lin <linkao@picb.ac.cn>
  }
}
\section{Changes in mboost version 1.0-3 (2008-11-07)}{
  \itemize{
    \item Date was malformed in DESCRIPTION
  }
}
\section{Changes in mboost version 1.0-2 (2008-11-05)}{
  \itemize{
    \item improved memory footprint in gamboost() and cvrisk()
    \item option to suppress saving of ensembles added to
         boost_control()
    \item bbs(), bns(), bspatial(): default number of knots changed to
         a fixed value (= 20)
    \item changed default for grid (now uses all iterations) in
         cvrisk() and changed plot.cvrisk()
    \item bols: works now for factors and can be set-up
         to use Ridge-estimation. Intercept can be omitted
         now (via center = TRUE).
    \item new btree() base-learner for gamboost() available
    \item fix inconsistencies in regression tests
    \item add coef.gamboost
    \item new generic \code{survFit}
    \item cosmetics for trace = TRUE
  }
}
\section{Changes in mboost version 1.0-1 (2007-12-09)}{
  \itemize{
    \item inst/mboost_Bioinf.R was missing from mboost 1.0-0
  }
}
\section{Changes in mboost version 1.0-0 (2007-11-13)}{
  \itemize{
    \item documentation updates
  }
}
\section{Changes in mboost version 0.9-0}{
  \itemize{
    \item tests update and release the new version on CRAN
    \item predict(..., allIterations = TRUE) returns the matrix
         of predictors for all boosting iterations
  }
}
\section{Changes in mboost version 0.6-2}{
  \itemize{
    \item move \pkg{mboost} to R-forge
    \item improvements in \code{gamboost}:
         \item P-splines as base learners available
         \item new formula interface for specifying the base learner
         \item new plot.gamboost
    \item add the number of selected variables as degrees of freedom
         (as mentioned in the discussion of Hastie to Buehlmann & Hothorn)
    \item status information during fitting is now available via
         boost_control(trace = TRUE) but is switched off by default
    \item acknowledge contributions by Thomas Kneib and Matthias Schmid
         in DESCRIPTION
  }
}
\section{Changes in mboost version 0.6-1}{
  \itemize{
    \item gamboost() now allows for user-specified base learners
         via the formula interface
    \item gamboost.matrix(x = x, ...) requires colnames being set
         for \code{x}
    \item na.action = na.omit fix for g\{al\}mboost()
  }
}
\section{Changes in mboost version 0.5-8 (2007-05-31)}{
  \itemize{
    \item gamboost(..., weights = w) was broken
  }
}
\section{Changes in mboost version 0.5-7 (2007-05-30)}{
  \itemize{
    \item extract response correctly in fitted.blackboost
    \item hatvalues (and thus AICs) for GLMs with centering of
         covariates may have been wrong since version 0.5-0
    \item add paper examples to tests
  }
}
\section{Changes in mboost version 0.5-6 (2007-05-07)}{
  \itemize{
    \item fix Rd problems
  }
}
\section{Changes in mboost version 0.5-5 (2007-04-25)}{
  \itemize{
    \item \code{westbc} regenerated
    \item LazyLoad: yes (no SaveImage: yes)
  }
}
\section{Changes in mboost version 0.5-4 (2007-04-18)}{
  \itemize{
    \item plot() method for \code{glmboost} objects visualizing the
         coefficient path (feature request by Axel Benner <benner@dkfz.de>).
    \item predict(newdata = <matrix>) was broken for gamboost(),
         thanks to Max Kuhn <Max.Kuhn@pfizer.com> for spotting this.
  }
}
\section{Changes in mboost version 0.5-3 (2007-03-23)}{
  \itemize{
    \item predict() for gamboost(..., dfbase = 1) was not working correctly
    \item small performance and memory improvements for glmboost()
  }
}
\section{Changes in mboost version 0.5-2 (2007-02-28)}{
  \itemize{
    \item some performance improvements for \code{glmboost()}
    \item blackboost() is now generic with formula and x, y interface
    \item plot() method for cvrisk() and AIC() output now allows for ylim
         specification without troubles
  }
}
\section{Changes in mboost version 0.5-1 (2007-02-02)}{
  \itemize{
    \item depends party 0.9-9
  }
}
\section{Changes in mboost version 0.5-0 (2007-01-30)}{
  \itemize{
    \item new \code{baselearner} argument to \code{gamboost} allowing to
         specify difference component-wise base-learners to
         be used. Currently implemented: "ssp" for smoothing splines
         (default), "bsp" for B-splines and "ols" for linear models.
         The latter two haven't been tested yet.
    \item The \code{dfbase} arguments now applies to each covariate and
         no longer to each column of the design matrix.
    \item cvrisk() for blackboost() was broken, totally :-(
    \item centered covariates were returned by glmboost() and gamboost()
    \item Poisson() used an incorrect offset
    \item check for y being positive counts when family = "Poisson()"[B
    \item checks for Poisson() logLik() and AIC() methods
    \item fire a warning when all u > 0 or u < 0
    \item update vignette \file{mboost_illustrations}
  }
}
\section{Changes in mboost version 0.4-17 (2007-01-15)}{
  \itemize{
    \item fix problem with \code{dfbase} in \code{gamboost}, spotted by
         Karin Eckel <Karin.Eckel@imbe.imed.uni-erlangen.de>
  }
}
\section{Changes in mboost version 0.4-16 (2007-01-12)}{
  \itemize{
    \item work around stats4:::AIC
  }
}
\section{Changes in mboost version 0.4-15 (2006-12-06)}{
  \itemize{
    \item fix plot problems in plot.cvrisk
    \item allow for centering of the numerical covariates in glmboost and
         gamboost
  }
}
\section{Changes in mboost version 0.4-14 (2006-10-27)}{
  \itemize{
    \item AIC(..., "classical") is now faster for non-Gaussian families
  }
}
\section{Changes in mboost version 0.4-13 (2006-10-04)}{
  \itemize{
    \item predict(..., newdata) can take a matrix now
  }
}
\section{Changes in mboost version 0.4-12 (2006-09-13)}{
  \itemize{
    \item predict(<blackboost-object>, type = "response") did not return
         factors when the response was actually a factor
    \item report offset in print methods
    \item add offset attribute to coef.glmboost
  }
}
\section{Changes in mboost version 0.4-11 (2006-09-07)}{
  \itemize{
    \item add \code{contrasts.arg} argument to \code{glmboost.formula}
    \item more meaningful default for \code{grid} in \code{cvrisk}
    \item R-2.4.0 fixes
  }
}
\section{Changes in mboost version 0.4-10 (2006-08-30)}{
  \itemize{
    \item add checks for CoxPH (against coefficients and logLik of CoxPH)
    \item add weights to CoxPH
    \item the ngradient function in Family objects needs to implement
         arguments (y, f, w), not just (y, f)
    \item check for meaningful class of the response for some families
  }
}
\section{Changes in mboost version 0.4-9 (2006-07-17)}{
  \itemize{
    \item some small speed improvements in \code{gamboost}
    \item handle factors in \code{gamboost} properly (via a linear model)
    \item the dfbase argument can take a vector now (in \code{gamboost})
    \item update and improve entries in DESCRIPTION
    \item documentation updates
  }
}
\section{Changes in mboost version 0.4-8 (2006-07-05)}{
  \itemize{
    \item Huber() is \sQuote{Huber Error}, not \sQuote{Huber Absolute Error}
    \item added \code{CoxPH} family object for fitting Cox models
    \item remove inst/LaTeX
    \item use NROW / NCOL more often (now that \code{y} may be a \code{Surv} object)
    \item implement \code{cvrisk}, a general cross-validation function for the
         empirical risk and a corresponding plot method
    \item unify risk computations in all three fitting functions
    \item unify names for \code{gb} objects
    \item allow for out-of-bag risk computations
    \item some cosmetics
    \item update keywords in Rd-files
    \item risk was always 0 in Huber()@risk when d was chosen adaptively
    \item pData(westbc)$nodal.y has levels \code{negative} and \code{positive}
         (lymph node status)
  }
}
\section{Changes in mboost version 0.4-7 (2006-06-19)}{
  \itemize{
    \item add src/Makevars (required for Windows builds)
    \item make sure objects that are modified at C-level get _copied_ in
         \code{blackboost}
  }
}
\section{Changes in mboost version 0.4-6 (2006-06-14)}{
  \itemize{
    \item some minor \code{codetools} fixes: removed unused variables
         and an out-dated function
    \item add \code{codetools} checks to regression tests
    \item fix xlab in plot.gbAIC
    \item \pkg{mboost} version 0.4-5 published on CRAN 2006-06-13
  }
}
