mboost
======

[![Build Status](https://travis-ci.org/hofnerb/mboost.svg?branch=master)](https://travis-ci.org/hofnerb/mboost)
[![](http://cranlogs.r-pkg.org/badges/mboost)](http://cran.rstudio.com/web/packages/mboost/index.html)

`mboost` implements boosting algorithms for fitting generalized linear, additive and interaction models
to potentially high-dimensional data.

This [GitHub repository](https://github.com/hofnerb/mboost) is essentially just
a copy of the r-forge repository which is hosted at
[R-forge](https://r-forge.r-project.org/projects/mboost).

## Using mboost

For installation instructions see below.

Instructions on how to use `mboost` can be found in various places:
- Have a look at the tutorials:
  - [mboost tutorial](http://cran.r-project.org/web/packages/mboost/vignettes/mboost_tutorial.pdf)
  - [mboost 2.0](http://cran.r-project.org/web/packages/mboost/vignettes/mboost.pdf)
- Visit the [project homepage](http://mboost.r-forge.r-project.org/) and see further tutorials and references there.

## Reporting Issues (etc.)

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/hofnerb/mboost/issues).

## Installation Instructions

- Current version (from CRAN):
  ```
  install.packages("mboost")
  ```

- Latest **patch version** (under development) from GitHub:
  ```
  library("devtools")
  install_github("hofnerb/mboost/pkg/mboostPatch")
  library("mboost")
  ```

- Latest **development version** from GitHub:
  ```
  library("devtools")
  install_github("hofnerb/mboost/pkg/mboostDevel")
  library("mboostDevel")
  ```

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  ```
  install.packages("devtools")
  ```

- Alternatively, both the current patch and development versions of `mboost` (or `mboostDevel` respectively)
  can be downloaded from R-forge if it was successfully built:
  ```
  install.packages("mboost", repos = "http://r-forge.r-project.org")
  ## or
  install.packages("mboostDevel", repos = "http://r-forge.r-project.org")
  ```
  However, currently these builds often don't succeed and furthermore are only available
  for recent versions of R.

[inst]: inst

