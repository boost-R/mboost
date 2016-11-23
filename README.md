mboost
======

[![Build Status (Linux)](https://travis-ci.org/boost-R/mboost.svg?branch=master)](https://travis-ci.org/boost-R/mboost)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/5mkvicgin1j6pfc6/branch/master?svg=true)](https://ci.appveyor.com/project/hofnerb/mboost-h73a1/branch/master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mboost)](https://CRAN.R-project.org/package=mboost)
[![Coverage Status](https://coveralls.io/repos/github/boost-R/mboost/badge.svg?branch=master)](https://coveralls.io/github/boost-R/mboost?branch=master)
[![](http://cranlogs.r-pkg.org/badges/mboost)](https://CRAN.R-project.org/package=mboost)

`mboost` implements boosting algorithms for fitting generalized linear, additive and interaction models
to potentially high-dimensional data.

## Using mboost

For installation instructions see below.

Instructions on how to use `mboost` can be found in various places:
- Have a look at the tutorials:
  - [mboost tutorial](https://CRAN.R-project.org/package=mboost/vignettes/mboost_tutorial.pdf)
  - [mboost 2.0](https://CRAN.R-project.org/package=mboost/vignettes/mboost.pdf)
- Visit the [project homepage](http://mboost.R-forge.R-project.org/) and see further tutorials and references there.

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/boost-R/mboost/issues).

## Installation Instructions

- Current version (from CRAN):
  ```r
  install.packages("mboost")
  ```

- Latest **patch version** (patched version of CRAN package; under development) from GitHub:
  ```r
  library("devtools")
  install_github("boost-R/mboost")
  library("mboost")
  ```

- Latest **development version** (version with new features; under development) from GitHub:
  ```r
  library("devtools")
  install_github("boost-R/mboost", ref = "devel")
  library("mboost")
  ```

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  ```r
  install.packages("devtools")
  ```

