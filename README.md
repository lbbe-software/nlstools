[![CRAN_Release_Badge](http://www.r-pkg.org/badges/version-ago/nlstools)](https://cran.r-project.org/package=nlstools)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/nlstools)](https://cran.r-project.org/package=nlstools)
[![R-CMD-check](https://github.com/lbbe-software/nlstools/workflows/R-CMD-check/badge.svg)](https://github.com/lbbe-software/nlstools/actions)

---------------------------

**Please note! Since February 2024, this repository has belonged to the *lbbe-software* organization.**
To avoid confusion, we strongly recommend updating any existing local clones to point to the new 
repository URL. You can do this by using `git remote` on the command line:

`git remote set-url origin git@github.com:lbbe-software/nlstools.git`

or 

`git remote set-url origin https://github.com/lbbe-software/nlstools.git`

---------------------------

The stable version of `nlstools` can be installed from CRAN using:
```r
install.packages("nlstools")
```

The development version of `nlstools` can be installed from GitHub (`remotes` needed):
```r
if (!requireNamespace("remotes", quietly = TRUE))
   install.packages("remotes")
   
remotes::install_github("lbbe-software/nlstools")
``` 

Finally load the package in your current R session with the following R command:
```r
library(nlstools)
```
