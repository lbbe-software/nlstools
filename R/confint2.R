#' Confidence intervals in nonlinear regression
#' 
#' Produces confidence intervals for the parameters in nonlinear regression
#' model fit. The intervals can either be based large sample results or on
#' profiling.
#' 
#' The profiling used is the method \code{\link[MASS]{confint.nls}.}
#' 
#' @param object object of class \code{\link{nls}}.
#' @param parm a vector character strings with names of the parameter for which
#' to calculate confidence intervals (by default all parameters).
#' @param level the confidence level required.
#' @param method method to be used: "asympotic" for large sample and "profile"
#' for profiling approach.
#' @param \dots additional argument(s) to pass on the method doing the
#' profiling.
#' 
#' @importFrom stats coef qt df.residual vcov profile confint
#' 
#' @return A matrix with columns giving lower and upper confidence limits for
#' each parameter.
#' @author Christian Ritz
#' @keywords models nonlinear
#' @examples
#' 
#' 
#' L.minor.m1 <- nls(rate ~ Vm*conc/(K+conc), data = L.minor, start = list(K=20, Vm=120))
#' 
#' confint2(L.minor.m1)
#' 
#' confint2(L.minor.m1, "K")
#' 
#' 
#' @export confint2
"confint2" <- function(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...)
{
    method <- match.arg(method)

    format.perc <- function(probs, digits)
      ## Not yet exported, maybe useful in other contexts:
      ## quantile.default() sometimes uses a version of it
      paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
            "%")
    
    ## Taken from confint.nls
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- seq_along(pnames)
    if (is.numeric(parm)) 
        parm <- pnames[parm]
    
    ## Taken from confint.default and modified slightly to use t-distribution
    asCI <- function(object, parm, level)
    {
        a <- (1 - level)/2
        a <- c(a, 1 - a)
#        pct <- stats:::format.perc(a, 3)
        pct <- format.perc(a, 3)
        fac <- qt(a, df.residual(object))
        
        parmInd <- match(parm, pnames)
        ci <- array(NA, dim = c(length(parmInd), 2), dimnames = list(parm, pct))
        ses <- sqrt(diag(vcov(object)))[parmInd]
        ci[] <- cf[parmInd] + ses %o% fac
        ci
    }

    ## Taken from confint.nls
    asProf <- function(object, parm, level)
    {
        message("Waiting for profiling to be done...")
        utils::flush.console()
        object <- profile(object, which = parm, alphamax = (1 - level)/4)
        confint(object, parm = parm, level = level, ...)    
    }

    switch(method, asymptotic = asCI(object, parm, level), profile = asProf(object, parm, level))
}

