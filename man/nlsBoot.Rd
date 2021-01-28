\name{nlsBoot}
\alias{nlsBoot}
\alias{plot.nlsBoot}
\alias{print.nlsBoot}
\alias{summary.nlsBoot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Bootstrap resampling }

\description{
Bootstrap resampling
}

\usage{
nlsBoot (nls, niter = 999)
\method{plot}{nlsBoot} (x, type = c("pairs", "boxplot"), 
    mfr = c(ceiling(sqrt(ncol(x$coefboot))), 
    ceiling(sqrt(ncol(x$coefboot)))), 
    ask = FALSE, \dots)
\method{print}{nlsBoot} (x, \dots)
\method{summary}{nlsBoot} (object, \dots)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
    \item{nls}{ an object of class 'nls' }
    \item{niter}{ number of iterations  }
    \item{x, object}{ an object of class 'nlsBoot' }
    \item{type}{ type of representation (options are "pairs" or "boxplot")}
    \item{mfr}{ layout definition (number of rows and columns in the graphics device) }
    \item{ask}{ if TRUE, draw plot interactively }
    \item{...}{ further arguments passed to or from other methods }
}

\details{
    Non-parametric bootstrapping is used. Mean centered residuals are bootstrapped. By default, 999 resampled data sets are created from which parameter estimates are obtained by fitting the model on each of these data sets. Whenever the fit fails to converge, a flag reports the number of non-convergences. If the fitting procedure fails to converge in more than 50\% of the cases, the procedure is interrupted with a flag and no result is given. The function \code{summary} returns the bootstrap estimates (mean and std. dev. of the bootstrapped estimates) and the median and 95 percent confidence intervals (50, 2.5, and 97.5 percentiles of the bootstrapped estimates). The bootstrapped estimate distributions can be visualized using the function \code{plot.nlsBoot} either by plotting the bootstrapped sample for each pair of parameters or by displaying the boxplot representation of the bootstrapped sample for each parameter. Notice that \code{nlsBoot} does not currently handle transformed dependent variables specified in the left side of the \code{nls} formula.  
}

\value{ 
    \code{nlsBoot} returns a list of 4 objects:
    \item{ coefboot }{ contains the bootstrap parameter estimates }
    \item{ bootCI }{ contains the bootstrap medians and the bootstrap 95\% confidence intervals }
    \item{ estiboot }{ contains the means and std. errors of the bootstrap parameter estimates  }
    \item{ rse }{ is the vector of bootstrap residual errors }
}

\references{ 
Bates DM and Watts DG (1988) Nonlinear regression analysis and its applications. Wiley, Chichester, UK.\cr\cr
Huet S, Bouvier A, Poursat M-A, Jolivet E (2003) Statistical tools for nonlinear regression: a practical guide with S-PLUS and R examples. Springer, Berlin, Heidelberg, New York.
}

\author{ 
Florent Baty \email{florent.baty@gmail.com}\cr
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr}
}

%\note{  }

\examples{
formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2res + (t > 5.883) * 
                        (VO2res + (VO2peak - VO2res) * 
                        (1 - exp(-(t - 5.883) / mu))))
O2K.nls1 <- nls(formulaExp, start = list(VO2res = 400, VO2peak = 1600, 
                mu = 1), data = O2K)
O2K.boot1 <- nlsBoot(O2K.nls1, niter = 200)
plot(O2K.boot1)
plot(O2K.boot1, type = "boxplot", ask = FALSE)
summary(O2K.boot1)
  }

\keyword{ nonlinear }% at least one, from doc/KEYWORDS
