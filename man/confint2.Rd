\name{confint2}

\alias{confint2}

\title{Confidence intervals in nonlinear regression}

\description{
  Produces confidence intervals for the parameters in nonlinear regression model fit. The intervals
  can either be based large sample results or on profiling.
}

\usage{
  confint2(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...)
}

\arguments{
  \item{object}{object of class \code{\link{nls}}.}
  \item{parm}{a vector character strings with names of the parameter for which to calculate 
  confidence intervals (by default all parameters).}
  \item{level}{the confidence level required.}
  \item{method}{method to be used: "asympotic" for large sample and "profile" for profiling approach.}
  \item{\dots}{additional argument(s) to pass on the method doing the profiling.}
}

\details{
  The profiling used is the method \code{\link[MASS]{confint.nls}.} 
}

\value{
  A matrix with columns giving lower and upper confidence limits for each parameter.
}

%\references{}

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

L.minor.m1 <- nls(rate ~ Vm*conc/(K+conc), data = L.minor, start = list(K=20, Vm=120))

confint2(L.minor.m1)

confint2(L.minor.m1, "K")

}

\keyword{models}
\keyword{nonlinear}

%\concept{}

