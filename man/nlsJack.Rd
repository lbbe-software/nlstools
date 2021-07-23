\name{nlsJack}
\alias{nlsJack}
\alias{plot.nlsJack}
\alias{print.nlsJack}
\alias{summary.nlsJack}
\title{Jackknife resampling}
\description{Jackknife resampling}

\usage{
  nlsJack (nls)
  \method{plot}{nlsJack} (x, mfr = c(nrow(x$reldif),1), ask = FALSE, \dots)
  \method{print}{nlsJack} (x, \dots)
  \method{summary}{nlsJack} (object, \dots)
}

\arguments{
	\item{nls}{an object of class 'nls'}
	\item{x, object}{an object of class 'nlsJack'}
	\item{mfr}{layout definition, default is k rows (k: number of parameters) and 1 column}
	\item{ask}{if TRUE, draw plot interactively}
	\item{...}{further arguments passed to or from other methods}
}

\value{ 
	\code{nlsJack} returns a list with 7 objects:
	\item{estijack}{a data frame with jackknife estimates and bias}
	\item{coefjack}{a data frame with the parameter estimates for each jackknife sample}
	\item{reldif}{a data frame with the absolute relative difference (in percent of the estimate) of the estimates with and without each observation}
	\item{dfb}{a data frame with dfbetas for each parameter and each observation}
	\item{jackCI}{a data frame with jackknife confidence intervals}
	\item{rse}{a vector with residual standard error for each jackknife sample}
	\item{rss}{residual a vector with residual sum of squares for each jackknife sample}
}

\details{
	A jackknife resampling procedure is performed. Each observation is sequentially removed from the initial data set using a leave-one-out strategy. A data set with \emph{n} observations provides thus \emph{n} resampled data sets of \emph{n-1} observations. The jackknife estimates with confidence intervals are calculated as described by Seber and Wild (1989) from the results of \emph{n} new fits of the model on the \emph{n} jackknife resampled data sets. The leave-one-out procedure is also employed to assess the influence of each observation on each parameter estimate. An observation is empirically defined as influential for one parameter if the difference between the estimate of this parameter with and without the observation exceeds twice the standard error of the estimate divided by \emph{sqrt(n)}. This empirical method assumes a small curvature of the nonlinear model. For each parameter, the absolute relative difference (in percent of the estimate) of the estimates with and without each observation is plotted. An asterisk is plotted for each influential observation.
}

\references{ 
	Seber GAF, Wild CJ (1989) Nonlinear regression. Wiley, New York.\cr\cr
}

\author{Florent Baty, Marie-Laure Delignette-Muller}

\examples{
formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
                        (VO2rest + (VO2peak - VO2rest) * 
                        (1 - exp(-(t - 5.883) / mu))))
O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, mu = 1), 
               data = O2K)
O2K.jack1 <- nlsJack(O2K.nls1)
plot(O2K.jack1)
summary(O2K.jack1)
}

\keyword{nonlinear}