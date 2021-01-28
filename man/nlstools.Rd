\name{nlstools}
\alias{nlstools}
\alias{preview}
\alias{plotfit}
\alias{overview}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Nonlinear least squares fit }

\description{
Tools to help the fit of nonlinear models with nls
}

\usage{
preview (formula, data, start, variable = 1)
plotfit (x, smooth = FALSE, variable = 1, xlab = NULL, ylab = NULL, 
         pch.obs = 1, pch.fit = "+", lty = 1, lwd = 1, col.obs = "black", 
         col.fit = "red", ...)
overview (x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{formula}{ formula of a non-linear model }
	\item{data}{ a data frame with header matching the variables given in the formula }
	\item{start}{ a list of parameter starting values which names match the parameters given in the formula }
	\item{variable}{ index of the variable to be plotted against the predicted values; default is the first independent variable as it appears in the orginal dataset }
	\item{x}{ an object of class 'nls' }
	\item{smooth}{ a logical value, default is FALSE. If smooth is TRUE, a plot of observed values is plotted as a function of 1000 values continuously taken in the range interval [min(variable),max(variable)]. This option can only be used if the number of controlled variables is 1. } 
	\item{xlab}{ X-label }
	\item{ylab}{ Y-label }
	\item{pch.obs}{ type of point of the observed values }
	\item{pch.fit}{ type of point of the fitted values (not applicable if smooth=TRUE)}
	\item{lty}{ type of line of the smoothed fitted values (if smooth=TRUE) }
	\item{lwd}{ thickness of line of the smoothed fitted values (if smooth=TRUE) }
	\item{col.obs}{ color of the observed points }
	\item{col.fit}{ color of the fitted values }
	\item{...}{ further arguments passed to or from other methods }

}
\details{
	The function \code{preview} helps defining the parameter starting values prior fitting the model. It provides a superimposed plot of observed (circles) and predicted (crosses) values of the dependent variable versus one of the independent variables with the model evaluated at the starting values of the parameters. The function \code{overview} returns the parameters estimates, their standard errors as well as their asymptotic confidence intervals and the correlation matrix (alternately, the function \code{confint} provides better confidence interval estimates whenever it converges). \code{plotfit} displays a superimposed plot of the dependent variable versus one the independent variables together with the fitted model.
}

%\value{ }

\seealso{
\code{nls} in the \code{stats} library and \code{confint.nls} in the package \code{MASS}
}

\references{ 
Baty F, Ritz C, Charles S, Brutsche M, Flandrois J-P, Delignette-Muller M-L (2015). A Toolbox for Nonlinear Regression in R: The Package nlstools. \emph{Journal of Statistical Software}, \bold{66}(5), 1-21.\cr\cr
Bates DM and Watts DG (1988) Nonlinear regression analysis and its applications. Wiley, Chichester, UK.
}
\author{ 
Florent Baty \email{florent.baty@gmail.com}\cr
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr}
}

%\note{  }

\examples{
formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
                        (VO2rest + (VO2peak - VO2rest) * 
                        (1 - exp(-(t - 5.883) / mu))))
preview(formulaExp, O2K, list(VO2rest = 400, VO2peak = 1600, mu = 1))
O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, 
                mu = 1), data = O2K)
overview(O2K.nls1)
plotfit(O2K.nls1, smooth = TRUE)
}

\keyword{ nonlinear }% at least one, from doc/KEYWORDS
