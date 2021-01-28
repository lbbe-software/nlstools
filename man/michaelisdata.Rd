\name{michaelisdata}
\alias{vmkm}
\alias{vmkmki}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Michaelis Menten data sets }

\description{
Michaelis Menten data sets
}

\usage{
data(vmkm)
data(vmkmki)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{vmkm} is a data frame with 2 columns (S: concentration of substrat, v: reaction rate)\cr
\code{vmkmki} is a data frame with 3 columns (S: concentration of substrat, I: concentration of inhibitor, v: reaction rate)
}

\source{ These datasets were provided by the French research unit INRA UMR1233. }

%\references{  }

\examples{
data(vmkm)
data(vmkmki)
plot(vmkm)
plot(vmkmki)
}

\keyword{ datasets }% at least one, from doc/KEYWORDS
