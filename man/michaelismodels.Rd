\name{michaelismodels}
\alias{michaelismodels}
\alias{michaelis}
\alias{compet_mich}
\alias{non_compet_mich}
\title{ Michaelis-Menten model and derived equations to model competitive and non-competitive inhibition }
\description{Formula of Michaelis-Menten model commonly used to describe enzyme kinetics, and derived formulas taking into account the 
effect of a competitive or a non-competitive inhibitor
}

\usage{
  michaelis
  compet_mich
  non_compet_mich
}

\details{
    These models describe the evolution of the reaction rate (v) as a function
    of the concentration of substrate (S) and the concentration of inhibitor (I) for \code{compet_mich} and \code{non_compet_mich}.
    
    \code{michaelis} is the classical Michaelis-Menten model (Dixon, 1979)  with two parameters (Km, Vmax) :
    \deqn{v = \frac{S}{S+K_m} V_{max}}{v = S/(S+Km)*Vmax} \cr\cr
    \code{compet_mich} is the Michaelis-Menten derived model with three parameters (Km, Vmax, Ki), describing
    a competitive inhibition : \deqn{v = \frac{S}{S + K_m (1+\frac{I}{K_i})} V_{max}}{v = S/(S + Km*(1+I/Ki) ) * Vmax} \cr\cr
    \code{non_compet_mich} is the Michaelis-Menten derived model with three parameters (Km, Vmax, Ki), describing
    a non-competitive inhibition : \deqn{v = \frac{S}{(S+K_m)(1+\frac{I}{Ki})} V_{max}}{v = S/( (S + Km)*(1+I/Ki) ) * Vmax} \cr\cr
}

\value{A formula}

\references{ 
  Dixon M and Webb EC (1979) \emph{Enzymes}, Academic Press, New York.
}

\author{Florent Baty, Marie-Laure Delignette-Muller}

\examples{
# Example 1

data(vmkm)
nls1 <- nls(michaelis,vmkm,list(Km=1,Vmax=1))
plotfit(nls1, smooth = TRUE)

# Example 2

data(vmkmki)
def.par <- par(no.readonly = TRUE)
par(mfrow = c(2,2))

nls2_c <- nls(compet_mich, vmkmki, list(Km=1,Vmax=20,Ki=0.5))
plotfit(nls2_c, variable=1)
overview(nls2_c)
res2_c <- nlsResiduals(nls2_c)
plot(res2_c, which=1)

nls2_nc <- nls(non_compet_mich, vmkmki, list(Km=1, Vmax=20, Ki=0.5))
plotfit(nls2_nc, variable=1)
overview(nls2_nc)
res2_nc <- nlsResiduals(nls2_nc)
plot(res2_nc, which=1)

par(def.par)
}

\keyword{models}