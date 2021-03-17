#' Michaelis-Menten model and derived equations to model competitive and
#' non-competitive inhibition
#' 
#' Formula of Michaelis-Menten model commonly used to describe enzyme kinetics,
#' and derived formulas taking into account the effect of a competitive or a
#' non-competitive inhibitor
#' 
#' 
#' These models describe the evolution of the reaction rate (v) as a function
#' of the concentration of substrate (S) and the concentration of inhibitor (I)
#' for \code{compet_mich} and \code{non_compet_mich}.
#' 
#' \code{michaelis} is the classical Michaelis-Menten model (Dixon, 1979) with
#' two parameters (Km, Vmax) : \deqn{v = \frac{S}{S+K_m} V_{max}}{v =
#' S/(S+Km)*Vmax} \cr\cr 
#' 
#' \code{compet_mich} is the Michaelis-Menten derived
#' model with three parameters (Km, Vmax, Ki), describing a competitive
#' inhibition : \deqn{v = \frac{S}{S + K_m (1+\frac{I}{K_i})} V_{max}}{v = S/(S
#' + Km*(1+I/Ki) ) * Vmax} \cr\cr 
#' 
#' \code{non_compet_mich} is the
#' Michaelis-Menten derived model with three parameters (Km, Vmax, Ki),
#' describing a non-competitive inhibition : \deqn{v =
#' \frac{S}{(S+K_m)(1+\frac{I}{Ki})} V_{max}}{v = S/( (S + Km)*(1+I/Ki) ) *
#' Vmax} \cr\cr
#' 
#' @aliases michaelismodels michaelis compet_mich non_compet_mich
#' @return A formula
#' @author Florent Baty \email{florent.baty@@gmail.com}\cr Marie-Laure
#' Delignette-Muller \email{ml.delignette@@vetagro-sup.fr}
#' @references Dixon M and Webb EC (1979) \emph{Enzymes}, Academic Press, New
#' York.
#' @keywords models
#' @examples
#' 
#' 
#' # Example 1
#' 
#' data(vmkm)
#' nls1 <- nls(michaelis,vmkm,list(Km=1,Vmax=1))
#' plotfit(nls1, smooth = TRUE)
#' 
#' # Example 2
#' 
#' data(vmkmki)
#' def.par <- par(no.readonly = TRUE)
#' par(mfrow = c(2,2))
#' 
#' nls2_c <- nls(compet_mich, vmkmki, list(Km=1,Vmax=20,Ki=0.5))
#' plotfit(nls2_c, variable=1)
#' overview(nls2_c)
#' res2_c <- nlsResiduals(nls2_c)
#' plot(res2_c, which=1)
#' 
#' nls2_nc <- nls(non_compet_mich, vmkmki, list(Km=1, Vmax=20, Ki=0.5))
#' plotfit(nls2_nc, variable=1)
#' overview(nls2_nc)
#' res2_nc <- nlsResiduals(nls2_nc)
#' plot(res2_nc, which=1)
#' 
#' par(def.par)
#' 
#' @export
#' 
michaelis <- as.formula(v ~ S/(S+Km) * Vmax)

#' @rdname michaelis
#' @export
compet_mich <- as.formula(v ~ S/(S + Km*(1+I/Ki) ) * Vmax)

#' @rdname michaelis
#' @export
non_compet_mich <- as.formula(v ~ S/( (S + Km)*(1+I/Ki) ) * Vmax)
