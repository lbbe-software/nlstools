#' Enzyme kinetics
#' 
#' Enzyme kinetics
#' 
#' 
#' @format A data frame with 8 observations on the following 2 variables.
#' \describe{
#'  \item{conc}{a numeric vector} 
#'  \item{rate}{a numeric vector} 
#' }
#' @source Cedergreen, N. and Madsen, T. V. (2002) Nitrogen uptake by the
#' floating macrophyte \emph{Lemna minor}, \emph{New Phytologist}, \bold{155},
#' 285--292.
#' @keywords datasets
"L.minor"

#' Michaelis Menten data sets
#' 
#' Michaelis Menten data sets
#' 
#' 
#' @aliases vmkmki michaelisdata
#' @format \code{vmkm} is a data frame with 2 columns (S: concentration of
#' substrat, v: reaction rate)\cr 
#' \code{vmkmki} is a data frame with 3 columns
#' (S: concentration of substrat, I: concentration of inhibitor, v: reaction
#' rate)
#' @source These datasets were provided by the French research unit INRA
#' UMR1233.
#' @keywords datasets
#' @examples
#' 
#' data(vmkm)
#' data(vmkmki)
#' plot(vmkm)
#' plot(vmkmki)
#' 
"vmkm"

#' Oxygen kinetics during 6-minute walk test data set
#' 
#' Oxygen uptake kinetics during a 6-minute walking test in a patient with
#' pulmonary disease. The first 5.83 minutes correspond to the resting phase
#' prior to exercise.
#' 
#' 
#' @format \code{O2K} is a data frame with 2 columns (t: time, VO2: oxygen
#' uptake)\cr
#' @source This data set was provided by the Cantonal Hospital St. Gallen,
#' Switzerland.
#' @keywords datasets
#' @examples
#' 
#' data(O2K)
#' plot(O2K)
#' 
"O2K"



