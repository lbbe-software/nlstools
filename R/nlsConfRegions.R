#' Confidence regions
#' 
#' Draws parameter values in the Beale's 95 percent unlinearized confidence
#' region
#' 
#' A sample of points in the 95 percent confidence region is computed according
#' to Beale's criterion (Beale, 1960). This region is also named the joint
#' parameter likelihood region (Bates and Watts, 1988). The method used
#' consists in a random sampling of parameters values in a hypercube centered
#' on the least squares estimate and rejecting the parameters values whose
#' residual sum of squares do not verify the Beale criterion. The confidence
#' region is plotted by projection of the sampled points in each plane defined
#' by a couple of parameters. Bounds of the hypercube in which random values of
#' parameters are drawn may be plotted in order to check if the confidence
#' region was totally included in the hypercube defined by default. If not the
#' hypercube should be expanded in order to obtain the full confidence region
#' 
#' @aliases nlsConfRegions plot.nlsConfRegions print.nlsConfRegions
#' @param nls an object of class 'nls'
#' @param length number of points to draw in the confidence region
#' @param exp expansion factor of the hypercube in which random values of
#' parameters are drawn
#' @param x an object of class 'nlsConfRegions'
#' @param bounds logical defining whether bounds of the drawing hypercube are
#' plotted
#' @param ask if TRUE, draw plot interactively
#' @param ...  further arguments passed to or from other methods
#' 
#' @importFrom stats coef formula residuals qf qt runif
#' 
#' @return \code{nlsConfRegions} returns a list of four objects: \item{ cr }{ a
#' data frame containing the sample drawn in the Beale's confidence region }
#' \item{ rss }{ a vector containing the residual sums of squares corresponding
#' to \code{cr} } \item{ rss95 }{ the 95 percent residual sum of squares
#' threshold according to Beale (1960) } \item{ bounds }{ lower and upper
#' bounds of the hypercube in which random values of parameters have been drawn
#' }
#' @author Florent Baty \email{florent.baty@@gmail.com}\cr Marie-Laure
#' Delignette-Muller \email{ml.delignette@@vetagro-sup.fr}
#' @seealso \code{ellipse.nls} in the \code{ellipse} library
#' @references Beale EML (1960) Confidence regions in non-linear estimations.
#' \emph{Journal of the Royal Statistical Society}, \bold{22B}, 41-88.\cr\cr
#' Bates DM and Watts DG (1988) Nonlinear regression analysis and its
#' applications. Wiley, Chichester, UK.
#' @keywords nonlinear
#' @examples
#' 
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
#'                         (VO2rest + (VO2peak - VO2rest) * 
#'                         (1 - exp(-(t - 5.883) / mu))))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, 
#'                 mu = 1), data = O2K)
#' O2K.conf1 <- nlsConfRegions(O2K.nls1, exp = 2, length = 200)
#' plot(O2K.conf1, bounds = TRUE)
#' 
#' @export nlsConfRegions
nlsConfRegions <- function(nls, length=1000, exp=1.5){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	data	<- eval(nls$call$data, sys.frame(0))
	np		<- length(coef(nls))
	nl		<- nrow(data)
	vardep <- all.vars(formula(nls)[[2]])
	varindep <- intersect(all.vars(formula(nls)[[3]]), colnames(data))
	"formula2function"<-function(formu){
		arg1		<- all.vars(formu)
		arg2		<- vector("list",length(arg1))
		names(arg2)	<- arg1
		Args		<- do.call("alist",arg2)
		fmodele		<- as.function(c(Args,formu))
		return(fmodele)
	}

	fmodele		<- formula2function(formula(nls)[[3]])
	scer <- sum(residuals(nls)^2)
	scer95 <- scer * (1 + (np/(nl-np)) * qf(p=0.95, df1=np, df2=(nl-np), lower.tail=TRUE))
	student95	<- qt(0.975, df=nl-np)
	bornes <- cbind(coef(nls)-student95*exp*(summary(nls)$parameters)[,2], coef(nls)+student95*exp*(summary(nls)$parameters)[,2])
	colnames(bornes) <- c("Lower", "Upper")
	tirage	<- vector(length=np)
	names(tirage)	<- row.names(summary(nls)$parameters)
	tab	<- matrix(ncol=np, nrow=0)
	rss	<- vector(length=0)

	cat("  ")
	while(nrow(tab)<length){
		tirage <- apply(bornes, 1, function(z) runif(n=1, min=z[1], max=z[2]))
		listparavar	<- c(tirage, data[varindep])
		predict	<- do.call("fmodele", listparavar)
		rss1	<- sum((predict-data[,vardep])^2)
		if(rss1 < scer95){
			tenth	<- floor(100*nrow(tab)/length)
			tab	<- rbind(tab,tirage)
			rss	<- c(rss,rss1)
			if(tenth!=floor(100*nrow(tab)/length)){
				if(tenth<11){cat("\b\b",tenth,"%",sep="")}
				else{cat("\b\b\b",tenth,"%",sep="")}
			}
		}
    }
    rownames(tab) <- 1:nrow(tab)
	cat("\b\b\b100%\a")
	cat("\n Confidence regions array returned \n")
	listcr <- list(cr=tab, rss=rss, rss95=scer95, bounds=bornes)
	class(listcr) <- "nlsConfRegions"
	return(listcr)
}

#' @rdname nlsConfRegions
#' @importFrom graphics par layout plot abline
#' @export
plot.nlsConfRegions <-function(x, bounds=FALSE, ask=FALSE, ...){
	if (!inherits(x, "nlsConfRegions"))
		stop("Use only with 'nlsConfRegions' objects")
	np <- ncol(x$cr)
	def.par <- par(no.readonly = TRUE)
	if(ask) par(ask=TRUE,mar=c(4,4,3,1))
	if(!ask){
		lay <- lower.tri(matrix(0,(np-1),(np-1)), TRUE)
		lay[which(lay, TRUE)] <- 1:choose(np,2)
		layout(lay)
		par(mar=c(5,4,0.2,0.2))
	}
	for(i in 1:(np-1))
		for(j in (i+1):np){
			if(!bounds) plot(x$cr[,i], x$cr[,j], pch="+", xlab=colnames(x$cr)[i], ylab=colnames(x$cr)[j])
			else{
				xrange <- range(c(x$cr[,i], x$bounds[i,]))
				yrange <- range(c(x$cr[,j], x$bounds[j,]))
				plot(x$cr[,i], x$cr[,j], pch="+", xlab=colnames(x$cr)[i], ylab=colnames(x$cr)[j], xlim=xrange, ylim=yrange)
				abline(h=x$bounds[j,], v=x$bounds[i,], col="red", lty=2)
			}
		}

	par(def.par)	
}

#' @rdname nlsConfRegions
#' @export
print.nlsConfRegions <- function (x, ...) {
	if (!inherits(x, "nlsConfRegions"))
		stop("Use only with 'nlsConfRegions' objects")
	cat("Beale's 95 percent confidence regions\n")
	cat("\n")
	sumry <- array("", c(2, 4), list(1:2, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$rss", length(x$rss), mode(x$rss), "Residual sums of squares")
	sumry[2, ] <- c("$rss95", length(x$rss95), mode(x$rss95), "95 percent RSS threshold")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(2, 4), list(1:2, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$cr", nrow(x$cr), ncol(x$cr), "Sample drawn in the confidence region")
	sumry[2, ] <- c("$bounds", nrow(x$bounds), ncol(x$bounds), "Bounds of the drawing hypercube")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}
