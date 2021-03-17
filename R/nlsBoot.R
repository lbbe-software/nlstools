#' Bootstrap resampling
#' 
#' Bootstrap resampling
#' 
#' Non-parametric bootstrapping is used. Mean centered residuals are
#' bootstrapped. By default, 999 resampled data sets are created from which
#' parameter estimates are obtained by fitting the model on each of these data
#' sets. Whenever the fit fails to converge, a flag reports the number of
#' non-convergences. If the fitting procedure fails to converge in more than
#' 50\% of the cases, the procedure is interrupted with a flag and no result is
#' given. The function \code{summary} returns the bootstrap estimates (mean and
#' std. dev. of the bootstrapped estimates) and the median and 95 percent
#' confidence intervals (50, 2.5, and 97.5 percentiles of the bootstrapped
#' estimates). The bootstrapped estimate distributions can be visualized using
#' the function \code{plot.nlsBoot} either by plotting the bootstrapped sample
#' for each pair of parameters or by displaying the boxplot representation of
#' the bootstrapped sample for each parameter. Notice that \code{nlsBoot} does
#' not currently handle transformed dependent variables specified in the left
#' side of the \code{nls} formula.
#' 
#' @aliases nlsBoot plot.nlsBoot print.nlsBoot summary.nlsBoot
#' @param nls an object of class 'nls'
#' @param niter number of iterations
#' @param x,object an object of class 'nlsBoot'
#' @param type type of representation (options are "pairs" or "boxplot")
#' @param mfr layout definition (number of rows and columns in the graphics
#' device)
#' @param ask if TRUE, draw plot interactively
#' @param ...  further arguments passed to or from other methods
#' 
#' @importFrom stats fitted resid formula update coef quantile sd
#' 
#' @return \code{nlsBoot} returns a list of 4 objects: \item{ coefboot }{
#' contains the bootstrap parameter estimates } \item{ bootCI }{ contains the
#' bootstrap medians and the bootstrap 95\% confidence intervals } \item{
#' estiboot }{ contains the means and std. errors of the bootstrap parameter
#' estimates } \item{ rse }{ is the vector of bootstrap residual errors }
#' @author Florent Baty \email{florent.baty@@gmail.com}\cr Marie-Laure
#' Delignette-Muller \email{ml.delignette@@vetagro-sup.fr}
#' @references Bates DM and Watts DG (1988) Nonlinear regression analysis and
#' its applications. Wiley, Chichester, UK.\cr\cr Huet S, Bouvier A, Poursat
#' M-A, Jolivet E (2003) Statistical tools for nonlinear regression: a
#' practical guide with S-PLUS and R examples. Springer, Berlin, Heidelberg,
#' New York.
#' @keywords nonlinear
#' @examples
#' 
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2res + (t > 5.883) * 
#'                         (VO2res + (VO2peak - VO2res) * 
#'                         (1 - exp(-(t - 5.883) / mu))))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2res = 400, VO2peak = 1600, 
#'                 mu = 1), data = O2K)
#' O2K.boot1 <- nlsBoot(O2K.nls1, niter = 200)
#' plot(O2K.boot1)
#' plot(O2K.boot1, type = "boxplot", ask = FALSE)
#' summary(O2K.boot1)
#'   
#' @export nlsBoot
nlsBoot <-function(nls, niter=999){

	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	data2 <- eval(nls$data, sys.frame(0))
	fitted1 <- fitted(nls)
	resid1 <- resid(nls)
	var1 <- all.vars(formula(nls)[[2]])
	
	l1 <- lapply(1:niter, function(i){
		data2[,var1] <- fitted1 + sample(scale(resid1, scale=FALSE), replace=TRUE);
		nls2 <- try(update(nls, start=as.list(coef(nls)), data=data2), silent=TRUE);
		if(inherits(nls2, "nls"))
			return(list(coef=coef(nls2), rse=summary(nls2)$sigma))
		})

	if(sum(sapply(l1, is.null)) > niter/2) stop(paste("Procedure aborted: the fit only converged in", round(sum(sapply(l1, is.null))/niter), "% during bootstrapping"))

	tabboot <- sapply(l1[!sapply(l1, is.null)], function(z) z$coef)
	rseboot <- sapply(l1[!sapply(l1, is.null)], function(z) z$rse)
	recapboot <- t(apply(tabboot, 1, quantile, c(.5, .025, .975))); colnames(recapboot) <- c("Median", "2.5%", "97.5%")
	estiboot <- t(apply(tabboot, 1, function(z) c(mean(z), sd(z)))); colnames(estiboot) <- c("Estimate", "Std. error")
	
	serr <- sum(sapply(l1, is.null))
	if(serr > 0) warning(paste("The fit did not converge", serr, "times during bootstrapping"))
	
	listboot <- list(coefboot = t(tabboot), rse = rseboot, bootCI = recapboot, estiboot = estiboot)
	class(listboot) <- "nlsBoot"
	return(listboot)
	
}

#' @rdname nlsBoot
#' @importFrom graphics par layout plot boxplot
#' @export
plot.nlsBoot <-function(x, type=c("pairs","boxplot"), mfr=c(ceiling(sqrt(ncol(x$coefboot))),ceiling(sqrt(ncol(x$coefboot)))),ask=FALSE, ...){
	if (!inherits(x, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	tab <- x$coefboot
	np <- ncol(tab)
 	def.par <- par(no.readonly = TRUE)	
	if(type[1] == "pairs"){
		if(ask) par(ask=TRUE, mar=c(4,4,3,1))
		if(!ask){
			lay <- lower.tri(matrix(0,(np-1),(np-1)), TRUE)
			lay[which(lay, TRUE)] <- 1:choose(np,2)
			layout(lay)
			par(mar=c(5,4,0.2,0.2))
		}
		for(i in 1:(np-1))
			for(j in (i+1):np)
				plot(tab[,i], tab[,j], xlab=colnames(tab)[i], ylab=colnames(tab)[j], pch="+")
	}
	if(type[1] == "boxplot"){ 
		if(ask) par(ask=TRUE, mar=c(4,4,3,1))
		if(!ask) par(mfrow=mfr, mar=c(4,4,3,1))
		for(i in 1:np){
			boxplot(tab[,i],main=colnames(tab)[i])
		}
	}
	par(def.par)
}

#' @rdname nlsBoot
#' @export
print.nlsBoot <- function (x, ...) {
	if (!inherits(x, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	cat("Bootstrap resampling\n")
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$rse", length(x$rse), mode(x$rse), "Bootstrap residual errors")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(3, 4), list(1:3, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$coefboot", nrow(x$coefboot), ncol(x$coefboot), "Bootstrap parameter estimates")
	sumry[2, ] <- c("$estiboot", nrow(x$estiboot), ncol(x$estiboot), "Bootstrap estimates and std. error")
	sumry[3, ] <- c("$bootCI", nrow(x$bootCI), ncol(x$bootCI), "Bootstrap medians and 95% CI")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

#' @rdname nlsBoot
#' @export
summary.nlsBoot <- function (object, ...) {
	if (!inherits(object, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	cat("\n------\n")
	cat("Bootstrap statistics\n")
	print(object$estiboot)
	cat("\n------\n")
	cat("Median of bootstrap estimates and percentile confidence intervals\n")
	print(object$bootCI)
	cat("\n")
}

