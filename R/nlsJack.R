#' Jackknife resampling
#' 
#' Jackknife resampling
#' 
#' 
#' A jackknife resampling procedure is performed. Each observation is
#' sequentially removed from the initial data set using a leave-one-out
#' strategy. A data set with \emph{n} observations provides thus \emph{n}
#' resampled data sets of \emph{n-1} observations. The jackknife estimates with
#' confidence intervals are calculated as described by Seber and Wild (1989)
#' from the results of \emph{n} new fits of the model on the \emph{n} jackknife
#' resampled data sets. The leave-one-out procedure is also employed to assess
#' the influence of each observation on each parameter estimate. An observation
#' is empirically defined as influential for one parameter if the difference
#' between the estimate of this parameter with and without the observation
#' exceeds twice the standard error of the estimate divided by \emph{sqrt(n)}.
#' This empirical method assumes a small curvature of the nonlinear model. For
#' each parameter, the absolute relative difference (in percent of the
#' estimate) of the estimates with and without each observation is plotted. An
#' asterisk is plotted for each influential observation.
#' 
#' @aliases nlsJack plot.nlsJack print.nlsJack summary.nlsJack
#' @param nls an object of class 'nls'
#' @param x,object an object of class 'nlsJack'
#' @param mfr layout definition, default is k rows (k: number of parameters)
#' and 1 column
#' @param ask if TRUE, draw plot interactively
#' @param ...  further arguments passed to or from other methods
#' 
#' @importFrom stats coef update residuals qt
#' 
#' @return \code{nlsJack} returns a list with 7 objects: \item{ estijack }{ a
#' data frame with jackknife estimates and bias } \item{ coefjack }{ a data
#' frame with the parameter estimates for each jackknife sample } \item{ reldif
#' }{ a data frame with the absolute relative difference (in percent of the
#' estimate) of the estimates with and without each observation } \item{ dfb }{
#' a data frame with dfbetas for each parameter and each observation } \item{
#' jackCI }{ a data frame with jackknife confidence intervals } \item{ rse }{ a
#' vector with residual standard error for each jackknife sample } \item{ rss
#' }{ residual a vector with residual sum of squares for each jackknife sample
#' }
#' @author Florent Baty \email{florent.baty@@gmail.com}\cr Marie-Laure
#' Delignette-Muller \email{ml.delignette@@vetagro-sup.fr}
#' @references Seber GAF, Wild CJ (1989) Nonlinear regression. Wiley, New
#' York.\cr\cr
#' @keywords nonlinear
#' @examples
#' 
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
#'                         (VO2rest + (VO2peak - VO2rest) * 
#'                         (1 - exp(-(t - 5.883) / mu))))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, mu = 1), 
#'                data = O2K)
#' O2K.jack1 <- nlsJack(O2K.nls1)
#' plot(O2K.jack1)
#' summary(O2K.jack1)
#' 
#' @export nlsJack
nlsJack <-function(nls){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	c1 <- nls$call
	data <- eval(c1$data, sys.frame(0))
	c1$start <- as.list(coef(nls))
	nl <- nrow(data)
	np <- length(coef(nls))
	
	l1 <- lapply(1:nl, function(z){
		nls2 <- update(nls, data=data[-z,], start=as.list(coef(nls)));
		return(list(coef=coef(nls2), sigma=summary(nls2)$sigma, rss=sum(residuals(nls2)^2), dfb=abs(coef(nls2)-coef(nls))/(summary(nls)$parameters[,2])))
	})

	tabjack <- t(sapply(l1, function(z) z$coef))
	rsejack <- sapply(l1, function(z) z$sigma)
	rssjack <- sapply(l1, function(z) z$rss)
	dfb <- t(sapply(l1, function(z) z$dfb))
	reldif	<- apply(tabjack, 1, function(x) 100*abs(x-coef(nls))/coef(nls))
	pseudo <- t(apply((nl-1) * tabjack, 1, function(z) nl * coef(nls) - z))
	estijack <- cbind.data.frame(Estimates=colSums(pseudo) / nl, Bias=coef(nls) - colSums(pseudo) / nl)
	sum1 <- crossprod(t(t(pseudo) - estijack$Estimates))
	varjack <- (1 / (nl * (nl - 1))) * sum1
	student95 <- qt(0.975, df = nl - np)
	ICjack <- cbind.data.frame(Esti = estijack$Estimates, Low = estijack$Estimates - student95 * sqrt(diag(varjack)), Up = estijack$Estimates + student95 * sqrt(diag(varjack)))

	listjack	<-list(estijack=estijack, coefjack=tabjack, reldif=reldif, rse=rsejack, rss=rssjack ,dfb=dfb, jackCI=ICjack)
	class(listjack) <- "nlsJack"
	return(listjack)
}

#' @rdname nlsJack
#' @importFrom graphics par plot text
#' @export
plot.nlsJack <- function(x, mfr=c(nrow(x$reldif),1), ask=FALSE, ...){
	if (!inherits(x, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	if(ask) par(ask=TRUE,mar=c(4,4,3,1))
	if(!ask) par(mfrow=mfr,mar=c(4,4,3,1))
	for(i in 1:nrow(x$reldif)){
		plot(x$reldif[i,],type="h",main=rownames(x$reldif)[i],xlab="Observation #",ylab="Rel Diff (%)",ylim=c(0,1.2*max(x$reldif[i,])))
		for(j in 1:nrow(x$dfb)){
			if(x$dfb[j,i]>(2/sqrt(nrow(x$dfb))))
				text(j,x$reldif[i,j],"*",col="red",cex=2)
		}
	}
	par(mfrow=c(1,1),ask=FALSE)
}

#' @rdname nlsJack
#' @export
print.nlsJack <- function (x, ...) {
	if (!inherits(x, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	cat("Jackknife resampling\n")
	cat("\n")
	sumry <- array("", c(3, 4), list(1:3, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$estijack", length(x$estijack), mode(x$estijack), "jackknife estimates and bias")
	sumry[2, ] <- c("$rse", length(x$rse), mode(x$rse), "residual errors")
	sumry[3, ] <- c("$rss", length(x$rss), mode(x$rss), "residual sum of squares")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$coefjack", nrow(x$coefjack), ncol(x$coefjack), "parameter estimates")
	sumry[2, ] <- c("$reldif", nrow(x$reldif), ncol(x$reldif), "relative diff. of parameter estimates")
	sumry[3, ] <- c("$dfb", nrow(x$dfb), ncol(x$dfb), "DF-beta")
	sumry[4, ] <- c("$jackCI", nrow(x$jackCI), ncol(x$jackCI), "jackknife confidence intervals")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

#' @rdname nlsJack
#' @export
summary.nlsJack <- function (object, ...) {
	if (!inherits(object, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	cat("\n------\n")
	cat("Jackknife statistics\n")
	print(object$estijack)
	cat("\n------\n")
	cat("Jackknife confidence intervals\n")
	print(object$jackCI[,c("Low","Up")])
	cat("\n------\n")
	cat("Influential values\n")
	inf1 <- which(object$dfb>(2/sqrt(nrow(object$dfb))), arr.ind=TRUE)
	for(i in 1:nrow(inf1))
		cat("* Observation", inf1[i,1], "is influential on", rownames(object$estijack)[inf1[i,2]], "\n")
	cat("\n")
}
