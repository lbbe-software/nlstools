#' Nonlinear least squares fit
#' 
#' Tools to help the fit of nonlinear models with nls
#' 
#' The function \code{preview} helps defining the parameter starting values
#' prior fitting the model. It provides a superimposed plot of observed
#' (circles) and predicted (crosses) values of the dependent variable versus
#' one of the independent variables with the model evaluated at the starting
#' values of the parameters. The function \code{overview} returns the
#' parameters estimates, their standard errors as well as their asymptotic
#' confidence intervals and the correlation matrix (alternately, the function
#' \code{confint} provides better confidence interval estimates whenever it
#' converges). \code{plotfit} displays a superimposed plot of the dependent
#' variable versus one the independent variables together with the fitted
#' model.
#' 
#' @param formula formula of a non-linear model
#' @param data a data frame with header matching the variables given in the
#' formula
#' @param start a list of parameter starting values which names match the
#' parameters given in the formula
#' @param variable index of the variable to be plotted against the predicted
#' values; default is the first independent variable as it appears in the
#' orginal dataset
#' @param x an object of class 'nls'
#' @param smooth a logical value, default is FALSE. If smooth is TRUE, a plot
#' of observed values is plotted as a function of 1000 values continuously
#' taken in the range interval [min(variable),max(variable)]. This option can
#' only be used if the number of controlled variables is 1.
#' @param xlab X-label
#' @param ylab Y-label
#' @param pch.obs type of point of the observed values
#' @param pch.fit type of point of the fitted values (not applicable if
#' smooth=TRUE)
#' @param lty type of line of the smoothed fitted values (if smooth=TRUE)
#' @param lwd thickness of line of the smoothed fitted values (if smooth=TRUE)
#' @param col.obs color of the observed points
#' @param col.fit color of the fitted values
#' @param ...  further arguments passed to or from other methods
#' 
#' 
#' @importFrom stats residuals coef qt
#' @importFrom graphics plot points
#' 
#' @seealso \code{nls} in the \code{stats} library and \code{confint.nls} in
#' the package \code{MASS}
#' 
#' @references Baty F, Ritz C, Charles S, Brutsche M, Flandrois J-P,
#' Delignette-Muller M-L (2015). A Toolbox for Nonlinear Regression in R: The
#' Package nlstools. \emph{Journal of Statistical Software}, \bold{66}(5),
#' 1-21.\cr\cr Bates DM and Watts DG (1988) Nonlinear regression analysis and
#' its applications. Wiley, Chichester, UK.
#' @keywords nonlinear
#' @examples
#' 
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
#'                         (VO2rest + (VO2peak - VO2rest) * 
#'                         (1 - exp(-(t - 5.883) / mu))))
#' preview(formulaExp, O2K, list(VO2rest = 400, VO2peak = 1600, mu = 1))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, 
#'                 mu = 1), data = O2K)
#' overview(O2K.nls1)
#' plotfit(O2K.nls1, smooth = TRUE)
#' 
#' @export

preview <- function(formula, data, start, variable=1){

	formula2function <- function(formu){
		arg1		<- all.vars(formu)
		arg2		<- vector("list",length(arg1))
		names(arg2)	<- arg1
		Args		<- do.call("alist",arg2)
		fmodele		<- as.function(c(Args,formu))
		return(fmodele)
	}

	f1 <- formula2function(formula(formula)[[3]])
	vardep <- all.vars(formula[[2]])
	varindep <- intersect(all.vars(formula[[3]]), colnames(data))
	predic <- do.call(f1, as.list(c(start, data[varindep])))
	rss1 <- signif(sum((predic-data[vardep])^2), 3)

	plot(data[c(variable, which(colnames(data)==vardep))], ylab="Predicted", main="", ylim=c(min(data[vardep],predic), max(data[vardep],predic)))
	points(cbind.data.frame(data[variable], predic), pch="+", col="red")
	cat("\nRSS: ", rss1,"\n")

}

#' @rdname preview
#' @importFrom graphics plot lines points
#' @importFrom stats formula predict
#' @export
plotfit <- function(x, smooth=FALSE, variable=1, xlab=NULL, ylab=NULL, pch.obs=1, pch.fit="+", lty=1, lwd=1, col.obs="black", col.fit="red", ...){
	if (!inherits(x, "nls"))
		stop("Use only with 'nls' objects")
	d <- eval(x$call$data, sys.frame(0))
	vardep <- all.vars(formula(x)[[2]])
	varindep <- intersect(all.vars(formula(x)[[3]]), colnames(d))
	variable1 <- which(varindep == colnames(d)[variable])
	if (smooth & length(varindep)!=1) 
        	stop("smooth option is only possible when the number of independent variables equals 1")
	if(smooth | smooth=="T"){
        w0 <- list(seq(min(d[,varindep]), max(d[,varindep]), len=1000))
		names(w0) <- varindep
		if(is.null(xlab)) xlab <- varindep
		if(is.null(ylab)) ylab <- vardep
		plot(d[c(varindep, vardep)], xlab=xlab, ylab=ylab, pch=pch.obs, col=col.obs, ...)
		lines(w0[[1]], predict(x,new=w0), col=col.fit, lty=lty, lwd=lwd)
	}
	else{
		if(is.null(xlab)) xlab <- varindep[variable1]
		if(is.null(ylab)) ylab <- vardep
		plot(d[,vardep] ~ d[,varindep[variable1]], xlab=xlab, ylab=ylab, pch=pch.obs, col=col.obs, ...)
		points(d[,varindep[variable1]], predict(x), pch=pch.fit, col=col.fit)
	}
}

#' @rdname preview
#' @export
overview <- function(x){
	if (!inherits(x, "nls"))
		stop("Use only with 'nls' objects")
	cat("\n------")
	print(summary(x))
	cat("------\n")
	cat("Residual sum of squares:", signif(sum(residuals(x)^2), 3),"\n\n")
	n <- length(residuals(x))
	np <- length(coef(x))
	esti <- summary(x)$parameters[,"Estimate"]
	ster <- summary(x)$parameters[,"Std. Error"]
	t95 <- qt(0.975, df=(n-np))
	binf <- esti - t95 * ster
	bsup <- esti + t95 * ster
	cat("------\n")
	cat("t-based confidence interval:\n")
	print(cbind.data.frame("2.5%" = binf, "97.5%" = bsup))
	cat("\n")
	cat("------\n")
	cat("Correlation matrix:\n")
	print(summary(x, correlation = TRUE)$correlation)
	cat("\n")
}
