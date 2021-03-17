#' NLS residuals
#' 
#' Provides several plots and tests for the analysis of residuals
#' 
#' 
#' Several plots and tests are proposed to check the validity of the
#' assumptions of the error model based on the analysis of residuals.\cr The
#' function \code{plot.nlsResiduals} proposes several plots of residuals from
#' the nonlinear fit: plot of non-transformed residuals against fitted values,
#' plot of standardized residuals against fitted values, plot of square root of
#' absolute value of standardized residuals against fitted values,
#' auto-correlation plot of residuals (i+1th residual against ith residual),
#' histogram of the non-transformed residuals and normal Q-Q plot of
#' standardized residuals.\cr \code{test.nlsResiduals} tests the normality of
#' the residuals with the Shapiro-Wilk test (shapiro.test in package stats) and
#' the randomness of residuals with the runs test (Siegel and Castellan, 1988).
#' The runs.test function used in \code{nlstools} is the one implemented in the
#' package \code{tseries}.
#' 
#' @aliases nlsResiduals plot.nlsResiduals test.nlsResiduals print.nlsResiduals
#' @param nls an object of class 'nls'
#' @param x an object of class 'nlsResiduals'
#' @param which an integer: \cr 0 = 4 graphs of residuals (types 1, 2, 4 and 6)
#' \cr 1 = non-transformed residuals against fitted values \cr 2 = standardized
#' residuals against fitted values \cr 3 = sqrt of absolute value of
#' standardized residuals against fitted values \cr 4 = auto-correlation
#' residuals (i+1th residual against ith residual) \cr 5 = histogram of the
#' residuals \cr 6 = qq-plot of the residuals
#' @param ...  further arguments passed to or from other methods
#' 
#' @importFrom stats fitted residuals qt resid coef
#' 
#' @return \code{nlsResiduals} returns a list of five objects: \item{ std95 }{
#' the Student value for alpha=0.05 (bilateral) and the degree of freedom of
#' the model } \item{ resi1 }{ a matrix with fitted values vs. non-transformed
#' residuals } \item{ resi2 }{ a matrix with fitted values vs. standardized
#' residuals } \item{ resi3 }{ a matrix with fitted values vs.
#' sqrt(abs(standardized residuals)) } \item{ resi4 }{ a matrix with ith
#' residuals vs. i+1th residuals }
#' @author Florent Baty \email{florent.baty@@gmail.com}\cr Marie-Laure
#' Delignette-Muller \email{ml.delignette@@vetagro-sup.fr}
#' @references Bates DM and Watts DG (1988) Nonlinear regression analysis and
#' its applications. Wiley, Chichester, UK.\cr\cr Siegel S and Castellan NJ
#' (1988) Non parametric statistics for behavioral sciences. McGraw-Hill
#' international, New York.
#' 
#' 
#' @keywords nonlinear
#' @examples
#' 
#' # Plots of residuals
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2rest + (t > 5.883) * 
#'                         (VO2rest + (VO2peak - VO2rest) * 
#'                         (1 - exp(-(t - 5.883) / mu))))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2rest = 400, VO2peak = 1600, mu = 1), 
#'                data = O2K)
#' O2K.res1 <- nlsResiduals(O2K.nls1)
#' plot(O2K.res1, which = 0)
#' 
#' # Histogram and qq-plot
#' plot(O2K.res1, which = 5)
#' plot(O2K.res1, which = 6)
#' 	
#' # Tests
#' test.nlsResiduals(O2K.res1)
#' 
#' @export nlsResiduals
nlsResiduals <- function(nls){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	#which = 1 residuals
	resi1		<- cbind(fitted(nls), residuals(nls))
	colnames(resi1)	<- c("Fitted values", "Residuals")

	#which = 2 residuals
	ecartres	<- summary(nls)$sigma
	nresiduals 	<- (residuals(nls)-mean(residuals(nls)))/ecartres
	std95		<- qt(0.975,df=length(resid(nls))-length(coef(nls)))
	resi2		<- cbind(fitted(nls), nresiduals)
	colnames(resi2)	<- c("Fitted values", "Standardized residuals")

	#which = 3 residuals
	ecartres	<- summary(nls)$sigma
	nresiduals	<- residuals(nls)/ecartres
	resi3		<- cbind(fitted(nls), sqrt(abs(nresiduals)))
	colnames(resi3)	<- c("Fitted values", "Sqrt abs. standardized residuals")

	#which = 4 residuals
	resiminus	<-vector()
	resiplus	<-vector()
       	for(i in 1:(length(residuals(nls))-1)){
               	resiminus[i]	<- residuals(nls)[i]
               	resiplus[i]	<- residuals(nls)[i+1]
        }
	resi4	<- cbind(resiminus, resiplus)
	colnames(resi4)	<- c("Residuals i", "Residuals i+1")

	listresi	<- list(resi1=resi1, resi2=resi2, resi3=resi3, resi4=resi4, std95=std95)	
	class(listresi)	<- "nlsResiduals"
	return(listresi)
}

#' @rdname nlsResiduals
#' @importFrom graphics hist boxplot par plot abline hist
#' @importFrom stats ppoints quantile qnorm qqnorm qqline
#' @export
plot.nlsResiduals <- function(x, which=0, ...){

	hist.nlsResiduals <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	hist(x$resi1[,"Residuals"], main="Residuals", xlab="Residuals")
	}

	boxplot.nlsResiduals <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	boxplot(x$resi1[,"Residuals"], main="Residuals")
	}

	qq.nlsResiduals <- function(x){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")

	qqploty	<- sort(x$resi2[,2])
	qqplotx	<- qnorm(ppoints(nrow(x$resi2)))
	qy	<- quantile(qqploty,c(0.25,0.75))
	qx	<- qnorm(c(0.25,0.75))
	slope	<- diff(qy)/diff(qx)
	ori	<- qy[1]-slope*qx[1]     
	
	qqnorm(x$resi2[,2], main="Normal Q-Q Plot of\n Standardized Residuals")
	qqline(x$resi2[,2])
	}

	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	if(!(which %in% 0:6)) 
		stop("\n Expected 'which':\n 1 = un-transformed\n 2 = normed\n 3 = sqrt absolute normed\n 4 = auto-correlation\n 5 = histogram\n 6 = qq-plot")
	if(which == 0){
		def.par <- par(no.readonly = TRUE)
        par(mfrow=c(2,2))
        plot.nlsResiduals(x, which=1)
        plot.nlsResiduals(x, which=2)
        plot.nlsResiduals(x, which=4)
        qq.nlsResiduals(x)
        par(def.par)
	}
	if(which == 1){
		plot(x$resi1, xlab="Fitted values", ylab="Residuals", main="Residuals")
		abline(h=0, lty=2)
	}
	if(which == 2){
		yrange	<- c(min(min(x$resi2[,2]), -x$std95), max(max(x$resi2[,2]), x$std95))
		plot(x$resi2, xlab="Fitted values", ylab="Standardized residuals", ylim=yrange, main="Standardized Residuals")
		abline(h=0,lty=2); abline(h=x$std95); abline(h=-x$std95)
	}
	if(which == 3){
		plot(x$resi3, xlab="Fitted values", ylab=expression(sqrt(abs("Standardized residuals"))), main="Sqrt abs residuals")
	}
	if(which == 4){
		plot(x$resi4, xlab="Residuals i", ylab="Residuals i+1", main="Autocorrelation")
		abline(h=0, lty=2)
	}
	if(which == 5){
		hist(x)
	}
	if(which == 6){
		qq.nlsResiduals(x)
	}
}

#' @rdname nlsResiduals
#' @importFrom stats pnorm shapiro.test
#' @export
test.nlsResiduals <- function(x){
	"runs.test" <- function (x, alternative = c("two.sided", "less", "greater")){
		if(!is.factor(x))
			stop("x is not a factor")
		if(any(is.na(x)))
			stop("NAs in x")
		if(length(levels(x)) != 2)
			stop("x does not contain dichotomous data")
		alternative <- match.arg(alternative)
		DNAME <- deparse(substitute(x))
		n <- length(x)
		R <- 1 + sum(as.numeric(x[-1] != x[-n]))
		n1 <- sum(levels(x)[1] == x)
		n2 <- sum(levels(x)[2] == x)
		m <- 1 + 2*n1*n2 / (n1+n2)
		s <- sqrt(2*n1*n2 * (2*n1*n2 - n1 - n2) / ((n1+n2)^2 * (n1+n2-1)))
		STATISTIC <- (R - m) / s
		METHOD <- "Runs Test"
		if(alternative == "two.sided")
			PVAL <- 2 * pnorm(-abs(STATISTIC))
		else if(alternative == "less")
			PVAL <- pnorm(STATISTIC)
		else if(alternative == "greater")
			PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
		else stop("irregular alternative")
		names(STATISTIC) <- "Standard Normal"
		structure(list(statistic = STATISTIC,
			alternative = alternative,
			p.value = PVAL,
			method = METHOD,
			data.name = DNAME),
			class = "htest")
	}

	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	#Shapiro-Wilk test
	stdres <- x$resi2[,2]
	shapi	<- shapiro.test(stdres)
	cat("\n------")
	print(shapi)

	#Test runs
	run	<- vector(length=nrow(x$resi1))
	for(i in 1:nrow(x$resi1)){
		if(x$resi2[i,2]<0) run[i] <- "N" 
        if(x$resi2[i,2]>0) run[i] <- "P"
	}
	cat("\n------")
	runs.test(as.factor(run))
}

#' @rdname nlsResiduals
#' @export
print.nlsResiduals <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	cat("Residuals\n")
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$std95", length(x$std95), mode(x$std95), "Student value for alpha = 0.05")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(4, 4), list(1:4, c("matrix", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$resi1", nrow(x$resi1), ncol(x$resi1), "fitted values vs. non-transf. resid")
	sumry[2, ] <- c("$resi2", nrow(x$resi2), ncol(x$resi2), "fitted values vs. standardized resid")
	sumry[3, ] <- c("$resi3", nrow(x$resi3), ncol(x$resi3), "fitted values vs. sqrt abs std resid")
	sumry[4, ] <- c("$resi4", nrow(x$resi4), ncol(x$resi4), "resid i vs. resid i+1")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}
