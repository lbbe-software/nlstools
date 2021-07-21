"nlsBootPredict" <- function (nlsBoot, newdata, 
                              interval = c("confidence", "prediction")) 
{
  if (!inherits(nlsBoot, "nlsBoot"))
    stop("Use only with 'nlsBoot' objects")
  if (missing(newdata))
  {
    newdata <- get(as.character(nlsBoot$nls$data))
  }
  nlsformula <- formula(nlsBoot$nls)
  nlsresid <- resid(nlsBoot$nls)
  param <- nlsBoot$coefboot
  bootparam <- nlsBoot$coefboot
  niter <- length(nlsBoot$rse)
  
  "formula2function"<-function(formu){
    arg1		<- all.vars(formu)
    arg2		<- vector("list",length(arg1))
    names(arg2)	<- arg1
    Args		<- do.call("alist",arg2)
    fmodele		<- as.function(c(Args,formu))
    return(fmodele)
  }
  f1 <- formula2function(formula(nlsformula)[[3]])
  vardep <- all.vars(nlsformula[[2]])
  varindep <- intersect(all.vars(nlsformula[[3]]), colnames(newdata))

  ## vector of mean predictions on newdata with one bootstrap sample  
  one.mean.pred <- function(i)
  {
    do.call(f1, as.list(c(param[i,], newdata[varindep])))
  }
  boot.mean.pred <- sapply(1:niter, one.mean.pred)
  
  if (interval == "confidence")
  {
    recap.boot.mean.pred <- t(apply(boot.mean.pred, 1, 
                                    quantile, c(.5, .025, .975))) 
    colnames(recap.boot.mean.pred) <- c("Median", "2.5%", "97.5%")
    return(recap.boot.mean.pred)
  } else
  {
    ## vector of individual predictions on newdata with one bootstrap sample
    one.indiv.pred <- function(i)
    {
      boot.mean.pred[, i] + sample(scale(nlsresid, scale=FALSE), 
                                   size = nrow(newdata), replace=TRUE)
    }
    boot.indiv.pred <- sapply(1:niter, one.indiv.pred)
    recap.boot.indiv.pred <- t(apply(boot.indiv.pred, 1, 
                                     quantile, c(.5, .025, .975))) 
    colnames(recap.boot.indiv.pred) <- c("Median", "2.5%", "97.5%")
    return(recap.boot.indiv.pred)
  }
}



