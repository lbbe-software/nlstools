#' Approximation of hat values for nls.
#'
#' @param model An object that inherits from class nls
#' @param ... further arguments for compatibility with S3 method.
#'
#' @return Vector of approximated hat values
#'
#' @details https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
#' @importFrom stats hatvalues
#' @export
#' @examples 
#' DNase1 <- subset(DNase, Run == 1)
#' ## using a selfStart model
#' fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#' hatvalues(fm1DNase1)
hatvalues.nls <- function(model, ...) {
  stopifnot(methods::is(model, 'nls'))
  V <- model$m$gradient()
  Q1 <- qr.Q(qr(V))
  rowSums(Q1*Q1)
}

#' Residuals of all types for nls models
#'
#' @param object An object that inherits from class nls
#' @param type the type of residuals which should be returned.
#' Can be abbreviated.
#' @param ... further arguments for compatibility with S3 method.
#'
#' @details
#' http://statweb.stanford.edu/~jtaylo/courses/stats306b/restricted/notebooks/quasilikelihood.pdf
#' @export
#' @examples 
#' DNase1 <- subset(DNase, Run == 1)
#' ## using a selfStart model
#' fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#' residuals_nls(fm1DNase1, type = "response")
#' residuals_nls(fm1DNase1, type = "pearson")
#' residuals_nls(fm1DNase1, type = "partial")
residuals_nls <- function(
  object,
  type = c("working", "response", "deviance",
           "pearson", "partial"), ...)
{
  type <- match.arg(type)
  r <- summary(object)$residuals
  res <- switch(type, working = , response = r, deviance = ,
                pearson = if (is.null(object$weights)) r else r * sqrt(object$weights),
                partial = r)
  res <- stats::naresid(object$na.action, res)
  if (type == "partial")
    res <- res + predict(object, type = "terms")
  res
}

#' weighted.residuals for nls
#'
#' @param obj An object that inherits from class nls
#' @inheritParams stats::weighted.residuals
#' @export
#' @examples 
#' DNase1 <- subset(DNase, Run == 1)
#' ## using a selfStart model
#' fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#' weighted.residuals.nls(fm1DNase1)
weighted.residuals.nls <- function(obj, drop0 = TRUE)
{
  w <- stats::weights(obj)
  # Les résidus de deviance sont les résidus classiques
  r <- residuals_nls(obj)
  if (drop0 && !is.null(w)) {
    if (is.matrix(r))
      r[w != 0, , drop = FALSE]
    else r[w != 0]
  }
  else r
}

#' Cooks.distance for nls
#'
#' @param model An object that inherits from class nls
#' @inheritParams stats::cooks.distance
#' @importFrom stats cooks.distance
#' @export
#' @examples 
#' DNase1 <- subset(DNase, Run == 1)
#' ## using a selfStart model
#' fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#' cooks.distance(fm1DNase1)
cooks.distance.nls <- function(
  model,
  # infl = lm.influence(model, do.coef = FALSE),
  res = weighted.residuals.nls(model),
  sd = sqrt(stats::deviance(model)/df.residual(model)),
  hat = hatvalues.nls(model), ...)
{
  # Rank is number of fitted parameters
  p <- length(coef(model)) # model$rank
  res <- ((res/(sd * (1 - hat)))^2 * hat)/p
  res[is.infinite(res)] <- NaN
  stats::setNames(res, seq_along(res))
}
