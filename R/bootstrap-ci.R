#' Function for bootstraping nls model and retrieve std errors
#'
#' @param dataset original dataset
#' @param i vector of indices to bootstrap in dataset
#' @param model nls (or lm) model on which to do bootstrap
#' @param newdata new dataset for which to provide std errors
#'
#' @details http://r.789695.n4.nabble.com/standard-errors-for-predict-nls-td873670.html
pred.bs <- function(dataset, i, model, newdata) {
  # Same data but randomised residuals
  # res <- scale(residuals(model), scale = FALSE) 
  # dataset$degradation = predict(model) + res[i] 
  ## intercept errors
  
  n1 <- try(stats::nls(formula(model), dataset[i,],
                       start = coef(model)), silent = TRUE)
  
  if (class(n1) == "try-error") {
    r <- rep(NA, length(newdata))
    # r <- rep(NA, length(coef(model)))
  } else {
    r <- stats::predict(n1, newdata = newdata)
    # r <- coef(n1)
  }
  r
}

#' Extract bootstrap confidence intervals
#' 
#' @inheritParams boot::boot.ci

getbootci <- function(boot.out, index, conf = 0.95) {
  a <- (1 - conf)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  
  b <- boot::boot.ci(boot.out, type = "norm",
               conf = conf, index = index)$normal
  b <- data.frame(low = b[2], high = b[3])
  names(b) <- pct
  b
}

#' Confidence interval for predictions with nls
#'
#' @param object nls model on which to apply
#' @param newdata vector of new observation
#' @param R  number of bootstrap to realise. Higher 
#' gives more precision but longer calculations
#' @param conf confidence interval to return
#' @export
#' @examples 
#' formulaExp <- as.formula(VO2 ~ (t <= 5.883) * VO2res + (t > 5.883) *
#'                            (VO2res + (VO2peak - VO2res) *
#'                               (1 - exp(-(t - 5.883) / mu))))
#' O2K.nls1 <- nls(formulaExp, start = list(VO2res = 400, VO2peak = 1600,
#'                                          mu = 1), data = O2K)
#' predict_boot_nls(O2K.nls1)
predict_boot_nls <- function(object, newdata, R = 1000, conf = 0.95) {
  if (!inherits(object, "nls"))
    stop("Use only with 'nls' objects")
  
  dataset <- eval(object$data, sys.frame(0))
  if (missing(newdata)) {newdata <- dataset}
  
  index <- 1:nrow(newdata)
  
  nls.boot <- boot::boot(dataset, pred.bs, R = R,
                         model = object,
                         newdata = newdata)
  # names(nls.boot)
  # object$m$predict()
  # object$m$predict(newdata)
  # object$m$form()
  # stats::predict(n1, newdata = newdata)
  # Comparaison with lm
  # predict(mod_lm, newdata = list(duree = newdata), se.fit = TRUE)
  
  rn <- paste0("t", index)
  op_list <- lapply(index, function(x) {
    op <- boot::imp.moments(nls.boot, index = x)$rat
    data.frame(op1 = op[1], op2 = op[2])
  })
  op <- do.call("rbind", op_list)
  
  op <- cbind(nls.boot$t0, op[, 1L] - nls.boot$t0, 
              (op[, 1L] - nls.boot$t0)/op[, 1L],
              sqrt(op[, 2L]),
              apply(nls.boot$t, 
                    2L, mean, na.rm = TRUE))
  dimnames(op) <- list(rn, c("fit", "boot.bias",
                             "boot.relative.bias",
                             "se.fit", "boot.mean"))

  bootvals_list <- lapply(
    index, 
    function(x, boot.out, conf) getbootci(boot.out, x, conf),
    boot.out = nls.boot, conf = conf
  )
  bootvals <- do.call("rbind", bootvals_list)
  
  cbind(op, bootvals)
}


