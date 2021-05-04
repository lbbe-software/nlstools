# Test nlsBootpredict on a linear model
library(nlstools)
set.seed(123)
n <- 30
x <- rnorm(n)
eps <- 0.3 * rnorm(n)
d <- data.frame(x = x, y = x + eps)
model <- lm(y ~ x, data = d)
new <- data.frame(x = seq(-3, 3, 0.1))
pred.plim <- predict(model, new, interval = "prediction")
pred.clim <- predict(model, new, interval = "confidence")
# matplot(new$x, cbind(pred.w.clim, pred.w.plim[,-1]),
#         lty = c(1,2,2,3,3), type = "l", ylab = "predicted y")
# ipred <- predict(model, interval = "prediction")
# iconf <- predict(model, interval = "confidence")
plot(y ~ x, data = d)
abline(model)
lines(new$x, pred.clim[, 2], pch = 19, col = "red")
lines(new$x, pred.clim[, 3], pch = 19, col = "red")
lines(new$x, pred.plim[, 2], pch = 19, col = "blue")
lines(new$x, pred.plim[, 3], pch = 19, col = "blue")

# using nls
formu <- as.formula(y ~ a + b *x)
nlsmodel <- nls(formu, start = list(a = 0, b = 1), data = d)
niter <- 200
# niter <- 2000
nlsboot <- nlsBoot(nlsmodel, niter = niter)
(pred.clim <- nlsBootpredict(nlsboot, newdata = new, interval = "confidence"))
(pred.plim <- nlsBootpredict(nlsboot, newdata = new, interval = "prediction"))

lines(new$x, pred.clim[, "2.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new$x, pred.clim[, "97.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new$x, pred.plim[, "2.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)
lines(new$x, pred.plim[, "97.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)
