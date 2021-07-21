### Test nlsBootPredict() on a linear model
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
(pred.clim <- nlsBootPredict(nlsboot, newdata = new, interval = "confidence"))
(pred.plim <- nlsBootPredict(nlsboot, newdata = new, interval = "prediction"))

lines(new$x, pred.clim[, "2.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new$x, pred.clim[, "97.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new$x, pred.plim[, "2.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)
lines(new$x, pred.plim[, "97.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)

### Use of nlsBootPredict() a non linear model 
data(vmkm)
nls1 <- nls(michaelis,vmkm,list(Km=1,Vmax=1))
plotfit(nls1, smooth = TRUE)

nlsb1 <- nlsBoot(nls1, niter = niter)
new1 <- data.frame(S = seq(0.3, 2, 0.1))
(pred.clim <- nlsBootPredict(nlsb1, newdata = new1, interval = "confidence"))
(pred.plim <- nlsBootPredict(nlsb1, newdata = new1, interval = "prediction"))

lines(new1$S, pred.clim[, "2.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new1$S, pred.clim[, "97.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new1$S, pred.plim[, "2.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)
lines(new1$S, pred.plim[, "97.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)

### Use of nlsBootPredict() on a model with two independent variables
data(vmkmki)
nls2 <- nls(compet_mich, vmkmki, list(Km=1,Vmax=20,Ki=0.5))
plotfit(nls2, variable=1)
nlsb2 <- nlsBoot(nls2, niter = niter)
# Define a new dataframe with a fixed value of one of the independent variable
new2 <- data.frame(S = seq(0, 200, length.out = 50), I = rep(20, 50))
(pred.clim <- nlsBootPredict(nlsb2, newdata = new2, interval = "confidence"))
(pred.plim <- nlsBootPredict(nlsb2, newdata = new2, interval = "prediction"))

plot(new2$S, pred.clim[, "Median"], type = "l", col = "black", lwd = 2, lty = 2)
lines(new2$S, pred.clim[, "2.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new2$S, pred.clim[, "97.5%"], pch = 19, col = "red", lwd = 2, lty = 2)
lines(new2$S, pred.plim[, "2.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)
lines(new2$S, pred.plim[, "97.5%"], pch = 19, col = "blue", lwd = 2, lty = 2)

