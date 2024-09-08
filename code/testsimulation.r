library(seasonal)
library(stats)

t <- 1:120
S <- rep(rnorm(12, 0, 2), 10)
M <- ((t/120)**3 - 2*(t/120)**2 + 2*(t/120))*50
I <- numeric(120)
theta <- 0.7
er <- rnorm(1)
for(i in 1:120){
  e <- rnorm(1)
  I[i] = theta*er + e
  er <- e
}
true_seas_adj <- ts(M + I, start = c(2011, 1), frequency = 12)
X <- ts(M + S + I, start = c(2011, 1), frequency = 12)
m <- seas(X)
plot(m, main = "Seasonal Adjustment", ylab = "Value")
lines(true_seas_adj, col = 'green')
legend("bottomright", legend = c("Original Data", "Seasonally Adjusted", "True Seasonal Adjustment"), 
       col = c("blue", "black", "green"), lty = 1)


est_trend <- m$series$s12
plot(est_trend, main = "Estimated Trend vs. True Trend", ylab = "Value", col = 'blue', type = 'l')
lines(ts(M, start = c(2011, 1), frequency = 12), col = 'red')
legend("bottomright", legend = c("Estimated Trend", "True Trend"), col = c("blue", "red"), lty = 1)

