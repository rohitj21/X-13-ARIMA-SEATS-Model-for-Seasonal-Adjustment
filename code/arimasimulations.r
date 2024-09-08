rARIMA <- function(n, freq = 12, Phi = NULL, d = 0, Theta = NULL, Sphi = NULL, D = 0, Stheta = NULL, sig2 = 1) {
  # Phi = [phi_1, phi_2, ...phi_p]
  # Theta = [theta_1, theta_2, ... theta_q]
  # Sphi = [seasonal phi_1, seasonal phi_2, ... seasonal phi_P]
  # Stheta = [seasonal theta_1, seasonal theta_2, ... seasonal theta_Q]
  
  p = length(Phi)
  q = length(Theta)
  P = length(Sphi)
  Q = length(Stheta)
  
  epsilon = if(q > 0) numeric(q) else NULL
  X_prev = if(p > 0) numeric(p) else NULL
  
  epsilon_s = if(Q > 0) numeric(Q) else NULL
  X_prev_s = if(P > 0) numeric(P) else NULL
  
  K = 1
  Xt = numeric(K + n)
  
  for (i in 1:length(Xt)) {
    epsilon_t = rnorm(1, mean = 0, sd = sqrt(sig2))
    Xt[i] = epsilon_t
    
    # Add the MA part
    if (q > 0) {
      Xt[i] = Xt[i] + sum(Theta * epsilon)
      epsilon = c(epsilon_t, epsilon[-q])
    }
    
    # Add the Seasonal MA part
    if (Q > 0 && i > freq) {
      Xt[i] = Xt[i] + sum(Stheta * epsilon_s)
      epsilon_s = c(Xt[i - freq], epsilon_s[-Q])
    }
    
    # Add the AR part
    if (p > 0) {
      Xt[i] = Xt[i] + sum(Phi * X_prev)
      X_prev = c(Xt[i], X_prev[-p])
    }
    
    # Add the Seasonal AR part
    if (P > 0 && i > freq) {
      Xt[i] = Xt[i] + sum(Sphi * X_prev_s)
      X_prev_s = c(Xt[i - freq], X_prev_s[-P])
    }
  }
  
  # Apply differencing (d)
  if (d > 0) {
    Xt <- diffinv(Xt, differences = d)
  }
  
  # Apply seasonal differencing (D)
  if (D > 0) {
    Xt <- diffinv(Xt, lag = freq, differences = D)
  }
  
  # Remove the initial K values
  N <- length(Xt)
  Xt <- Xt[-c(1:(N - n))]
  
  # Return the final time series with the correct frequency
  return(ts(Xt, frequency = freq))
}

# Example usage:
phi = c(0.2, 0.1)
theta = c(0.8, 0.8)
Sphi = c(0.1)  # Seasonal AR component
Stheta = c(0.4)  # Seasonal MA component
d = 1
D = 1  # Seasonal differencing
y <- rARIMA(1e2, freq = 12, Phi = phi, d = d, Theta = theta, Sphi = Sphi, D = D, Stheta = Stheta, sig2 = 0.1)
plot(y, type = 'l')

Visualise <-  function(Phi , d, Theta, sig2 = 1){
  
  y <- rARIMA(1e2,Phi = phi, d=d, Theta = theta, sig2 = 1)
  
  plot(0,0, xlim = c(1, 1e2), ylim = c(0.8*min(y), 1.2*max(y)), type = 'n')
  
  for(i in 1:10){
    y <- rARIMA(1e2,Phi = phi, d=d, Theta = theta, sig2 = 1)
    lines(y, col = 'grey')
  }
  
}

#Visualise(Phi = phi, d = d, Theta = theta)


#spectrum_result <- spectrum(y, plot = TRUE)

