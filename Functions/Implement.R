library(MASS)

FAME_mse <- c()

snr.all <- c()
for (loop in c(1:50)){
  print(loop)
  set.seed(loop)
  n = 2000; p = 1000; rho = 0.5
  mu <- rnorm(p, 0, 4)            # Mean
  sigma <- diag(rep(1-rho, p)) + matrix(rho, p, p) 
  # sigma <- diag(rep(1-rho, p))
  
  # Set the first five coefficients as true coefficients, others as 0
  beta <- rep(0, p)
  tr <- 5
  index_beta <- c(1:tr)
  beta[index_beta] <- abs(rnorm(tr, 0, 10))          # Sampling beta with 5 non-zero coefficients
  
  # Generate the response variable and predictor variables
  x = mvrnorm(n, mu, sigma)      
  error = rnorm(n, 0, 4)  
  y = x%*%beta + error  
  snr.tmp <- var(x%*%beta)/var(y - x%*%beta)
  snr.all <- c(snr.all, snr.tmp)
  data = cbind(y, x)
  
  # FAME
  FAME_fit <- FAME(0.7, x, y, 200, 5, loop, index_beta, tr)
  loop_MSE <- FAME_fit$mse
  FAME_mse <- c(FAME_mse, loop_MSE)
}