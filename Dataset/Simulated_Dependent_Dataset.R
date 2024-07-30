# Generate Dependent Data Normal --------------------------------------------------------------------
library(MASS)
set.seed(12300)
n = 2000; p = 1000; rho = 0.5
mu <- rnorm(p, 0, 4)            # Mean
sigma <- diag(rep(1-rho, p)) + matrix(rho, p, p)        

# Set the first five coefficients as true coefficients, others as 0
beta <- rep(0, p)
tr <- 5
index_beta <- c(1:tr)
beta[index_beta] <- abs(rnorm(tr, 0, 10))            # Sampling beta with 5 non-zero coefficients

# Generate the response variable and predictor variables
x = mvrnorm(n, mu, sigma)      
error = rnorm(n, 0, 4)                        
y = x%*%beta + error          
data = cbind(y, x)
