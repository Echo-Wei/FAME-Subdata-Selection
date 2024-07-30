library(glmnet)
ScoreIndex <- function(x, y, family = 'continuous'){
  ## Given original dataset, return the order of predictors ranked by score index in descending order 
  # Param x: continous predictor variables
  # Param y: continous response variable
  # Param family: default value is continuous
  
  if (family == 'continuous'){
    w = c()
    for(i in c(1 : dim(x)[2])){
      model = lm(y ~ x[,i])
      model_summary = summary(model)
      model_coefficient = model_summary$coefficients
      t_value_x = model_coefficient[2, 3]
      w <- c(w, abs(t_value_x))
    }
    w.order <- order(w, decreasing = TRUE)
  } else if(family == 'experimental'){
    p = dim(x)[2]
    w <- matrix(NA, p, 1)
    for (index_col in 1:p){
      dat <- data.frame(Cov_value = y, Group = x[,index_col])
      w_tmp <- aggregate(Cov_value ~ Group, data = dat, mean)
      w[index_col] <- abs(w_tmp[2, 2] - w_tmp[1, 1])
    }
    w.order <- order(w, decreasing = TRUE)     
  } else if(family == 'calssification'){
    p = dim(x)[2]
    w <- matrix(NA, p, 1)
    for (index_col in 1:p){
      dat <- data.frame(Cov_value = x[,index_col], Group = y[,])
      w_tmp <- aggregate(Cov_value ~ Group, data = dat, mean)
      w[index_col] <- abs(w_tmp[2, 2] - w_tmp[1, 1])
    }
    w.order <- order(w, decreasing = TRUE)
  }
  
  return(w.order) 
}