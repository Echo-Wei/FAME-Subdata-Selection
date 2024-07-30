SIS <- function(sample_size, x, y, dim, random_seed, index_beta, tr){
  ## Given original dataset, return MSE of linear regression on subdata selected by SIS method
  # Param sample_size: a numeric value on the range of 0 to 1 exclusive
  # Param x: continous predictor variables
  # Param y: continous response variable
  # Param dim: an integer representing the number of columns in the subdata
  # Param random_seed: the seed for generating the random number
  # Param index_beta: the location of the true predictors
  # Param tr: the number of true coefficients
  
  # standardize predictor variables
  x = scale(x); 
  n = nrow(x); p = ncol(x)
  smp_siz = floor(sample_size*n) 
  
  # split data into training and testing dataset
  set.seed(random_seed)
  train_ind = sample(seq_len(n),size = smp_siz)         
  x.train <- x[train_ind,]                   
  x.train <- cbind(x.train, c(1:smp_siz))
  x.test <- x[-train_ind,]
  y.train <- y[train_ind]
  y.test <- y[-train_ind]
  de <- ncol(x.train)
  
  # Multivariate linear regression and prediction
  w = abs(t(x.train[,-de])%*%y.train)              
  w.order <- order(w, decreasing = TRUE)     
  z.sis <- x.train[,w.order[1:dim]]
  t1 <- Sys.time()
  lm.sis <- lm(y.train ~ z.sis) 
  z.sis <- x.test[,w.order[1:dim]]
  predict.sis <- predict(lm.sis, as.data.frame(z.sis))  
  t2 <- Sys.time()
  
  # measurement
  mse.sis <- mean((y.test - predict.sis)^2)  
  choose_beta <- w.order[which(lm.sis$coefficients != 0)]
  tp <- sum(choose_beta %in% index_beta); fp <- length(which(lm.sis$coefficients != 0))- tp
  fn <- tr - tp; tn <- p - tr - fp
  tp.sis <- tp/(tp + fn)
  tn.sis <- tn/(tn+fp)                     # fp and fn rate 
  t.sis <- t2 - t1
  
  var_estimate.sis <- var(predict.sis)
  var.sis <- var(y.test)
  snr.sis <- var(predict.sis)/var(y.test - predict.sis)
  
  output_list <- list("mse" = mse.sis, "lasso_model" = lm.sis, "train_index" = train_ind,
                      "tp" = tp.sis, "tn" = tn.sis, "time" = t.sis, "prediction_variance" = var_estimate.sis,
                      "test_variance" = var.sis, 'snr' = snr.sis)
  
  return(output_list)
}