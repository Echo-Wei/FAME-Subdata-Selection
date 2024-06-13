SIS <- function(sample_size, x, y, subdata_num, dim, random_seed, index_beta){
  ## Given original dataset, return MSE of sparse regression on subdata selected by SIS method
  # Param sample_size: a numeric value on the range of 0 to 1 exclusive
  # Param x: continous predictor variables
  # Param y: continous response variable
  # Param subdata_num: an integer representing the number of data points in the subdata
  # Param dim: an integer representing the number of columns in the subdata
  # Param random_seed: the seed for generating the random number
  # Param index_beta: the location of the true predictors
  
  # standardize predictor variables and response variable
  x = scale(x); y = scale(y)
  n = nrow(x); p = ncol(x)
  smp_siz = floor(sample_size*n) 
  
  # parameter testing
  if ((subdata_num >= n) | (dim >= p)){
    print('The size of the subdata should be smaller than that of the original data.')
    break
  }
  
  if ((sample_size >= 1) | (sample_size <= 0)){
    print('The sample_size should be a value between 0 and 1 exclusively.')
    break
  }
  
  # split data into training and testing dataset
  set.seed(random_seed)
  train_ind = sample(seq_len(n),size = smp_siz)         
  x.train <- x[train_ind,]                   
  x.train <- cbind(x.train, c(1:smp_siz))
  x.test <- x[-train_ind,]
  y.train <- y[train_ind]
  y.test <- y[-train_ind]
  de <- ncol(x.train)
  
  # Sparse regression and prediction
  t1 <- Sys.time()
  w = abs(t(x.train[,-de])%*%y.train)              
  w.order <- order(w, decreasing = TRUE)     
  z.sis <- x.train[,w.order[1:dim]]
  cv.sis <- cv.glmnet(z.sis, y.train, alpha = 1, lambda = NULL)    
  lambda.sis <- cv.sis$lambda.min                              
  lasso.sis <- glmnet(z.sis, y.train, alpha = 1, lambda = lambda.sis)  
  predict.sis <- predict(lasso.sis, s = lambda.sis, newx = x.test[,w.order[1:dim]])  
  t2 <- Sys.time()
  
  # measurement
  mse.sis <- mean((y.test - predict.sis)^2)  
  choose_beta <- w.order[which(lasso.sis$beta != 0)]
  tp <- sum(choose_beta %in% index_beta); fp <- length(which(lasso.sis$beta != 0))- tp
  fn <- tr - tp; tn <- p - tr - fp
  tp.sis <- tp/(tp + fn)
  tn.sis <- tn/(tn+fp)                     # fp and fn rate 
  t.sis <- t2 - t1
  
  var_estimate.sis <- var(predict.sis)
  var.sis <- var(y.test)
  snr.sis <- var(predict.sis)/var(y.test - predict.sis)
  
  output_list <- list("mse" = mse.sis, "lasso_model" = lasso.sis, "lasso_lambda" = lambda.sis, "train_index" = train_ind,
                      "tp" = tp.sis, "tn" = tn.sis, "time" = t.sis, "prediction_variance" = var_estimate.sis,
                      "test_variance" = var.sis, 'snr' = snr.sis)
  
  return(output_list)
}