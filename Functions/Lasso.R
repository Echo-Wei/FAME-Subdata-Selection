LASSO <- function(sample_size, x, y, subdata_num, dim, random_seed, index_beta){
  ## Given original dataset, return MSE of sparse regression on subdata selected by IBOSS method
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
  cv.lasso <- cv.glmnet(x.train[,-de], y.train, alpha = 1, lambda = NULL)     # Cross Validation Lasso
  lambda.lasso <- cv.lasso$lambda.min                               # Best Lambda
  lassolm <- glmnet(x.train[,-de], y.train, alpha = 1, lambda = lambda.lasso)   # Rebuild the model
  predict.lasso <- predict(lassolm, s = lambda.lasso, newx = x.test)       # Prediction
  t2 <- Sys.time()
  
  # measurement
  mse.lasso <- mean((y.test - predict.lasso)^2)        # MSE
  choose_beta <- which(lassolm$beta != 0)
  tp <- sum(choose_beta %in% index_beta); fp <- length(which(lassolm$beta != 0))- tp
  fn <- tr - tp; tn <- p - tr - fp
  tp.lasso <- tp/(tp + fn)
  tn.lasso <- tn/(tn+fp)                    # fp and fn rate  
  t.lasso <- t2 - t1
  
  var_estimate.lasso <- var(predict.lasso)
  var.lasso <- var(y.test)
  snr.lasso <- var(predict.lasso)/var(y.test - predict.lasso)
  
  # output
  output_list <- list("mse" = mse.lasso, "lasso_model" = lassolm, "lasso_lambda" = lambda.lasso, "train_index" = train_ind,
                      "tp" = tp.lasso, "tn" = tn.lasso, "time" = t.lasso, "prediction_variance" = var_estimate.lasso,
                      "test_variance" = var.lasso, 'snr' = snr.lasso)
  
  return(output_list) 
}