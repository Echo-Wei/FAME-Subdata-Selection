IBOSS_LASSO <- function(sample_size, x, y, subdata_num, dim, random_seed, index_beta, tr){
  ## Given original dataset, return MSE of sparse regression on subdata selected by IBOSS method
  # Param sample_size: a numeric value on the range of 0 to 1 exclusive
  # Param x: continous predictor variables
  # Param y: continous response variable
  # Param subdata_num: an integer representing the number of data points in the subdata
  # Param dim: an integer representing the number of columns in the subdata
  # Param random_seed: the seed for generating the random number
  # Param index_beta: the location of the true predictors
  # Param tr: the number of true coefficients
  
  # standardize predictor variables
  x = scale(x); 
  n = nrow(x); p = ncol(x)
  smp_siz = floor(sample_size*n) 
  
  # Calculate the number of data points extracted from each column
  r = ceiling(subdata_num/dim)                      
  if(r %% 2 == 0){
    r_top <- r /2
    r_bot <- r/2
  } else {
    r_top <- (r+1)/2
    r_bot <- r - r_top
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
  
  # Generate subdata
  x.index <- c()
  z.index <- c()
  set.seed(random_seed)
  dim.iboss = sample(1:p, dim, replace=F)
  for (i in 1:dim){
    if (i == 1){
      x.temp <- x.train[ , c(dim.iboss[i], de)]                
    } else {
      x.temp <- x.train[-x.index, c(dim.iboss[i], de)]         
    }
    x.order <- order(x.temp[,1], decreasing = TRUE)
    x.length <- length(x.order)
    x.index <- c(x.index, c(x.temp[x.order[1:r_top], 2], x.temp[x.order[(x.length-r_bot+1):x.length], 2]))  # Rows index matching original dataset
  }                          
  
  # Sparse regression and prediction
  z.iboss <- x.train[x.index, -de]      
  y.iboss <- y.train[x.index]
  t1 <- Sys.time()
  cv.iboss <- cv.glmnet(z.iboss, y.iboss, alpha = 1, lambda = NULL)     
  lambda.iboss <- cv.iboss$lambda.min                               
  lasso.iboss <- glmnet(z.iboss, y.iboss, alpha = 1, lambda = lambda.iboss)  
  predict.iboss <- predict(lasso.iboss, s = lambda.iboss, newx = x.test)
  t2 <- Sys.time()
  
  # measurement
  mse.iboss <- mean((y.test - predict.iboss)^2)  
  choose_beta <- which(lasso.iboss$beta != 0)
  tp <- sum(choose_beta %in% index_beta); fp <- length(which(lasso.iboss$beta != 0))- tp
  fn <- tr - tp; tn <- p - tr - fp
  tp.iboss <- tp/(tp + fn)
  tn.iboss <- tn/(tn+fp)                    # fp and fn rate
  t.iboss <- t2 - t1
  
  var_estimate.iboss <- var(predict.iboss)
  var.iboss <- var(y.test)
  snr.iboss <- var(predict.iboss)/var(y.test - predict.iboss)
  
  output_list <- list("mse" = mse.iboss, "lasso_model" = lasso.iboss, "lasso_lambda" = lambda.iboss, "train_index" = train_ind,
                      "tp" = tp.iboss, "tn" = tn.iboss, "time" = t.iboss, "prediction_variance" = var_estimate.iboss,
                      "test_variance" = var.iboss, 'snr' = snr.iboss)
  
  return(output_list) 
}