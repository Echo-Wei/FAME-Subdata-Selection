library(glmnet)
FAME <- function(sample_size, x, y, subdata_num, dim, random_seed, index_beta, tr){
  ## Given original dataset, return MSE of sparse regression on subdata selected by the proposed FAME algorithm
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
  
  # parameter testing
  if ((subdata_num >= n) | (dim >= p)){
    print('Warning: The size of the subdata should be smaller than that of the original data.')
  }

  if ((sample_size >= 1) | (sample_size <= 0)){
    print('Warning: The sample_size should be a value between 0 and 1 exclusively.')
  }
  
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
  
  # Calculate the score index
  w = c()
  for(i in c(1 : (dim(x.train)[2]-1))){
    model = lm(y.train ~ x.train[,i])
    model_summary = summary(model)
    model_coefficient = model_summary$coefficients
    t_value_x = model_coefficient[2, 3]
    w <- c(w, abs(t_value_x))
  }
  w.order <- order(w, decreasing = TRUE)
  
  for (i in 1:dim){
    if (i == 1){
      x.temp <- x.train[ , c(w.order[i], de)]                
    } else {
      x.temp <- x.train[-x.index, c(w.order[i], de)]         
    }
    x.order <- order(x.temp[,1], decreasing = TRUE)
    x.length <- length(x.order)
    x.index <- c(x.index, c(x.temp[x.order[1:r_top], 2], x.temp[x.order[(x.length-r_bot+1):x.length], 2]))  
  }                                 
  z.fame <- (x.train[x.index, w.order[1:dim]])      
  y.fame <- y.train[x.index]
  
  # Sparse regression and prediction
  t1 <- Sys.time()
  cv.fame <- cv.glmnet(z.fame, y.fame, alpha = 1, lambda = NULL)   
  lambda.fame <- cv.fame$lambda.min                              
  lasso.fame <- glmnet(z.fame, y.fame, alpha = 1, lambda = lambda.fame)   
  predict.fame <- predict(lasso.fame, s = lambda.fame, newx = x.test[,w.order[1:dim]])   
  t2 <- Sys.time()
  
  # measurement
  mse.fame <- mean((y.test - predict.fame)^2)
  choose_beta <- w.order[which(lasso.fame$beta != 0)]
  tp <- sum(choose_beta %in% index_beta); fp <- length(which(lasso.fame$beta != 0))- tp
  fn <- tr - tp; tn <- p - tr - fp
  tp.fame <- tp/(tp + fn)
  tn.fame <- tn/(tn+fp)                 # fp and fn rate
  t.fame <- t2 - t1
  
  var_estimate.fame <- var(predict.fame)
  var.fame <- var(y.test)
  snr.fame <- var(predict.fame)/var(y.test - predict.fame)
  
  output_list <- list("mse" = mse.fame, "lasso_model" = lasso.fame, "lasso_lambda" = lambda.fame, "train_index" = train_ind,
                      "tp" = tp.fame, "tn" = tn.fame, "time" = t.fame, "prediction_variance" = var_estimate.fame,
                      "test_variance" = var.fame, 'snr' = snr.fame)

  return(output_list) 
}