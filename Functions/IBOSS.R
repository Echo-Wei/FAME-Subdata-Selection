IBOSS <- function(sample_size, x, y, subdata_num, dim, sim_num){
  ## Given original dataset, return MSE of sparse regression on subdata selected by IBOSS method
  # Param sample_size: a numeric value on the range of 0 to 1 exclusive
  # Param x: continous predictor variables
  # Param y: continous response variable
  # Param subdata_num: an integer representing the number of data points in the subdata
  # Param dim: an integer representing the number of columns in the subdata
  # Param sim_num: the number of simulations
  
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
  
  # Calculate the number of data points extracted from each column
  r = ceiling(subdata_num/dim)                      
  if(r %% 2 == 0){
    r_top <- r /2
    r_bot <- r/2
  } else {
    r_top <- (r+1)/2
    r_bot <- r - r_top
  }
  
  mse.iboss <- c()
  for (loop in c(1:sim_num)){
    # split data into training and testing dataset
    set.seed(loop)
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
    set.seed(100 + loop)
    dim.iboss = sample(1:(p-1), dim, replace=F)
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
      
    z.iboss <- x.train[x.index, -de]      
    y.iboss <- y.train[x.index]
    cv.iboss <- cv.glmnet(z.iboss, y.iboss, alpha = 1, lambda = NULL)     
    lambda.iboss <- cv.iboss$lambda.min                               
    lasso.iboss <- glmnet(z.iboss, y.iboss, alpha = 1, lambda = lambda.iboss)  
    predict.iboss <- predict(lasso.iboss, s = lambda.iboss, newx = x.test)       
      
    mse.iboss <- c(mse.iboss, mean((y.test - predict.iboss)^2))        
    }
  return(mean(mse.iboss)) 
}