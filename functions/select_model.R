
 # Description:   Performs variable selection with lagged variables for time series data by evaluating models based on the adjusted R-squared value

 # Attributes:
 #  data            -- data frame containing the time series data
 #  variable_names  -- vector of variable names for which lagged variables will be created and evaluated
 #  max_lag         -- maximum number of lags to consider for each variable
 
 function(data = zx_data,
 variable_names = features,
 max_lag = lag_max){
 models <- NULL
 N <- length(variable_names)
 M <- max_lag + 1L
 for(index in 1:N){
 step <- c(variable_names[index], paste(variable_names[index],".lag",1:max_lag, sep = ""))
 models <- cbind(models, rep(step, each = M^(N-index), times =M^(index-1L) ))
 }
 adj_r2 <- apply(models,1, function(x){
 summary(lm(z~., data = data[c("z", x)]))$adj.r.squared})
 selection <- which.max(adj_r2)
 return(lm(z~., data = data[c("z", models[selection,])]))
 }
