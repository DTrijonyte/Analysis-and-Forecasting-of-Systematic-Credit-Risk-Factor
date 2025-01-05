
 # Description:   Re-fit final selected model via OLS (Post-Lasso with ARMA errors)
 #                auto.arima(z, xreg = as.matrix(LASSO_selected_variables), 
 #                         stationary = TRUE, seasonal = FALSE)
 
 # Attributes:
 #  zx_data -- combined dataset with z dependent variables and x features
 #  selected_variables -- prescreened list of variables or any set of variables < sample size
 #  t_threshold -- 1 for data mining 2 for validation of final models
 #  max_vars -- not more than variables to be in the final model to prevent overfitting
 #  trace   -- set to TRUE to print the elimination sequence
 
 function(zx_data, 
 selected_variables,
 t_threshold = 1,
 max_vars = round(dim(zx_data)[1]/8),
 trace = FALSE){
 t_exclude <- 0
 
 if(length(selected_variables) == 1){
 xreg_matrix <- as.matrix(zx_data[, selected_variables, drop = FALSE])
 arima_model <- auto.arima(zx_data$z, xreg = xreg_matrix, stationary = TRUE, 
 seasonal = FALSE, allowmean = FALSE, max.p = 4L, max.q = 2L)
 }else{
 while(((t_exclude < t_threshold) & 
 (length(selected_variables) > 1)) | 
 (length(selected_variables) > max_vars)){
 selected_variables <- setdiff(selected_variables, names(t_exclude))
 xreg_matrix <- as.matrix(zx_data[, selected_variables, drop = FALSE])
 arima_model <- auto.arima(zx_data$z, xreg = xreg_matrix, stationary = TRUE, 
 seasonal = FALSE, allowmean = FALSE, max.p = 4L, max.q = 2L)
 arima_model$var.coef %>% diag() %>% sqrt() -> coef_se #get s.e. of errors
 t_stat <- (coef(arima_model) / coef_se) %>% #  make t-stat
 {.[selected_variables]} %>%         # ignore non xreg
 abs()
 t_exclude <- t_stat[which.min(t_stat)]  # exclude if < 1 (or any selected threshold) refit arima
 if(trace) print(t_exclude)
 }
 }
 
 return(arima_model)
 }
