
 # Description:   Estimate bivariate regression with up to ARIMA(2, 0, 1) errors
 
 # Arguments:
 # combined_data   -- data frame with all data in one place
 
 function(combined_data) {
 variables <- setdiff(colnames(combined_data), "z")
 y <- combined_data[, "z"]
 results_list_arima <- list()
 for (i in seq_along(variables)) {
 x <- combined_data[, variables[i]]
 arima_model <- auto.arima(y, xreg = x, stationary = TRUE, seasonal = FALSE, 
 max.p = 2L, max.q = 1L)
 results_list_arima[[i]] <- data.frame(
 Variable = variables[i],
 Coefficient_ar1 = arima_model$coef["ar1"],
 Coefficient_ar2 = arima_model$coef["ar2"],
 Coefficient_ma1 = arima_model$coef["ma1"],
 Coefficient_xreg = arima_model$coef["xreg"],
 row.names = i
 )
 }
 arima_regressions <- do.call(rbind, results_list_arima)
 return(arima_regressions)
 }
