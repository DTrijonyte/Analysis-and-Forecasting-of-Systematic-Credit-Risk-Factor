 # Description:   Apply bivariate analysis for LM and regression with ARMA errors

 # Arguments:
 # eda_output   -- information with exploratory analysis for LE or PI
 # lag_max      -- lag_max parameter with default set to 5 quarters
 # segment      -- string name of the segment either "le" or "pi"
 
 function(eda_output, lag_max = 5L, segment = "le") {
 # Find max columns for LE and PI
 eda_output$max_column <- max_column(eda_output)
 
 # Taken from prescreen_features():
 # Create lagged data
 all_scaled_dataframe <- as.data.frame(all_scaled_data[,-c(1:2)])
 lagged_data <- data.frame(lapply(all_scaled_dataframe, create_lagged_vars, 
 lag_min = 0, lag_max = lag_max))
 dates <- time(all_scaled_data)
 lagged_data <- cbind.data.frame(dates, lagged_data)
 
 # Prepare Z data
 z_var_name <- paste0("Z_", toupper(segment))
 z <- na.trim(all_scaled_data[, z_var_name], sides = "both") 
 z_start_date <- start(z)
 dates <- time(z)
 z_with_dates <- cbind.data.frame(dates, z)
 
 # Combine Z data with lagged data
 combined_data <- z_with_dates %>%
 right_join(lagged_data, by = "dates") %>%
 arrange(dates) %>%
 select(-1)
 
 # Select columns for LE and PI
 variables <- gsub("z.lag0", "z", eda_output$max_column)
 selected_columns <- combined_data[, variables]
 
 # Perform LM and ARIMA analysis
 lm_regressions_results    <- lm_analysis(selected_columns)
 arima_regressions_results <- arima_analysis(selected_columns)
 
 # Combine LM and ARIMA results
 bivariate_regression <- merge(arima_regressions_results, lm_regressions_results, by = "Variable")
 rownames(bivariate_regression) <- bivariate_regression$Variable
 bivariate_regression <- bivariate_regression[, -1]
 
 return(bivariate_regression)
 }
