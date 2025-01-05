
 # Description:   Screens the features based on one of the approaches, then reduces with LASSO.
 
 # Attributes:
 #  input_data -- all scaled data frame
 #  segment    -- either "LE" for Legal entities, or "PI" for private individuals
 #  features   -- any reasonable subset of features, default is G1 group with predictions
 #  lag_max    -- the maximum number of lags for the FLI effects of features
 #  end_date   -- end date for screening
 
 #  Returns list of features by different approaches:
 #   forward for block forward, best-subset only if (max_lag + 1)^p (disabled),
 #   where p is the number of features is below 10^5, correlation for correlation based and
 #   all-in to allow for the same variable lags in LASSO selection
 
 function(input_data = all_scaled_data,
 segment = "PI",
 features  = c(G1),      
 lag_max = 3L,
 end_date = c(2019, 4)){
 # Take the raw dependent variable
 P <- length(features)
 # Note: Smoothing of z disconnects it from macro cycles
 z_var_name <- paste0("Z_", segment)
 z <- na.trim(input_data[, z_var_name], sides = "both") 
 z_start_date <- start(z)
 
 # 1. Add lags to features up to lag_max quarters
 dates <- time(input_data)
 selected_data <- input_data[, features] 
 all_scaled_df <- as.data.frame(selected_data)
 
 lagged_data <- data.frame(lapply(all_scaled_df, create_lagged_vars, lag_max = lag_max))
 combined_data <- cbind.data.frame(dates, all_scaled_df, lagged_data)
 
 # 2. Pre-screening window(z ends at 2019 Q4 or any relevant end_date):
 # Select subsets:
 
 z_subset <- data.frame(z = window(z, start = z_start_date, end = end_date))  
 scaled_data_subset <- data.frame(window(selected_data, start = z_start_date, end = end_date))
 start_yearqtr <- paste(z_start_date, collapse = " Q")
 end_yearqrt   <- paste(end_date, collapse = " Q")
 combined_data_subset <- subset(combined_data[, -1], as.yearqtr(dates) >= start_yearqtr &
 as.yearqtr(dates) <= end_yearqrt)
 
 zx_data <- cbind(z_subset, combined_data_subset)
 colnames(zx_data) <- c("z", colnames(combined_data_subset))
 
 # A. Hard-thresholding: use feature's lag with max(abs()) correlation
 
 prescreened_by_correlation <- NULL  
 search_names <- names(combined_data_subset)
 for (i in 1:P) {
 variable_group <- search_names[grepl(paste0("^", features[i]), search_names)]
 cor(cbind(z_subset, combined_data_subset[,variable_group]), 
 use = "pairwise.complete.obs")[1,-1] %>%
 abs() %>%
 which.max() %>%
 names() -> prescreened_by_correlation[i]
 }
 
 # B. Block-forward selection, at each step find one feature including lags
 #       that increases r2 the most, then remove this variable with all its
 #       lags from search space and move to next variable selection lm(z ~ x1.lag1 + ...)
 
 search_space <- search_names
 prescreened_by_block_forward <- NULL
 for(i in 1:P){
 P_lags <- length(search_space)
 if(i==1){step <- matrix(search_space, nrow = P_lags)
 }else{step <- cbind(matrix(rep(prescreened_by_block_forward, P_lags),
 ncol = i - 1L, byrow = TRUE), search_space)}
 step %>%
 apply(1, function(x){
 summary(lm(z ~ ., zx_data[,c("z", x)]))$r.squared
 }) %>%
 which.max() -> selection
 selected_variable <- step[selection, dim(step)[2]]
 prescreened_by_block_forward[i] <- selected_variable
 exclude_from_search <- search_space[grepl(paste0("^", unlist(strsplit(selected_variable, "[.]"))[1]), 
 search_space)]
 search_space <- setdiff(search_space, exclude_from_search)
 }
 
 ## C. Best subset version of selection
 
 #if((1 + lag_max)^P <= 10^4){selected_model <- select_model() # slow to compute!
 #prescreened_best_subset <- names(selected_model$coefficients)[-1]}else{
 #  prescreened_best_subset <- NULL
 #}
 
 # 3. up to 24 pre-screened variables are used then for LASSO approach, where lambda is set
 #    via cross-validation -- selection step
 
 correlation   <- run_lasso(combined_data_subset, z_subset, prescreened_by_correlation)
 block_forward <- run_lasso(combined_data_subset, z_subset, prescreened_by_block_forward)
 #if(!is.null(prescreened_best_subset)){
 #  best_subset <- run_lasso(combined_data_subset, z_subset, prescreened_best_subset)}else{
 #  best_subset <- NULL
 #}
 
 all_in        <- run_lasso(combined_data_subset, z_subset, search_names)
 
 return(list(correlation = correlation, 
 block_forward = block_forward, 
 #best_subset = best_subset,
 ps_cor = prescreened_by_correlation,
 ps_bf = prescreened_by_block_forward,
 all_in = all_in,
 z_start_date = z_start_date,
 z = z,
 zx_data = zx_data,
 combined_data = combined_data,
 dates = dates))
 }  
