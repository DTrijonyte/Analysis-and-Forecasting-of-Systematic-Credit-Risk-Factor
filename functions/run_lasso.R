
 # Description:   Performs variable selection using LASSO regression with Leave-One-Out Cross-Validation (LOOCV) on a subset of the data
 
 # Arguments
 #  data_subset           -- a data frame or matrix containing the subset of the data from which variables will be selected
 #  z_subset              -- a vector or matrix containing the response variable for the subset of the data
 #  prescreened_variables -- a vector of variable names that have been pre-screened for inclusion in the LASSO regression
 
 function(data_subset, z_subset, prescreened_variables) {
 prescreened_data <- data_subset[, prescreened_variables]
 prescreened_data <- as.matrix(prescreened_data)
 z_subset <- as.matrix(z_subset)
 
 # Identify complete cases in prescreened_data
 cc <- complete.cases(prescreened_data)
 
 # Subset prescreened_data based on complete cases
 prescreened_data <- prescreened_data[cc, , drop = FALSE]
 z_subset <- z_subset[cc, , drop = FALSE]  # Ensure z_subset is also subset accordingly
 
 N <- length(z_subset) # Apply LOOCV, no seed dependency
 lasso_model <- cv.glmnet(x = prescreened_data, y = z_subset, alpha = 1, 
 type.measure = "mse", nfolds = N, foldid = seq(N),
 grouped = FALSE)
 lasso_coefficients <- coef(lasso_model)
 selected_variables <- which(lasso_coefficients[, 1] != 0)
 selected_variable_names <- rownames(lasso_coefficients)[selected_variables]
 
 return(selected_variable_names[-1])
 }
