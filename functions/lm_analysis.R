
 # Description:   Estimate bivariate regression

 # Arguments:
 # combined_data   -- data frame with all data in one place
 
 function(combined_data) {
 variables <- setdiff(colnames(combined_data), "z")
 y <- combined_data[, "z"]
 results_list_lm <- list()
 for (i in seq_along(variables)) {
 x <- combined_data[, variables[i]]
 model <- lm(y ~ x)
 results_list_lm[[i]] <- data.frame(
 Variable = variables[i],
 Coefficient = coef(model)[2],
 r_squared = summary(model)$r.squared,
 p_value = summary(model)$coefficients[2, "Pr(>|t|)"],
 stringsAsFactors = FALSE
 )
 }
 lm_regressions <- do.call(rbind, results_list_lm)
 return(lm_regressions)
 }
