
 # Description:   Computes Z factors cyclical parameters, performs MCMC simulations and checks stationarity and normality

 
 rm(list = ls())
 gc()
 
 options(max.print = 200, warn = 1)
 
 path <- dirname(rstudioapi::getSourceEditorContext()$path)
 setwd(path)
 
 # Load required libraries
 source("@CCFunctions.R")
 
 set.seed(2024) # un-comment (Ctrl + Shift + C) to replicate simulation
 CN    <- "LT"
 nsim  <- 10000L  # number of time series to simulate
 
 PL    <- 6   # lower limit in quarters
 PU    <- 80  # upper limit in quarters up to 20 years cycles (~ ECB WGEM paper (2018))
 SN    <- 'LE' # segment
 YQ <- "2024 Q3"
 results_path <- paste0("../4_Results/",YQ)
 output_file_name <- paste("CC_macro_cycles_", CN,".Rdata", sep = "")
 setwd(results_path)
 load(output_file_name)
 setwd(path)
 
 all_scaled_data[, paste0("Z_", SN)] %>% na.omit() -> y
 sample_size <- length(y) # the length of the data as in input dataset
 
 stochastic.cycle(y, fixed = FALSE, method = "Bayes", cyc = 40, thin = 4, 
 n.sample = nsim, burn = nsim/10) -> temp_res
 
 parameters <- data.frame(
 Parameter = c("$\\rho$", "$\\lambda$", "$\\rho_c$", "$
lpha$", "$\\sigma^2$", "$T$, years"),
 CI_90_lower = c(round(temp_res$interval[1, ], 3), round(pi/temp_res$interval[2, 2]/2, 3)),
 Mean = c(round(temp_res$params[1, 1:5], 3),  round(pi/temp_res$params[1, 2]/2, 3)),
 CI_90_upper = c(round(temp_res$interval[2, ], 3), round(pi/temp_res$interval[1, 2]/2, 3)),
 row.names = NULL
 )
 
 names(parameters)[names(parameters) == "CI_90_lower"] <- "CI 90% lower"
 names(parameters)[names(parameters) == "CI_90_upper"] <- "CI 90% upper"
 # 
 # 
 
 results_path <- "E:/user_DTrijonyte/Credit Cycle Model"
 file_name <- file.path(results_path, paste0("parameters_", CN, "_", SN, ".xlsx"))
 write.xlsx(parameters, file = file_name, rowNames = FALSE, encoding = "UTF-8")
 
 
 
 # Simulate time series
 sim_data <- simulate.cts(n = 710, plot = TRUE, seed = NULL,
 rho = temp_res$params[1,1],
 lambda = temp_res$params[1,2],
 rhoc = temp_res$params[1,3],
 alfa = temp_res$params[1,4],
 eps_irr_sd_2 = temp_res$params[1,5])
 matrix(c(temp_res$params), nrow = 2)[1,] -> params
 
 # Simulate all 10000 time series
 1:nsim %>% lapply(function(x)simulate.cts(n = 710, plot = FALSE, seed = NULL,
 rho = temp_res$params[1,1],
 lambda = temp_res$params[1,2],
 rhoc = temp_res$params[1,3],
 alfa = temp_res$params[1,4],
 eps_irr_sd_2 = temp_res$params[1,5])) -> temp_list
 
 
 
 
 
 # Define the test functions for stationarity and normality
 all_tests <- list(
 ADF = function(x) tseries::adf.test(x)$p.value,
 PP = function(x) tseries::pp.test(x)$p.value,
 KPSS = function(x) tseries::kpss.test(x)$p.value,
 AD = function(x) nortest::ad.test(x)$p.value,
 SW = function(x) stats::shapiro.test(x)$p.value,
 EPPS = function(x) nortsTest::epps.test(x)$p.value  # Ensure we handle p.value for Epps test correctly
 )
 
 # Initialize an empty data frame to store all results
 final_result <- data.frame()
 
 ### Analysis for the Actual Time Series
 actual_cycle <- y  
 
 # Calculate p-values for the actual time series
 p_values <- sapply(all_tests, function(test) {
 result <- test(actual_cycle)
 
 # Check if the result is a p-value or a list, and handle accordingly
 if (is.list(result) && "p.value" %in% names(result)) {
 return(result$p.value)
 } else if (is.list(result) && "epps" %in% names(result)) {
 return(result$epps)  # Extract the relevant part of the Epps test result
 } else {
 return(result)  # For atomic results (like numeric values)
 }
 })
 
 # Skewness and kurtosis for the actual time series
 actual_skewness <- JarqueBera.test(actual_cycle)[[2]]$statistic
 actual_kurtosis <- JarqueBera.test(actual_cycle)[[3]]$statistic - 3  # Excess kurtosis adjustment
 
 # Add actual time series results as the first row and round accordingly
 actual_result <- c(
 n = "Actual Time Series",
 round(p_values, 3),  # P-values for the actual time series rounded to 3 decimal places
 Skewness = round(actual_skewness, 3),
 Kurtosis = round(actual_kurtosis, 3)
 )
 
 # Append the actual time series result to the final result data frame first
 final_result <- rbind(final_result, actual_result)
 
 ### Simulation results for different sample sizes
 # Define the different sample sizes to iterate through
 n_values <- c(71, 142, 213, 355, 710)
 
 # Loop through each sample size in n_values for the simulations
 for (n in n_values) {
 # Extract the corresponding simulation cycle data for current n
 sim_cycle <- lapply(temp_list, function(x) x[1:n, 'cycle'])
 
 # Calculate percentages for stationarity and normality tests
 percentages <- lapply(seq_along(all_tests), function(i) {
 if (names(all_tests)[i] %in% c("PP", "ADF")) {
 # For PP and ADF, calculate percentage where p < 0.05
 return(100 - (mean(sapply(sim_cycle, function(x) all_tests[[i]](x) > 0.05)) * 100))
 } else {
 # For other tests, calculate percentage where p > 0.05
 return(mean(sapply(sim_cycle, function(x) all_tests[[i]](x) > 0.05)) * 100)
 }
 })
 
 # Skewness and kurtosis calculation
 skewness_values <- sapply(sim_cycle, function(x) JarqueBera.test(x)[[2]]$statistic)
 kurtosis_values <- sapply(sim_cycle, function(x) JarqueBera.test(x)[[3]]$statistic)
 
 # Compute average skewness and excess kurtosis
 avg_skewness <- mean(skewness_values)
 avg_kurtosis <- mean(kurtosis_values) - 3  # Excess kurtosis adjustment
 
 # Round the results and ensure percentages are formatted
 formatted_percentages <- sprintf("%.2f", round(unlist(percentages), 2))
 
 # Round skewness and kurtosis to 3 decimal places
 sim_result <- c(
 n = n,
 formatted_percentages,  # Formatted percentages
 Skewness = round(avg_skewness, 3),
 Kurtosis = round(avg_kurtosis, 3)
 )
 
 # Append the simulation result to the final result data frame
 final_result <- rbind(final_result, sim_result)
 }
 
 # Assign meaningful column names
 colnames(final_result) <- c("n", "ADF", "PP", "KPSS", "AD", "SW", "EPPS", "Skewness", "Kurtosis")
 
 file_name <- file.path(results_path, paste0("MCsimulations_", CN, "_", SN, ".xlsx"))
 write.xlsx(final_result, file = file_name, rowNames = TRUE, encoding = "UTF-8")
 
 
 
 
 
 
 
 
