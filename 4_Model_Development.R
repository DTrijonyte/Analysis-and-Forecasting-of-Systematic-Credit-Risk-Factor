
 # Description:   - Identification of the model for a selected country
 
 # The general modelling idea:
 #  1. In pre-Covid sample 2006Q4 (or earlier for EE) -- 2019Q4, pre-screen features
 #  2. Select the model using Regression with ARMA errors fits
 #  3. Using dummy variables for Covid-19 and nDoD, refit the model
 #      or rotate the curves and refit so that covid prediction errors would be mimicked
 #  4. Use dummies to correct the levels of Z and produce Z-adjusted, refit the model
 #  5. Check statistical validity of model and cross-validate accuracy 
 
 # Scenarios currently expect the Z correlation with WAG, YER, HPI is positive
 # and with HIC, URX, E6M(?) is negative
 
 rm(list = ls())
 gc()
 t4 <- Sys.time()
 options(max.print = 200, warn = -1) 
 #to reduce the amount of ggplot warnings warn = -1, use warn = 0 to see all warnings
 
 path <- dirname(rstudioapi::getSourceEditorContext()$path)
 setwd(path)
 
 # Load required libraries, install missing ones:
 source("@CCFunctions.R")
 
 # Take the data from the most recent folder for results:
 CN <- "LT"
 SN <- "PI"
 PL <- 6L  # lower limit for smoothing the noise 
 YQ <- "2024 Q4"#format(as.yearqtr(Sys.Date()), "%Y Q%q")
 lag_max <- 5L # 1 quarter after shock for 90 dpd, up to 1 year delay due to FLI (1-5)
 # Load prepared data set with features and dependent variables
 results_path <- paste0("../4_Results/",YQ)
 output_file_name <- paste("CC_macro_cycles_", CN,".Rdata", sep = "")
 setwd(results_path)
 load(output_file_name)
 setwd(path)
 # Split variables into groups:
 # G1 is the group for which the projections exist at the moment of model's development
 # G1 is the primary group, try the best pre-Covid fit using only G1, and challenge with G2-G3
 G1 <- c("E6M", "URX", "HIC", "HPI", "WAR", "WAG", "YER", "d4HIC", "d4HPI", "d4WAG", "d4WAR")
 # G3 is soft information that reflects economic agents expectations (could be good to fix Covid data)
 G3 <- c("CCI", "ESI", "FCI", "ICI", "RCI", "SCI")
 # G2 additional information economists can potentially predict, good challenger features
 G2 <- setdiff(colnames(all_scaled_data), c(G1, G3, "Z_PI", "Z_LE"))
 
 # Produce single factor analysis table and explore it:
 if(SN == "PI"){eda_output_pi <- tcd_eda(all_scaled_data, segment = "PI", lag_max = lag_max, mincycle = 8L)
 ztable(eda_output_pi)}
 
 
 if(SN == "LE"){eda_output_le <- tcd_eda(all_scaled_data, segment = "LE", lag_max = lag_max, mincycle = 8L)
 ztable(eda_output_le)}
 # LE conclusions: WAR better but not much than WAG, d4WAR is more reasonable
 #                 E6M impact small and in both cases counter-intuitive -- exclude
 
 G1_pi <- c("URX", "HIC", "d4WAG", "d4HPI", "YER") 
 G1_le <- c("URX", "d4WAG", "HIC", "YER")  # exclude d4HPI from LE impacts, E6M has positive correlation
 
 # Expert judgement working on groups relevant for PI and LE
 PI_relevant <- c("E6M", "d4HPI", "URX", "d4WAG", "FCI", "ESI",
 "CCI", "YER", "MTR", "PCR", "JVN")
 LE_relevant <- c("HIC", "URX", "ITR", "OMX", "JVN", "YER", "XTR", "ESI", "d4WAG")
 
 # Basic model includes only G1 variables
 # Challenger model includes also G2 variables
 # G3 information for short run predictions and COVID period with GG support 2020Q1-2021Q2/Q3
 # ------------------------------------------------------
 #  1 step. Pre-screen features in pre-Covid period
 # ------------------------------------------------------
 # training sample end depends on Covid fix starting point, assume fix starts from Q1 2020
 raw_features <- if(SN == "PI") G1_pi else G1_le
 prescreened_features <- prescreen_features(all_scaled_data, SN, raw_features, lag_max, c(2019, 4)) 
 if("prescreened_features" %in% search())detach(prescreened_features)
 attach(prescreened_features)
 
 # --------------------------------------------------------------------
 # 2 step. Select statistically valid regression with ARMA errors
 # --------------------------------------------------------------------
 final_model_cor    <- run_arima_with_lasso(zx_data, correlation, 1)
 final_model_bf     <- run_arima_with_lasso(zx_data, block_forward, 1)
 final_model_all_in <- run_arima_with_lasso(zx_data, all_in, 2, 6)
 
 summary(final_model_cor) 
 summary(final_model_bf) 
 summary(final_model_all_in)
 
 
 cbind("Regression Errors" = residuals(final_model_cor, type="regression"),
 "ARIMA errors" = residuals(final_model_cor, type="innovation")) %>%
 autoplot(facets=TRUE)
 
 checkresiduals(final_model_cor)
 
 # Prediction
 H = 28
 forecast_oos(final_model_cor, H)                          
 forecast_oos(final_model_bf, H)                                 
 forecast_oos(final_model_all_in, H)                          
 
 # --------------------------------------------------------------------
 # 3 step. Refit including data after COVID-19 
 # --------------------------------------------------------------------
 # convert to mts
 #temp_names <- colnames(combined_data)
 #combined_data <- cbind(dates, combined_data[,-1])
 #colnames(combined_data) <- temp_names
 
 start_yearqtr <- paste(z_start_date, collapse = " Q")
 end_yearqrt   <- paste(c(2024, 1), collapse = " Q")  # end of available training data
 combined_data_subset <- subset(combined_data[, -1], as.yearqtr(dates) >= start_yearqtr &
 as.yearqtr(dates) <= end_yearqrt)
 
 z_df <- data.frame(z = window(z, end = c(2024, 1)))
 
 #selected_features <- intersect(names(coef(final_model_all_in)), all_in)
 selected_features <- intersect(names(coef(final_model_bf)), block_forward) 
 #selected_features <- intersect(names(coef(final_model_cor)), correlation)
 #selected_features <- intersect(names(coef(final_model_bs)), best_subset)
 #selected_features <- ps_bf
 #selected_features <- ps_cor
 
 lagged_data      <- data.frame(lapply(as.data.frame(optimistic), create_lagged_vars, lag_max = lag_max))
 optimistic_lags  <- cbind.data.frame(dates = combined_data$dates, optimistic, lagged_data)
 lagged_data      <- data.frame(lapply(as.data.frame(pessimistic), create_lagged_vars, lag_max = lag_max))
 pessimistic_lags <- cbind.data.frame(dates = combined_data$dates, pessimistic, lagged_data)
 

 relaxed <- find_model(setdiff(selected_features, "YER.lag4"), "relaxed", h = 15, order = NULL, 
 print_plots = FALSE)
 pca     <- find_model(ps_cor, "pca", h = 15, order = c(1, 0, 5), print_plots = TRUE, 
 eci_start = c(2007, 3), eci_end = c(2022, 1))
 #pls     <- find_model(ps_cor, "pls", h = 15, order = NULL, print_plots = FALSE, 
 #                      eci_start = c(2007, 2), eci_end = c(2019, 2))
 # LT_PI 2007Q3 2024Q2, LV_PI 2007Q2 2024Q1 
 
 # for PCA method (main chosen method) make cross-validation
 # cv_outputs <- cv_pca(model = pca$model, error_order = c(1, 0, 5), segment = SN, H = 15,
 #                     print_plots = TRUE, cv_start = 2016, cv_end = 2023, lag_max = lag_max, 
 #                     eci_start = c(2007, 3), eci_end = c(2024, 2))
 # cv_relaxed_out <- cv_relaxed(model = relaxed$model, error_order = c(1, 0, 2), segment = SN, H = 15,
 #                          print_plots = TRUE, cv_start = 2022, cv_end = 2023, lag_max = lag_max)
 
 
 # Choose the model to save
 cc_model <- pca
 objects_to_save <- c("cc_model") #, "cv_outputs")
 setwd(results_path) 
 output_file_name <- paste("CC_model_", CN, "_", SN,".Rdata", sep = "")
 save(list = objects_to_save, 
 file = output_file_name)
 setwd(path)
 
 
 cat("Single country model specification completed...\n")
 round(Sys.time() - t4, 1) %>% print()
