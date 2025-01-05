
 # Description:   - Cross-validation of the model, relaxed approach
 
 # Residual analysis and cross validation of prediction accuracies at Z1Y Z2Y Z3Y
 
 # Rolling origin approach:
 # - for every window starting from 2016 repeat the ECI procedure:
 # - pre_screen lags via correlation
 # - estimate weights and parameters for eci_model  (estimate initial model)
 # - predict h = 4, 6, 12 and store the accuracy measures
 # - show their distribution and averages
 # - repeat until the end of sample (operational outputs are till T-4 quarter)
 # - draw the multiple prediction paths for baseline / optimistic and pessimistic (separately)
 
 # Attributes:
 #    error_order  -- ARMA(p, q) order of regression residuals c(p, 0, q)
 #    segment      -- either PI or LE
 #    H            -- forecast horizon
 #    print_plots     -- print plots flag
 #    cv_start     -- cv start year
 #    cv_end       -- cv end year
 #    lag_max      -- maximum lag for feature precreening (up to 5L is reasonable)
 #    ps_method    -- lag prescreening approach (cor, bf, allin)
 
 function(model = relaxed$model,     
 error_order = c(1,0,5), 
 segment = "PI", 
 H = 15L,
 print_plots = TRUE,
 cv_start = 2016,
 cv_end   = 2023,
 lag_max  = 3L,
 ps_method = "cor") {
 model$xreg %>% 
 colnames %>% 
 strsplit("\\.") %>% 
 lapply(function(x)x[1]) %>% 
 unlist() -> feature_names
 segment <- toupper(segment)
 if(!segment %in% c("PI", "LE")) stop(paste("Unknown segment -", segment,
 "- expected one of \"PI\", \"LE\"."))
 
 coef_all <- se_all <- LB_all <- fit_macro_all <- fit_total_all <-              ###### eci_all, eci_weights_all REMOVED
 b_all <- u_all <- l_all <- accuracy_all <- accuracy_z1y <- accuracy_z2y <- 
 accuracy_z3y <- selected_features_all <- NULL  
 stability_set <- list()
 
 df <- length(coef(model))
 
 cv_ends <- cbind(year = rep(cv_start:cv_end, each = 4), quarter = rep(1:4, times = cv_end - cv_start + 1))
 
 G1_segment <- if (segment == "PI") G1_pi else G1_le
 
 for(index in 1:nrow(cv_ends)){
 # print(cv_ends[index,])
 prescreened_features <- prescreen_features(all_scaled_data, segment, c(G1_segment), lag_max, cv_ends[index,])         
 
 # Limitation: lag structure may change but not main feature composition 
 # -- the sign assessment is harder to automate
 
 ps_approach <- switch (ps_method,
 cor   = prescreened_features$ps_cor,
 bf    = prescreened_features$ps_bf,
 allin = G1_segment
 )
 
 selected_features <- NULL
 for(feature in feature_names)selected_features <- c(selected_features, 
 ps_approach[grepl(paste0("^", feature), ps_approach)])
 
 selected_features_all <- rbind(selected_features_all, selected_features)
 start_yearqtr <- paste(z_start_date, collapse = " Q")
 end_yearqrt   <- paste(cv_ends[index,], collapse = " Q") # end of available training data
 combined_data_subset <- subset(combined_data[, -1], as.yearqtr(dates) >= start_yearqtr & as.yearqtr(dates) <= end_yearqrt)
 data <- combined_data_subset[, selected_features, drop = FALSE]    
 
 # eliminate variable with not logical direction of impact
 z_df <- data.frame(z = window(z, end = cv_ends[index,]))
 zx_data <- cbind.data.frame(z = z_df, combined_data_subset)
 zx_data_all <- cbind.data.frame(z = z_df, data)
 
 features_data <- cbind.data.frame(dates = dates[dates >= dates[as.numeric(rownames(data)[1])]], 
 combined_data[dates >=  dates[as.numeric(rownames(data)[1])], 
 selected_features, drop = FALSE])
 stability_model <- run_arima_with_lasso(zx_data, ps_approach)
 stability_set[[index]] <- coef(stability_model)
 updated_model <- Arima(zx_data_all$z, order = error_order, include.mean = FALSE, xreg = as.matrix(data))
 coef_all <- rbind(coef_all, coef(updated_model))
 se_all <- rbind(se_all, sqrt(diag(updated_model$var.coef))) # for statistical validity
 LB_all <- c(LB_all, Box.test(zoo::na.approx(resid(updated_model)), fitdf = df, lag = df + 3, type = "Ljung-Box")$p.value)
 
 fit_macro <- z - residuals(updated_model, type="regression")  # shows the macro fitted part
 fit_total <- z - residuals(updated_model, type="innovation")
 
 fit_macro_all <- ts.union(fit_macro_all, fit_macro)
 fit_total_all <- ts.union(fit_total_all, fit_total)
 
 b <- forecast_oos(updated_model, H, z_main = z, features_data = features_data, do_plot = F)
 
 hist_b <- c(b$x, b$mean) %>% ts(start = tsp(z)[1], frequency = tsp(z)[3])
 hist_u <- c(b$x, b$upper[, 1]) %>% ts(start = tsp(z)[1], frequency = tsp(z)[3])
 hist_l <- c(b$x, b$lower[, 1]) %>% ts(start = tsp(z)[1], frequency = tsp(z)[3])
 b_all <- ts.union(b_all, hist_b)
 u_all <- ts.union(u_all, hist_u)
 l_all <- ts.union(l_all, hist_l)
 f <- b$mean
 a <- window(z, start = tsp(f)[1], end = min(tsp(z)[2], tsp(f)[2]))
 accuracy_all <- rbind(accuracy_all, accuracy(f, a)[c("RMSE", "MAE")])
 #move by one lag due to data lag in Z: Q2 based on Q1 point, Q4 on Q3
 accuracy_z1y <- rbind(accuracy_z1y, c(f[5], a[5]))
 accuracy_z2y <- rbind(accuracy_z2y, c(f[9], a[9]))
 accuracy_z3y <- rbind(accuracy_z3y, c(f[13], a[13]))
 }
 accuracy_z1y <- accuracy(accuracy_z1y[,1],accuracy_z1y[,2])[c("RMSE", "MAE")]
 accuracy_z2y <- accuracy(accuracy_z2y[,1],accuracy_z2y[,2])[c("RMSE", "MAE")]
 accuracy_z3y <- accuracy(accuracy_z3y[,1],accuracy_z3y[,2])[c("RMSE", "MAE")]
 
 colnames(fit_macro_all) <- colnames(fit_total_all) <- colnames(b_all) <- 
 rownames(selected_features_all) <- names(stability_set) <-
 colnames(u_all) <- colnames(l_all) <- apply(cv_ends, 1, paste0, collapse = "Q")
 fit_inertia_all <- fit_total_all - fit_macro_all
 errors_all <- z - fit_total_all
 
 plot_list <- list()
 # plot contribution factors for the last index only == length(cv_ends)
 bar_data <- ts.union(macro = fit_macro_all[,index], inertia = fit_inertia_all[,index], 
 errors = errors_all[,index])
 z_data <- data.frame(Time = time(z)[1:nrow(bar_data)], z = z[1:nrow(bar_data)])
 bar_data_pos_neg <- as.data.frame(bar_data) 
 bar_data_pos_neg$Time <- time(z)[1:nrow(bar_data_pos_neg)]
 bar_data_long <- pivot_longer(bar_data_pos_neg, cols = -Time, 
 names_to = "Series", values_to = "Value")
 
 
 contrib1 <- ggplot() +
 geom_bar(data = bar_data_long, aes(x = Time, y = Value, fill = Series), 
 stat = "identity", color = "black") +
 geom_line(data = z_data, aes(x = Time, y = z), color = "black", size = 1.2) +
 scale_fill_manual(values = col.custom, 
 labels = c("Residual", "Inertia", "Macroeconomic data")) +
 scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +  
 labs(x = "", y = "") +
 theme(legend.position = "top") +
 geom_hline(yintercept = 0, color = "black", size = 0.5)
 plot_list[["contrib1"]] <- contrib1
 
 # bar plot with each variable contribution
 macro_contributions <- coef_all[index, (sum(error_order) + 1L):ncol(coef_all)]
 
 total_data <- t(t(data) * macro_contributions)
 bar_data <- ts.union(errors_all[,index], fit_inertia_all[,index], total_data)
 colnames(bar_data) <- c("errors", "inertia", colnames(total_data))
 # if at least one value NA set all values to NA
 bar_data %>% apply(1, function(x)sum(is.na(x))) %>% {which(.>0)} -> na_values
 if(length(na_values) >0) bar_data[na_values,] <- NA
 
 bar_data %>% as.data.frame() -> bar_data_frame
 bar_data_frame$Time <- time(z)[1:nrow(bar_data_frame)]
 
 bar_data_long <- pivot_longer(bar_data_frame, cols = -Time, names_to = "Series", 
 values_to = "Value")
 
 contrib2 <- ggplot() +
 geom_bar(data = bar_data_long, aes(x = Time, y = Value, fill = fct_inorder(Series)), 
 stat = "identity", color = "black") +
 geom_line(data = z_data, aes(x = Time, y = z), color = "black", size = 1.2) +
 scale_fill_manual(values = col.lg.rev(ncol(bar_data))) +
 scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
 labs(x = "", y = "") +
 theme(legend.position = "top") +
 geom_hline(yintercept = 0, color = "black", size = 0.5)
 plot_list[["contrib2"]] <- contrib2
 
 # Quality of total fits including inertia
 total_fits <- autoplot(fit_total_all, xlab = "", ylab = "", color = col.custom[1], alpha = 0.2, ylim = c(-3, 3)) + 
 autolayer(z, lwd = 1, color = "black") +
 geom_hline(yintercept = 0, color = "black", size = 0.5)
 plot_list[["fits"]] <- total_fits
 
 # Forecasts with upper and lower bonds should not deviate from red-blue interval
 uncertainty <- autoplot(b_all, xlab = "", ylab = "", color = col.custom[1], alpha = 0.5, ylim = c(-3, 3)) + 
 autolayer(u_all, color = col.custom[2], alpha = 0.5) +
 autolayer(l_all, color = col.custom[3], alpha = 0.5) +
 autolayer(z, lwd = 1, color = "black") +
 geom_hline(yintercept = 0, color = "black", size = 0.5)
 plot_list[["uncertainty"]] <- uncertainty
 
 # Coefficients checks
 coef_all %>% ts() %>% autoplot() -> coef_plot
 plot_list[["coef"]] <- coef_plot
 
 if(print_plots)lapply(plot_list, print)
 
 cat("\n\nCheck percent insignificant t statistics:\n")
 t_all <- coef_all/se_all
 print(sum(abs(t_all) < 2, na.rm = T) / length(t_all) * 100)
 cat("\n\nCheck percent insignificant t statistics just for macro:\n")
 t_macro <- t_all[,(sum(error_order) + 1L):ncol(coef_all)]
 print(sum(abs(t_macro) < 2, na.rm = T) / length(t_macro) * 100)  
 
 print(summary(LB_all))
 
 cat("\n\n Accuracy of predictions over all quarters:\n")
 print(summary(accuracy_all))
 cat("\n\n Accuracy of predictions for Z1y, Z2y, and Z3y:\n")
 print(accuracy_z1y)
 print(accuracy_z2y)
 print(accuracy_z3y)
 
 return(list(feature_names  = feature_names,
 coef = coef_all,
 LB   = LB_all,
 fit_macro = fit_macro_all,
 fit_total = fit_total_all,
 b = b_all,
 l = l_all,
 u = u_all,
 accuracy_all = accuracy_all,
 accuracy_z1y = accuracy_z1y,
 accuracy_z2y = accuracy_z2y,
 accuracy_z3y = accuracy_z3y,
 selected_features = selected_features_all,
 stability_set = stability_set,
 plots = plot_list))
 
 }
