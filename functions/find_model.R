
 # Description:   Estimates the regression model over pre-screened features
 
 #   selected_features -- which columns we will choose
 #   method            -- c("relaxed" (the first option), "pca", "pls")
 #   h                 -- prediction horizon h = 11
 #   print_plots       -- TRUE if to print plots
 #   order             -- NULL if automatic ARMA order is good, c(p, 0, q) if overridden
 #   eci_start         -- start date for ECI full cycle period
 #   eci_end           -- end date for ECI full cycle period
 
 # Returns list with model outputs, prints diagnostics and plots 
 
 function(selected_features, 
 method = "relaxed", 
 h = 11,
 print_plots = TRUE,
 order = NULL,
 eci_start = c(2007, 3),
 eci_end   = c(2024, 2)) {
 if(!method %in% c("relaxed", "pca", "pls")) {
 stop(paste("Unknown method -", method,
 "- expected one of \"relaxed\", \"pca\", \"pls\"."))}
 # Add projections for selected model from different scenarios
 data  <- combined_data_subset[, selected_features]
 plots <- list()
 
 if(method == "pca") {
 data %>% na.omit() %>% prcomp() -> pr_out
 eci_weights <- pr_out$rotation[, 1] 
 } 
 if(method == "pls") {
 pls_data <- cbind.data.frame(z = z_df, data)
 pls_fit <- plsr(z ~., data = pls_data, cnomp = 1)
 eci_weights <- pls_fit$loadings[,1] 
 }
 
 if(method %in% c("pca", "pls")){
 starting_nas <- apply(data, 2, function(x)(length(x) - length(na.trim(x,"left")))) %>% max() 
 # Create an Economic Cycle Indicator
 eci_raw <- as.matrix(combined_data[dates >= dates[as.numeric(rownames(data)[1])], 
 names(eci_weights)]) %*% eci_weights %>% 
 na.omit() %>% {ts(c(rep(NA, starting_nas), .), 
 start = tsp(zx_data[,1])[1], frequency = tsp(zx_data[,1])[3])}
 if(sign(cor(na.omit(ts.union(eci_raw, z)))[1, 2]) == -1){
 eci_raw <- -1 * eci_raw  # ECI must correlate with z positively
 eci_weights <- -1 * eci_weights
 }
 m <- mean(window(eci_raw, eci_start, eci_end))
 s <- sd(window(eci_raw, eci_start, eci_end))
 eci <- ts_scaler(eci_raw, m, s)
 bbq_output <- BBQ(na.omit(eci), 8, 4)
 cat("Peaks and troughs periods:\n")
 print(bbq_output)
 
 plots[["bbq"]] <- function(){
 bbq_output %>% 
 BCDating::plot(na.omit(eci), las = 1, lwd = 2, col = col.custom[1]) %>% 
 abline(h = 0, lty = 2)
 }
 
 # Optimistic and pessimistic work only for G1 group
 if(names(eci_weights) %>% 
 lapply(function(x)any(grepl(substr(x, 1, 3), G1))) %>% 
 unlist() %>% all()){
 eci_optimistic <- as.matrix(optimistic_lags[dates >= dates[as.numeric(rownames(data)[1])], 
 names(eci_weights)]) %*% eci_weights %>% na.omit() %>% c()
 eci_pessimistic <- as.matrix(pessimistic_lags[dates >= dates[as.numeric(rownames(data)[1])], 
 names(eci_weights)]) %*% eci_weights %>% na.omit() %>% c()
 eci_o <- ts(c(rep(NA, starting_nas), ts_scaler(eci_optimistic, m, s)), 
 start = tsp(zx_data[,1])[1], frequency = tsp(zx_data[,1])[3])
 eci_p <- ts(c(rep(NA, starting_nas), ts_scaler(eci_pessimistic, m, s)), 
 start = tsp(zx_data[,1])[1], frequency = tsp(zx_data[,1])[3])
 pessimistic_lags <- data.frame(dates =  dates[dates >= tsp(eci_p)[1]], eci = eci_p)
 optimistic_lags  <- data.frame(dates =  dates[dates >= tsp(eci_o)[1]], eci = eci_o)
 } else {eci_o <- eci_p <- NULL}
 # Select regression with ARMA model
 zx_data_all <- cbind.data.frame(z = z_df, eci = eci[1:nrow(z_df)])
 selected_features <- "eci"
 features_data <- data.frame(dates =  dates[dates >= tsp(eci)[1]], eci = eci)
 }else{
 features_data <- cbind.data.frame(dates = dates[dates >= dates[as.numeric(rownames(data)[1])]], 
 combined_data[dates >=  dates[as.numeric(rownames(data)[1])], 
 selected_features])
 zx_data_all <- cbind.data.frame(z = z_df, data)
 eci_weights <- NULL
 eci <- eci_o <- eci_p <- NULL
 }
 
 # The same for every method
 cat('\nExclude statistically insignificant variables, the last variable is a stop point:\n')
 updated_model <- run_arima_with_lasso(zx_data_all, selected_features)
 # If auto run choice of ARMA errors was insufficient for checkresiduals, add more MA terms
 if(!is.null(order)){
 updated_model <- Arima(zx_data_all$z, order = order, include.mean = FALSE, 
 xreg = updated_model$xreg)
 }
 cat('\nSummary of the final model:\n')
 summary(updated_model)
 cat('\nP-values of t-tests:\n')
 df <- updated_model$nobs - length(coef(updated_model))
 
 2 * pt(abs(coef(updated_model)/sqrt(diag(updated_model$var.coef))), df, 
 lower.tail = FALSE) %>% round(3) %>% print()
 
 norder <- sum(updated_model$arma[1:4])
 # Check if all signs are as expected (initially vars are in the search space if correlation is correct):
 expected_signs <- cor(zx_data_all, use = "pairwise.complete.obs")[1, -1, drop = F]
 modeled_signs <- coef(updated_model)[(1 + norder):length(coef(updated_model))]
 model_var_selection <- names(modeled_signs)
 wrong_signs <- sign(expected_signs[, model_var_selection]) != sign(modeled_signs[model_var_selection])
 if(any(wrong_signs))cat("Signs for variables:", names(modeled_signs)[which(wrong_signs)], 
 "are not as expected, consider removal!\n" )
 
 # Correlation plot
 zx_data_all[,c("z", intersect(names(coef(updated_model)), selected_features))] %>%
 as.matrix() %>% data.frame() %>%
 GGally::ggpairs() +
 theme(axis.text.x = element_blank(),
 axis.text.y = element_blank())   -> plots[["corr"]]
 
 
 # Unit roots check:
 autoplot(updated_model) + theme(aspect.ratio = 1) -> plots[["ur"]]
 
 # Add the plotting function to the list
 plots[["resid"]] <- function() {checkresiduals(updated_model)}
 plots[["foos"]] <- function() {
 forecast_oos(updated_model, h, features_data = features_data, z_main = z)
 }
 
 checkresiduals(updated_model, plot = FALSE)
 
 print(shapiro.test(resid(updated_model)))
 
 # different scenarios predictions:
 b <- forecast_oos(updated_model, h, z_main = z, features_data = features_data, do_plot = F)
 
 z0 <- c(b$x, b$mean) %>% ts(start = tsp(z)[1], frequency = tsp(z)[3])
 if(!is.null(eci_o)|(method == "relaxed")){
 p <- forecast_oos(updated_model, h, z_main = z, features_data = pessimistic_lags, do_plot = F)
 o <- forecast_oos(updated_model, h, z_main = z, features_data = optimistic_lags, do_plot = F)
 z_out <- ts.union(z_baseline = z0, z_optimistic = o$mean, z_pessimistic = p$mean)
 } else {z_out <- z0}
 
 
 forec <- autoplot(z0) +
 ylim(-3, 3) + 
 labs(y = "", x = "") +
 geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
 if(!is.null(eci_o)|(method == "relaxed")){
 forec <- forec + 
 autolayer(o$mean, color = "red") +
 autolayer(p$mean, color = "blue") 
 }
 plots[["forec"]] <- forec
 
 if(print_plots){
 if(!is.null(plots[['bbq']]))plots[['bbq']]()
 print(plots[['corr']])
 print(plots[['ur']]) 
 plots[['resid']]()
 plots[['foos']]()
 print(plots[['forec']])
 }
 return(list(model = updated_model, b = b, o = o, p = p, z_out = z_out, method = method, 
 weights = eci_weights, eci = eci,  plots = plots, zx_data_all = zx_data_all,
 eci_start = eci_start, eci_end = eci_end))
 }
