
 # Description:   Produce baseline predictions.

 # Attributes:
 #   model    -- regression with ARIMA errors output
 #   h        -- prediction horizon
 #   z_main   -- actual Z factor data
 #   features_data -- data set with values for prediction from scenarios
 #   do_plot  -- logical to make a picture or not
 
 function(model, h = 12, z_main = z, features_data = combined_data, do_plot = TRUE){
 selected_variables <- colnames(model$xreg)
 end_date_history <- tsp(model$x)[2]%>%as.yearqtr()
 start_date_prediction <- tsp(model$x)[2] + 0.25
 features_prediction <- subset(features_data[, -1, drop = FALSE], 
 as.yearqtr(features_data[,1]) > end_date_history) %>%
 {as.matrix(.[1:h, selected_variables, drop = FALSE])} %>%
 {ts(., start = start_date_prediction, frequency = 4)}
 forecast_outputs <- forecast(model, h = h, xreg = features_prediction)
 
 z_data <- data.frame(date = time(z_main), z = z_main)
 fit_macro <- model$x - residuals(model, type="regression")  # shows the macro fitted part
 fit_total <- model$x - residuals(model, type="innovation") 
 fits_data <- data.frame(date = time(model$x), macro = fit_macro, total = fit_total)
 
 if(do_plot){
 p <- autoplot(forecast_outputs, ylab = "Z-factor", flwd = 1, fcol = col.custom[2]) +
 autolayer(forecast_outputs$mean, series = "Forecast") +
 geom_line(aes(x = date, y = z, color = "Actual"), data = z_data, lwd = 1)  +
 geom_line(aes(x = date, y = macro, color = "Macro fit"), data = fits_data, lwd = 1, lty = 2) +
 geom_line(aes(x = date, y = total, color = "Total fit"), data = fits_data, lwd = 1, lty = 2)  +
 labs(color = "Legend") + xlab("Year") +
 scale_color_manual(values = c(col.custom[1:2], "blue", col.custom[3]))  +
 guides(col = guide_legend(nrow = 1)) + # To add fixed number of rows or columns if needed
 theme(panel.background = element_rect(fill = NA),
 panel.ontop = TRUE)+
 geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5)
 print(p)
 }
 return(invisible(forecast_outputs))
 }
