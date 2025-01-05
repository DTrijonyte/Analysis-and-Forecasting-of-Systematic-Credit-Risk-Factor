
 # Description:   Summarizes cycle peaks and troughs and prepares data for the full cycle definition.
 
 # Attributes:
 #  data     -- time series data matrix
 #  max_date -- cut point, since some data is forecasted (not realized) cycle
 
 # Returns list of BBQ analysis outputs for full cycle:
 #  interval_type  -- which full cycle was longer peak to peak or trough to trough
 #  all_data       -- includes min max of sample, p2p, t2t, and selected longest full cycle period
 #  scaling_period -- just selected longest full cycle period
 #  bbq_result     -- full output from BBQ approach if additional analysis would be needed
 
 function(data, max_date = 2024){
 results <- list()
 
 for (i in 1:ncol(data)) {
 ts_series <- na.trim(data[, i], sides = "both")
 ts_series <- window(ts_series, end = max_date)
 ts_times  <- time(ts_series)
 ts_range  <- range(ts_times)
 bbq_result  <- BBQ(ts_series, 8, 4)
 # Use plot for expert judgement, override manually
 bbq_result %>% 
 BCDating::plot(ts_series, main = colnames(data)[i], col = col.custom[1], las = 1, lwd = 2) %>% 
 abline(h = 0, lty = 2)
 
 peaks_range <- ts_times[range(bbq_result@peaks)]
 troughs_range <- ts_times[range(bbq_result@troughs)]
 peaks_diff    <- peaks_range[2] - peaks_range[1] 
 troughs_diff  <- troughs_range[2] - troughs_range[1]
 chosen_range  <- if(peaks_diff > troughs_diff || is.na(troughs_diff)) peaks_range else troughs_range
 
 all_data <- c(ts_range, peaks_range, troughs_range, chosen_range)
 results[[colnames(data)[i]]] <- list(
 interval_type = if (peaks_diff > troughs_diff || is.na(troughs_diff)) "Peak to Peak" else "Trough to Trough",
 all_data = all_data,
 scaling_period = chosen_range,
 bbq_result = bbq_result
 )
 }
 
 return(results)
 }
