
 # Description:   Summarizes cycle peaks and troughs and prepares data for the full cycle definition.
 
 # Attributes:
 #  data            -- time series data matrix
 #  scaling_periods -- data frame with start and end points of full cycles
 
 # Returns data with properly scaled data over full cycle
 function(data, 
 scaling_periods) {
 scaled_center <- scaled_scale <- c()
 
 for (i in 1:nrow(scaling_periods)) {
 variable   <- scaling_periods$Variable[i]
 start_time <- as.numeric(scaling_periods$Start[i])
 end_time   <- as.numeric(scaling_periods$End[i])
 
 ts_data <- na.trim(data[, variable], sides = "both")
 ts_fullcycle <- window(ts_data, start = start_time, end = end_time)
 m <- mean(ts_fullcycle)
 s <- sd(ts_fullcycle)
 scaled_center[variable] <- m
 scaled_scale[variable] <- s
 time_indices <- which(time(data) %in% time(ts_data))
 data[time_indices, variable] <- ts_scaler(ts_data, m, s)
 }
 
 return(list(data = data, scaled_center = scaled_center, scaled_scale = scaled_scale))
 }
