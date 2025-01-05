
 # Description:   Trims leading and trailing NA values from a time series or vector and then applies Kalman smoothing to the remaining data if there are any internal NAs

 # Arguments:
 #              x - vector or time series that may contain NA values
 
 function(x) {
 x_trimmed <- na.trim(x, "both")
 gap_na <- sum(is.na(x_trimmed))
 
 # If there are NAs, apply Kalman smoothing
 if (gap_na > 0L) {
 x_trimmed <- na_kalman(x_trimmed)
 }
 
 return(x_trimmed)
 }
