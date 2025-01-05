
 # Description:   Generates lagged versions of a given column in a data frame
 
 # Arguments:
 #     column - the input column from a data frame for which lagged variables will be created
 #     lag_min - the starting point for lagging; by default, it is set to 1
 #     lag_max - the ending point for lagging
 
 function(column, lag_min = 1, lag_max) {
 lagged_columns <- lapply(lag_min:lag_max, function(lag) lag(column, lag))
 return(setNames(data.frame(lagged_columns), paste0(names(column), "lag", lag_min:lag_max)))
 }
