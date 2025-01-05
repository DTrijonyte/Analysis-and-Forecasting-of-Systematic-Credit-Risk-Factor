
 # Description:  Hodrick-Prescott filter

 # Arguments:
 #   x      - time series
 #   lambda - the penalty parameter
 
 
 function(x, lambda = 100){
 eye    <- diag(length(x))
 trend  <- solve(eye + lambda * crossprod(diff(eye, lag = 1, d = 2)), x)
 cycle  <- x - trend
 # Return the data
 data   <- ts.union(x, trend, cycle)
 return(data)
 }
