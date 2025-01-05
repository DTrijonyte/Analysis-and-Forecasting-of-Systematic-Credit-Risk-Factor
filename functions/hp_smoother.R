
 # Description:   Hodrick-Prescott smoothing function of time series y

 # Arguments:
 #   y        - time series data
 #   wl       - 2/T, where denominator T shows around how many time points to optimize
 
 function(y, wl = 2/PL){
 error <- optimize(function(x){hp.q(y, x, w1 = wl, nfreq = 1e4, order = NULL)}, 
 c(0,100))$minimum
 return(hp.filter(y, error)[,"trend"])
 }
