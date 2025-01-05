
 # Description:    HP filter: lambda optimization statistic

 # Arguments:
 #   y      - time series object with frequency
 #   lambda - smoothing constant
 #   w1     - cutoff frequency, see notes for more details
 #   nfreq  - number of points in spectral density to search through
 #   order  - the order of AR(p) to estimate for spectral density analysis
 #            set to NULL for AIC based automatic choice
 
 # Notes:
 #   Q is the weighted absolute value of the difference between the ideal
 #   and the distortionary filter
 #   Hhp is the power transfer function of the HP filter
 #   H is the power transfer function of the ideal filter
 #   w1 cutoff frequency for ideal filter, specified according
 #   2/"number of periods to pass through"
 #   to the view of the duration of business cycles 
 
 function(y, lambda, w1 = 1/20, nfreq = 10000, order = 2){
 hp.spec <- spec.ar(y, n.freq = nfreq, order = order, plot = FALSE)
 
 ## this estimates the spectral density
 
 w <- hp.spec$freq * 2/frequency(y) # transform to good frequencies [0,1]*pi
 S <- hp.spec$spec
 deltaw <- (w[2] - w[1])*pi
 Hhp <- abs((4 * lambda * (1 - cos(pi*w))^2) /
 (4 * lambda * (1 - cos(pi*w))^2 + 1))^2
 H <- ifelse(w < w1, 0, 1)
 
 NormSum <- 2 * sum(S) * deltaw
 v <- 2 * S * deltaw / NormSum
 
 Q <- sum(abs(H - Hhp)*v)
 Q
 
 }
