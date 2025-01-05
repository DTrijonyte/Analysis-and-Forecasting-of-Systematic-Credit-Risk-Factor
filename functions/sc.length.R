
 # Description:   Length of stochastic cycle

 # Arguments
 #  y    - time series
 #  cyc0 - intitial cycle length
 #  rho0 - initial dampening parameter
 #  breakpoint - break point
 #  nd   - level of the trend (usually in the range of 0 - 2)
 #  nc   - level of the 1st cyclical process (usually in the range of 0 - 2)
 #  s    - the season: 4 for quarterly data, 12 for monthly data, etc., s = 1 (default), no season
 #  tol  - tolerance for iterations
 #  max.it - max number of iterations 
 
 function(y, 
 cyc0   = 20, 
 rho0   = 0.975, 
 breakpoint = NULL, 
 nd     = 1, 
 nc     = 1, 
 s      = 1, 
 tol    = 0.001,
 max.it = 100){
 # Mohr (2005) iterations for c = 1 only...
 # 1. Compute the cycle with initial values for mu and rho
 N    <- length(y)
 
 if(is.null(breakpoint)){
 DU <- NULL
 }else{
 DU   <- matrix(c(rep(0, breakpoint - 1), rep(1, N - breakpoint + 1)), ncol = 1)
 }
 tcs_cycle <- tcs(y, nd = nd, nc = nc, cyc = cyc0, rho = rho0, du = DU, s = s) %>% {.$data[,"cycle"]} 
 #browser()
 
 # 2. Estimate $\alpha(L)cycle_t = \beta(L)\nu_t$ for c = 1 only
 # tricky here, first iterate with restriction on \theta_1 = -\phi_1/2
 i <- 0
 rhon <- rho0
 cycn <- cyc0
 rho0 <- 0  # will be above tol for initial iteration
 cyc0 <- 120
 while((rhon - rho0)^2 + (2*pi/cycn - 2*pi/cyc0)^2 > tol & i < max.it){
 s_cycle <- stochastic.cycle(tcs_cycle, tol = tol, max.it = max.it, cyc = cycn, rho = rhon)
 cyc0 <- cycn
 rho0 <- rhon
 cycn <- s_cycle$cyc
 rhon <- s_cycle$rho
 tcs_cycle <- tcs(y, nd = nd, nc = nc, cyc = cycn, rho = rhon, du = DU, s = s) %>% {.$data[,"cycle"]}
 i <- i + 1
 }
 return(list(i = i, cyc = cycn, rho = rhon, cycle = tcs_cycle))
 }
