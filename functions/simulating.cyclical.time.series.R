
 # Description:   Simulates time series data from the Unobserved components model ECB version. The simulated model follows the model specified in 
 # Grinderslev et al. 2017 box 6.
 
 # Arguments:
 #   n             - length of time series to be generated
 #   lambda        - lambda governs the cycle length 2*pi/lambda is the cycle length so a lambda of 0.2 results in a cycle 
 #                   length of 31 quarters corresponding to approx 7.5-8 years
 #   rho           - rho governs the dampening parameter in the cycle state equation, a high rho yields persistent cycles.
 #   rhoc          - secondary dampening parameter, see Grinderslev et al 2017 Box 6
 #   alfa          - alfa is the loading parameter governing the impact of the cycle on the series.
 #   a             - a parameter we don't estimate but which could introduce long waves in the levels if a~[0.9, 1]
 #   eps_kappa_sd  - the standard deviation of the kappa error term, i.e the cycles error see Grinderslev et al 2017 box 6
 #   eps_zeta_sd   - the standard deviation of the zeta error term, i.e. the slope/trend error see Grinderslev et al 2017 box 6
 #   eps_eta_sd    - the standard deviation of the eta error term, i.e. the level error see Grinderslev et al 2017 box 6
 #   eps_irr_sd    - the standard deviation of the eps error term, i.e. the observed component error 
 #                   see Grinderslev et al 2017 box 6
 #   start_time    - starting year of time series
 #   freq          - frequency of time seires
 #   m0            - intital state of cycles (vector with 4 components)
 #   level         - initial state of trend's level
 #   slope         - initial state of trend's slope
 #   breakpoint    - when the structural break occurs
 #   break_level   - breaks are permanent level shifts times break_level
 #   break_drift   - break in drift instead of level's shift occurs
 #   seed          - random seed
 #   plot          - logical if plots should be produced
 #   main          - title's part used in plots
 
 #example: simulating.cyclical.time.series(n=200,lambda=0.2, rho=0.9, rhoc=0.9, alfa=0.5, eps_kappa_sd=0.1, eps_zeta_sd=0.1,
 #eps_eta_sd=0.1, eps_irr_sd=0.1, start_time=1990, freq=4)
 
 #Returns the simulated series and the cyclical component of the cycle
 
 
 function(n       = 67,
 lambda  = 0.137357598         , 
 rho     = 0.760121746    ,
 rhoc    = 0.645080307           ,
 alfa    = 0.356089465        , # sd for the cycle is actually driven by alfa
 a       = 1,    
 eps_kappa_sd = 1, # untouchable, alpha controls for that
 eps_irr_sd_2 = 0.028604708,   # noise
 start_time = c(2006, 4),
 freq       = 4,
 nc         = 1,
 m0         = rep(1, 2 + 2*nc),
 seed       = NULL,
 plot       = FALSE,
 main       = "Simulated series"){
 eps_irr_sd <- sqrt(eps_irr_sd_2)
 #Generating each random series add + 1 to the seed, it insures the series could be prolonged increasing n
 if(!is.null(seed))set.seed(seed)
 # Generate irreducible error
 error <- ts(rnorm(n, sd = eps_irr_sd), start = start_time, frequency = freq)
 # Generate cycle
 if(!is.null(seed))set.seed(seed + 1)
 kappa      <- rnorm(n, sd = eps_kappa_sd)
 if(!is.null(seed))set.seed(seed + 2)
 kappa_star <- rnorm(n, sd = eps_kappa_sd)
 
 # low level cycles:
 T_matrix <- rho*matrix(c(cos(lambda), sin(lambda),
 -sin(lambda), cos(lambda)), 
 nrow = 2, byrow = TRUE)
 lm <- cbind(rbind(matrix(0, 2, 2*(nc - 1)), diag(2*nc - 2)), matrix(0, 2*nc, 2))
 GG <- kronecker(diag(nc), T_matrix) + lm
 
 omegas <- cbind(m0[3:length(m0)], matrix(0, nrow = 2*nc, ncol = n))
 kappas <- rbind(kappa, kappa_star, matrix(0, nrow = 2*(nc - 1), ncol = n))
 
 for (t in 2:(n + 1)) { 
 omegas[, t] <- GG %*% omegas[, t - 1, drop = F] + kappas[, t - 1, drop = F]
 }
 omega <- omegas[(nc - 1)*2 + 1,]
 omega_star <- omegas[(nc - 1)*2 + 2,]
 
 # High level cycles
 C1 <- c(m0[1], rep(0,n))
 C2 <- c(m0[2], rep(0,n))
 
 for (t in 2:(n + 1)) { 
 
 C1[t] <- omega[t - 1] + rhoc * C1[t - 1]
 C2[t] <- omega_star[t - 1] + rhoc * C2[t - 1]
 
 }
 # Only C1 participates in cycle creation
 cycle <- ts(alfa * C1[-1], start = start_time, frequency = freq)
 series <- cycle + error  
 
 if(plot){
 plot(series,las = 1,col = "#63003D", lwd = 2)
 lines(cycle, col = "#6EAAD2", lwd = 2)
 legend("topleft", legend = c("Actual","Cycle"), border = FALSE,
 bty = "n", col = c("#63003D", "#6EAAD2"), lwd = 2)
 }
 # Return the data
 data <- ts.union(series, cycle, error)
 options(warn = 0)
 invisible(data)
 }
