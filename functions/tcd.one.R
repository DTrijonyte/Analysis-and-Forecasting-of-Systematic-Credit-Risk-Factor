
 # Description:   Wrapping function for running single variable TCS decomposition. Light version.

 # Arguments:
 #   var_nm         - variable name
 #   upper_limit    - the maximum search limit in quarters
 #   lower_limit    - the minimum search limit in quarters
 #   show_plots     - should the plots be printed (recommended for single-run analysis)
 #   use_all_tcs    - use all 3 TCS variants with orders of stochastic cycle from 1 to 3
 #   no_cycle_thr   - test if the not scaled cycle significantly deviates from 0 in amplitude
 
 function(var_nm         = "",
 upper_limit    = PU,    # 80 is reasonable for business cycles in Baltics
 lower_limit    = PL,    # 6 classical choice, 8 often alternative from annual data 2 years 
 show_plots     = FALSE,
 use_all_tcs    = FALSE,
 no_cycle_thr   = 1.5
 ){
 # 0) Extract processed data from the first step
 y  <- data[, var_nm]
 y  <- na.omit(y) 
 # In I(0) case we allow drift to be time-varying, i.e. the series are always drift adjusted
 ND <- max(1, ndiffs(y))
 # Allow for a single structural break in drift (maybe also sacrificed)
 breakpoint <- attributes(ur.za(y))$bpoint 
 
 cat("Working on", var_nm, "\n")
 N     <- length(y)
 START <- tsp(y)[1]
 FREQ  <- frequency(y)
 if(!is.null(breakpoint)) DU <- matrix(c(rep(0, breakpoint - 1), 
 rep(1, N - breakpoint + 1)), 
 ncol = 1) else DU <- NULL
 
 # 1) Asymmetric Christiano-Fitzgerald filter (least dependent on assumptions)
 # as in Grinderslev et al. assume I(1) with drift and MA(q) errors see CF paper (1999/2002)
 mod_cf <- auto.arima(y, max.p = 0, seasonal = FALSE, d = 1,
 max.q = 12, max.order = 12, stationary = FALSE)
 theta_length <- mod_cf$arma[2]
 root_cf <- 1 #(ND > 0)
 thetas  <- 1
 if(theta_length > 0) thetas <- c(thetas, mod_cf$coef[1:theta_length])
 # upper and lower limits are also from the paper but 6 for lower is used in short cycles
 # full cycle due to fast decay of waves ~ short + long cycle
 bpf <- cffilter(y, pl = lower_limit, pu = upper_limit, drift = TRUE, root = root_cf, theta = thetas) 
 # short-run  cycle 1.5-10 years (Regulatory documents often are about BC cycle lengths, 3+ cycles per sample)
 bps <- cffilter(y, pl = lower_limit, pu = 40, drift = TRUE, root = root_cf, theta = thetas) 
 # long-run cycle 11-30 years (~1 full cycle in the sample, harder to detect since the spectral density is narrow)
 bpl <- cffilter(y, pl = 40, pu = upper_limit, drift = TRUE, root = root_cf, theta = thetas)
 # compare the amplitudes of cycles
 # use bpf for peaks-troughs analysis
 cycle_f <- ts(bpf$cycle, start = START, freq = FREQ)
 cycle_s <- ts(bps$cycle, start = START, freq = FREQ)
 cycle_l <- ts(bpl$cycle, start = START, freq = FREQ)
 cycles  <- ts.union(cycle_f, cycle_s, cycle_l)
 if(show_plots){
 autoplot(cycles*100, lwd = 1) +
 labs(y = "", x = "") +
 theme_custom() +
 scale_color_manual(values = col.custom,
 labels = c("Full", "Short", "Long")) + 
 theme(legend.position = "right") -> picture
 print(picture)
 }
 # Relative importance (similar to Hyndman and Athanasopoulos (2018)):
 rel_importance_s <- max(0, (var(cycle_s) + cov(cycle_s, cycle_l))/var(cycle_f))
 rel_importance_l <- max(0, (var(cycle_l) + cov(cycle_s, cycle_l))/var(cycle_f))
 # Should sum up to 1
 cat("Relative importances, short:", round(rel_importance_s * 100, 1), 
 "% and long:", round(rel_importance_l * 100, 1), "%\n")
 
 # 3 years auto spectrum makes the series sharp enough around the mode point of the band
 max_cycle_freq <- which.max(spec.ar(cycle_f, 1e5, order = 12, plot = FALSE)$spec) 
 
 if(max_cycle_freq != 1){
 max_cycle_length <- 1e5*pi/max_cycle_freq # CF expected cycle upper limit
 } else { 
 max_cycle_length <- upper_limit # restrict to upper_limit/4 years max
 }
 
 if(FREQ == 4)cat("Max cycle length is", round(max_cycle_length, 1), 
 "quarters or", round(max_cycle_length/4, 1), "years\n")
 # this is a new upper cut-off frequency for CF (it is data driven yet depends on initial max cut-off!)
 bpu   <- cffilter(y, pl = lower_limit, pu = max_cycle_length, drift = TRUE, 
 root = root_cf, theta = thetas)
 cycle <- ts(bpu$cycle, start = START, freq = FREQ) # main CF output
 mid_range <- (max_cycle_length - lower_limit)/2 + lower_limit # mid-range point of the interval [6, max_cycle_length]
 1:3 %>% lapply(function(x){stochastic.cycle.nc(cycle, cyc  = mid_range, nc = x,
 m0 = 0, rho = 0.975, C0 = 1e-4)}) -> ini_list
 # Drop with errors
 ini_list %>% lapply(function(x)x$rho) %>% unlist() %>% is.na() %>% which() -> which_na
 if(length(which_na) > 0L) ini_list <- ini_list[-which_na]
 
 # 2) HP filter (in simulated examples works close to band-pass filter [6, mean_cycle_length])
 hp_lambda <- optimize(function(x){hp.q(y, x, w1 = 2/max_cycle_length, nfreq = 1e5, order = NULL)}, c(0,500000))$minimum
 hp_cycle  <- hp.filter(y, hp_lambda)[,"cycle"]
 # remove residual spikes upto 6 quarters and take trend as smoothed cycle:
 hp_error  <- optimize(function(x){hp.q(hp_cycle, x, w1 = 2/lower_limit, nfreq = 1e4, order = NULL)}, c(0,500000))$minimum
 hp_cycle_smoothed  <- hp.filter(hp_cycle, hp_error)[,"trend"]
 
 
 # 3) TCS filter by Mohr
 N_list <- length(ini_list)
 ini_list %>% lapply(function(x){
 tcs(y, nd = ND, nc = x$nc, cyc = x$cyc, rho = x$rho, du = DU, s = 1)$data[,"cycle"]}) %>% 
 unlist() %>% 
 matrix(ncol = N_list) %>% 
 ts(start = START, frequency = FREQ) -> tcs_cycles
 
 # Pick the one which is closest to CF cycle in correlation distance (default) terms 
 # or min MSE/MAPE and (possible) use all to catch the broader uncertainty ranges
 # initial values for UCM
 choice <- which.max(cor(cbind(cycle, tcs_cycles))[1, 2:(N_list + 1)])
 ini <- ini_list[[choice]]
 cyc <- ini$cyc
 rho <- ini$rho
 nc  <- ini$nc
 if(!use_all_tcs){
 tcs_cycle <- tcs_cycles[, choice]
 cycles <- ts.union(tcs_cycle*100,
 hp_cycle_smoothed*100,
 cycle*100)
 tcs_names <- "TCS"
 }else{
 cycles <- ts.union(tcs_cycles*100,
 hp_cycle_smoothed*100,
 cycle*100)
 tcs_names <- paste("TCS", 1:3, sep = "")
 }
 
 
 # combine 1)-3)
 
 # 4) Suite of 4 model combined into final cycle
 #     mean(range(x, na.rm = TRUE)) -- mid range (useful if particular class of 
 #                                                models dominates, e.g. TCS)
 #     mean(x, na.rm = TRUE)        -- arithmetic mean
 
 cycles %<>% {ts.union(., apply(., 1, function(x)mean(range(x, na.rm = TRUE))))}
 colnames(cycles) <- c(tcs_names, "HP", "CF", "Suite")
 
 if(show_plots){
 suite_range <- ts.union(center = cycles[, ncol(cycles)],
 min = apply(cycles, 1, min, na.rm = TRUE),
 max = apply(cycles, 1, max, na.rm = TRUE))
 x_coords <- xy.coords(suite_range[, 1])$x 
 
 polygon_coords <- data.frame(x = c(x_coords, rev(x_coords)),
 y = c(suite_range[, 2], rev(suite_range[, 3])))
 ggplot() +
 geom_polygon(aes(x = x, y = y), data = polygon_coords, fill = "#6eaad2", alpha = 0.3) +  
 autolayer(cycles, lwd = 1) +
 ggtitle("Suite of models, %") + 
 ylab("") + xlab("") + 
 theme(legend.position = "right") + 
 scale_color_manual(values = c(col.custom[1:2], 1, col.custom[3])) -> picture
 print(picture)
 }
 #Correlation between the cycles, suspicious flag if UCM is different from the rest of the methods.
 manual_check <- min(cor(cycles)) < 0.5
 cycle_check  <- max(abs(range(cycles[, ncol(cycles)]))) > no_cycle_thr
 #Trends:
 trends <- y - cycles / 100
 colnames(trends) <- colnames(cycles)
 
 #Store the objects potentially useful for outputs
 out <- list(cycles    = cycles,
 trends    = trends,
 manual    = manual_check,
 sig_cycle = cycle_check,
 rel_s     = rel_importance_s,
 rel_l     = rel_importance_l,
 thetas    = thetas,
 max_cycle_length = max_cycle_length,
 hp_lambda = hp_lambda,
 hp_error  = hp_error,
 ini       = ini)
 
 invisible(out)
 }
