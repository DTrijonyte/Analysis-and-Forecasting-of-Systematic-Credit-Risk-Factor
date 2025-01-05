
 # Description:   Wrapping function for running single TCS decomposition for simulated data
 #                From Celov and Comunale (2021)
 
 # Arguments:
 #   y    - simulated time series
 #   y_c  - simulated stochastic cycle or any cycle to compare against
 #   cn   - geographic location
 #   ND   - prior information on the number of differences
 #   breakpoint - prior information on the position of the structural break
 #   og   - make output gap transformation
 
 # By default any time series process is decomposed into the sum of trend + cycle + residual
 # The output gap is defined as actual ((data - residual)/trend - 1) so the transformation is (exp(cycle) - 1)*100  
 
 function(y, y_c, cn = "", ND = 1L, breakpoint = NULL, og = TRUE){
 N     <- length(y)
 START <- tsp(y)[1]
 FREQ  <- frequency(y)
 if(!is.null(breakpoint)) DU <- matrix(c(rep(0, breakpoint - 1), 
 rep(1, N - breakpoint + 1)), 
 ncol = 1) else DU <- NULL
 # ----------------------------------------------------------------------------------------------
 # A. Polynomial trends
 du <- c(DU)
 poly_res <- poly.trend(y, degree = 1L, du = du, do_plot = FALSE)
 POLY  <- poly_res[,"cycle"]
 POLY1 <- hp_smoother(POLY)*100
 poly_res <- poly.trend(y, degree = 2L, du = du, do_plot = FALSE)
 POLY  <- poly_res[,"cycle"]
 POLY2 <- hp_smoother(POLY)*100
 poly_res <- poly.trend(y, degree = 3L, du = du, do_plot = FALSE)
 POLY  <- poly_res[,"cycle"]
 POLY3 <- hp_smoother(POLY)*100
 # ----------------------------------------------------------------------------------------------
 # B. First differencees (acceleration cycles)
 SFD  <- 100*diff(y, lag = 1)  %>% hp_smoother
 AFD  <- 100*diff(y, lag = 4)  %>% hp_smoother
 LFD  <- 100*diff(y, lag = 16) %>% hp_smoother
 # Example to demonstrate how wrong could the growth rates be
 # Looses first lag = . of observations
 # ----------------------------------------------------------------------------------------------
 # C. One-sided filters
 # Beveridge-Nelson filter
 
 bnf_dm <-  try({bnf(y = y*100, d = 1, p = 12, auto_delta = T, d0 = 0.01, demean = "dm", wind = 40)
 })
 if(class(bnf_dm) == "try-error")BNF  <- NULL else BNF <- bnf_dm$cycle %>% hp_smoother
 # Hamilton (never-HP) filter
 y_n <- as.xts(y*100)
 dimnames(y_n) <- list(NULL, "GDP")
 HPnever   <- yth_filter(y_n, h = 8, p = 4, output = "cycle") %>% c() %>%
 ts(start = START + 2.75, frequency = FREQ) %>% hp_smoother
 HPnever_a <- yth_filter(y_n, h = 16, p = 4, output = "cycle") %>% c() %>%
 ts(start = START + 4.75, frequency = FREQ) %>% hp_smoother
 # ----------------------------------------------------------------------------------------------
 # D. Asymmetric two-sided filters
 # Christiano-Fitzgerald
 
 mod_cf <- auto.arima(y, max.p = 0, seasonal = FALSE,
 max.q = 12, max.order = 12, stationary = FALSE)
 theta_length <- mod_cf$arma[2]
 root_cf <- 1 #(ND > 0)
 thetas  <- 1
 if(theta_length > 0) thetas <- c(thetas, mod_cf$coef[1:theta_length])
 # upper and lower limits are also from the paper but 6 for lower is used in short cycles
 # upper threshold could be changes
 # full cycle due to fast decay of waves ~ short + long cycle
 bpf <- cffilter(y, pl = PL, pu = PU, drift = TRUE, root = root_cf, theta = thetas) 
 # short-run  cycle 1.5-11 years (Regulatory documents often are about BC cycle lengths, 3+ cycles per sample)
 bps <- cffilter(y, pl = PL, pu = 44, drift = TRUE, root = root_cf, theta = thetas) 
 # long-run cycle 11-30 years (~1 full cycle in the sample, harder to detect since the spectral density is narrow)
 bpl <- cffilter(y, pl = 44, pu = PU, drift = TRUE, root = root_cf, theta = thetas)
 # compare the amplitudes of cycles
 # use bpf for peaks-troughs analysis
 CFf <- ts(bpf$cycle, start = START, freq = FREQ)*100
 CFs <- ts(bps$cycle, start = START, freq = FREQ)*100
 CFl <- ts(bpl$cycle, start = START, freq = FREQ)*100
 
 # Relative importance (similar to Hyndman and Athanasopoulos (2018)):
 rel_importance_s <- max(0, (var(CFs) + cov(CFs, CFl))/var(CFf))
 rel_importance_l <- max(0, (var(CFl) + cov(CFs, CFl))/var(CFf))
 
 # Heuristic rule to switch to drift = TRUE option, for longer UCM cycles
 drift <- rel_importance_l > 2/3
 
 # 5 years auto spectrum makes the series sharp enough around the mode point of the band
 max_cycle_freq_s <- which.max(spec.ar(CFs, 1e5, order = 20, plot = FALSE)$spec) 
 max_cycle_freq_l <- which.max(spec.ar(CFl, 1e5, order = 20, plot = FALSE)$spec)
 # the weighted average is the expected cycle length
 max_cycle_freq   <- rel_importance_s * max_cycle_freq_s + rel_importance_l * max_cycle_freq_l
 if(max_cycle_freq != 1){
 max_cycle_length <- 1e5*pi/max_cycle_freq # CF expected cycle upper limit
 } else { 
 max_cycle_length <- upper_limit # restrict to upper_limit/4 years max
 }
 
 # this is a new upper cut-off frequency for CF (it is data driven yet depends on initial max cut-off!)
 bpu   <- cffilter(y, pl = PL, pu = max_cycle_length, drift = TRUE, 
 root = root_cf, theta = thetas)
 CFo <- ts(bpu$cycle, start = START, freq = FREQ)*100 # main CF output
 mid_range <- (max_cycle_length - PL)/2 + PL # mid-range point of the interval [6, max_cycle_length]
 1:3 %>% lapply(function(x){stochastic.cycle(CFo, cyc  = mid_range, nc = x,
 m0 = 0, rho = 0.975, C0 = 1e-4, 
 fixed_cyc = FALSE, silent = TRUE)}) -> ini_list
 # Drop with errors
 ini_list %>% lapply(function(x)x$rho) %>% unlist() %>% is.na() %>% which() -> which_na
 if(length(which_na) > 0L) ini_list <- ini_list[-which_na]
 
 
 mid_range_f <- (120 - PL)/2 + PL # mid-range point of the interval [6, max_cycle_length]
 1:3 %>% lapply(function(x){stochastic.cycle(CFf, cyc  = mid_range, nc = x,
 m0 = 0, rho = 0.975, C0 = 1e-4, 
 fixed_cyc = FALSE, silent = TRUE)}) -> ini_list_f
 # Drop with errors
 ini_list_f %>% lapply(function(x)x$rho) %>% unlist() %>% is.na() %>% which() -> which_na
 if(length(which_na) > 0L) ini_list_f <- ini_list_f[-which_na]
 
 # Hodrick-Prescott filters
 # Standard penalty
 HPs <- hp_smoother(hp.filter(y, 1600)[, "cycle"])
 # Longer penalty
 HPl <- hp_smoother(hp.filter(y, 51200)[, "cycle"])
 # Optimized (in simulated examples works close to band-pass filter [8, mean_cycle_length])
 hp_lambda <- optimize(function(x){hp.q(y, x, w1 = 2/max_cycle_length, nfreq = 1e5, order = NULL)}, c(0,500000))$minimum
 HPo <- hp_smoother(hp.filter(y, hp_lambda)[, "cycle"])*100
 
 #boosted HP
 bx_ADF <- try({BoostedHP(c(y), lambda = hp_lambda, iter= TRUE, stopping = "adf")})
 if(class(bx_ADF) == "try-error") HPob <- NULL else{
 HPb  <- bx_ADF$cycle %>% c()
 HPob <- hp_smoother(ts(HPb, start = START, frequency = FREQ))*100}
 
 hp_lambda <- optimize(function(x){hp.q(y, x, w1 = 2/120, nfreq = 1e5, order = NULL)}, c(0,500000))$minimum
 HPof <- hp_smoother(hp.filter(y, hp_lambda)[, "cycle"])*100
 
 # Butterworth filter
 bwf_model <- tcs(y, nd = 3, nc = 0, cyc = hp_lambda)
 BWF <- bwf_model$data[, "cycle"] %>% hp_smoother  %>% {.*100}
 
 # Trigonometric filter
 #tr_model <- trfilter(y[-1], pl = 6, pu = max_cycle_length, drift = TRUE)
 #TRIGO <- tr_model$cycle*100
 
 # Waveletes
 # drift adjustment
 drift <- (y[N] - y[1])/(N-1)
 y_drift_adj <- y - (0:(N-1))*drift
 
 wave_mod <- try({mra(y_drift_adj, J = 7, boundary = "reflection")})
 if(class(wave_mod) == "try-error") WAVE <- NULL else {
 # D1 and D2 are ommitted 2 and 4 quarter cycles
 # D3 + D4 + D5 + D6 are 4-8-16-32-64-128 cycles
 WAVE <- (wave_mod$D3 + wave_mod$D4 + wave_mod$D5 + wave_mod$D6+ wave_mod$D7)  %>% 
 ts(start = START, frequency = FREQ)*100}
 
 # ----------------------------------------------------------------------------------------------
 # E. Structural models
 # TCS filter by Mohr
 N_list <- length(ini_list)
 ini_list %>% lapply(function(x){
 tcs(y, nd = ND, nc = x$nc, cyc = x$cyc, rho = x$rho, du = DU, s = 1)$data[,"cycle"]}) %>% 
 unlist() %>% 
 matrix(ncol = N_list) %>% 
 ts(start = START, frequency = FREQ) -> tcs_cycles
 
 # Pick the one which is closest to CF cycle in correlation distance terms or min MSE/MAPE
 # or (possible) use all to catch the broader uncertainty ranges
 # initial values for UCM
 choice <- which.max(cor(cbind(CFo, tcs_cycles))[1, 2:(N_list + 1)])
 ini <- ini_list[[choice]]
 cyc <- ini$cyc
 rho <- ini$rho
 nc  <- ini$nc
 
 TCS <- tcs_cycles[,choice]*100
 # UCM model
 UCM <- uc.ecb(y, method = "MLE", fixed = FALSE, cyc = cyc, rho = rho, 
 rhoc = 0.5, m0 = rep(0, 4), print_out = FALSE, plot = FALSE)$data[,3]*100
 
 N_list <- length(ini_list_f)
 ini_list_f %>% lapply(function(x){
 tcs(y, nd = ND, nc = x$nc, cyc = x$cyc, rho = x$rho, du = DU, s = 1)$data[,"cycle"]}) %>% 
 unlist() %>% 
 matrix(ncol = N_list) %>% 
 ts(start = START, frequency = FREQ) -> tcs_cycles
 
 # Pick the one which is closest to CF cycle in correlation distance terms or min MSE/MAPE
 # or (possible) use all to catch the broader uncertainty ranges
 # initial values for UCM
 choice <- which.max(cor(cbind(CFf, tcs_cycles))[1, 2:(N_list + 1)])
 
 ini <- try(ini_list_f[[choice]])
 if(class(ini) == "try-error"){TCSf <- TCS}else{
 cyc <- ini$cyc
 rho <- ini$rho
 nc  <- ini$nc
 
 TCSf <- tcs_cycles[,choice]*100}
 
 # UCM model
 #UCMf <- uc.ecb(y, method = "MLE", fixed = FALSE, cyc = cyc, rho = rho, 
 #              rhoc = 0.5, m0 = rep(0, 4), print_out = FALSE, plot = FALSE)$data[,3]*100
 
 
 
 
 # ----------------------------------------------------------------------------------------------
 # F. Suite of models
 # Short list
 SUITEs  <- ts.union(CFo, HPo, TCS) %>% 
 {ts.union(., SoM = apply(., 1, function(x)mean(range(x, na.rm = TRUE))))} %>% {.[,"SoM"]}
 SUITEsa <- ts.union(CFf, HPo, TCS) %>%
 {ts.union(., SoM = apply(., 1, function(x)mean(range(x, na.rm = TRUE))))} %>% {.[,"SoM"]}
 # Long list (no acceleration cycle and CF subcomponents)
 SUITEf <- ts.union( POLY1, POLY2, POLY3,
 BNF, HPnever, HPnever_a,
 CFf, CFo,
 HPo, HPob,
 BWF, WAVE, #TRIGO
 TCS, UCM
 ) %>%
 {ts.union(., SoM = apply(., 1, function(x)mean(range(x, na.rm = TRUE))))} %>% {.[,"SoM"]}
 SUITEsf <- ts.union(CFf, HPof, TCSf) %>%
 {ts.union(., SoM = apply(., 1, function(x)mean(range(x, na.rm = TRUE))))} %>% {.[,"SoM"]}
 SUITEsfa <- ts.union(CFf, HPof, TCSf, POLY2) %>%
 {ts.union(., SoM = apply(., 1, function(x)mean(range(x, na.rm = TRUE))))} %>% {.[,"SoM"]}
 if(is.null(y_c))y_c <- SUITEsf #Default choice for comparison
 cycles_all <- ts.union(y_c,
 POLY1, POLY2, POLY3,
 SFD, AFD, LFD,
 BNF, HPnever, HPnever_a,
 CFf, CFo, CFs, CFl, 
 HPo, HPob, HPof, HPs, HPl,
 BWF, WAVE, #TRIGO
 TCS, TCSf, UCM, #UCMf,
 SUITEs, SUITEsa, SUITEf, SUITEsf, SUITEsfa
 )
 dates <- as.yearqtr(time(y_c))
 if(og) cycles_all %<>% {exp(./100)*100 - 100} 
 rownames(cycles_all) <- dates
 plot(cycles_all, plot.type = "single", col = alpha("black", 0.4), las = 1, ylab = "", xlab = "")
 abline(h = 0, lwd = 2, lty = 2)
 title(cn)
 lines(y_c, col = "#2958a3", lwd = 3)
 y_c_bbq <- BBQ(y_c, mincycle = 8, minphase = 4) 
 f=file()
 sink(file = f)
 
 N_cases <- ncol(cycles_all)
 
 cycles_all %>% data.frame() %>% 
 lapply(function(x){x <- ts(x, start = START, frequency = FREQ)
 try(BBQ(na.omit(x), mincycle = 8, minphase = 4)) -> y
 if(class(y) == "try-error")return(NULL)
 overlap <- ts.union(attributes(y_c_bbq)$states, attributes(y)$states) %>%
 apply(1, prod) %>% {. == 1} %>% na.omit() %>% sum()
 c(BCDating::summary(y)) -> z
 
 z <- try(c(z, z[3]+z[4], overlap/length(x)))
 return(z)}) -> bbq_res
 sink()
 close(f)
 w_null <- lapply(bbq_res, function(x)!is.null(x)) %>% unlist()
 no_null <- sum(!w_null)
 bbq_res %<>% unlist %>% matrix(nrow = N_cases - no_null, byrow = TRUE) 
 eda <- cbind(cor(cycles_all[,w_null], use = "pairwise.complete.obs")[,1], bbq_res,
 c(0, dist(bbq_res)[1:(N_cases - 1 - no_null)]),
 apply(cycles_all[,w_null], 2, function(x)sum(is.na(x))))
 colnames(eda) <- c("Correlaiton","Amp.Exp", "Amp.Rec", "Dur.Exp", "Dur.Rec", "Cycle.length", 
 "Overlap", "Dissim.","NA")
 
 res <- list(cycles = cycles_all,
 eda   = eda)
 return(res)
 }
