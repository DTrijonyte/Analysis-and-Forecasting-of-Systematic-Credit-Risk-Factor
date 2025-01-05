
 # Description:   Unobserved components model ECB version MLE/Bayes

 
 # Arguments:
 #   y         - time series
 #   one.sided - real time or ex post solution
 #   plot      - plot the figures
 #   cyc       - average cycle length (prior)
 #   rho       - dampening parameter for the stochastic AR-cycle model (based on Mohr (2005))
 #   rhoc      - dampening parameter for the long cycle
 #   fixed     - are the cycle parameters estimated, fixed = TRUE means the cyc and rhos are assumed
 #   bp        - breakpoint in trend for MLE only
 #   method    - use MLE or Bayes
 #   a.y       - prior mean of signal equation's variance
 #   b.y       - prior variance of signal equation's variance
 #   a.theta   - prior mean of state equations' variances
 #   b.theta   - prior variance of state equations' variances
 #   n.sample  - the number of MCMC samples
 #   thin      - number of iterations to be discarded -- number of iterations (1+thin)*n.sample
 #   burn      - the size of burn-in sample
 #   ma        - should the MA part be included (see Mohr (2004))
 #   print_out - print the diagnostic output
 #   seed      - set random seed
 #   m0_cycle  - prior beliefs on the credit cycle starting point (default is 0)
 
 function(y, 
 one.sided = FALSE, 
 plot      = TRUE,
 main      = "UCM ECB",
 cyc       = 40,
 rho       = 0.9,
 rhoc      = 0.5,
 fixed     = TRUE,
 drift     = FALSE,
 bp        = NULL,
 method    = c("MLE", "Bayes"),
 a.y       = 0.01, #prior mean
 b.y       = 1000, #prior variance
 a.theta   = 0.01,
 b.theta   = 1000,
 n.sample  = 1100,  # increase n-sample latter, just for testing 1100+
 thin      = 1,     # discard thin iterations for every saved iteration
 burn      = 100,
 print_out = TRUE,
 seed      = 20210101,
 m0_cycle  = c(rep(0, 4)),
 drift_res = 1
 ){
 options(warn = -1)
 set.seed(seed)
 if(!method %in% c("MLE", "Bayes")) stop("method should be \"MLE\" or \"Bayes\".")
 if(!"ts" %in% class(y)) stop("series is not a \"ts\" object.")
 N <- length(y)
 PI <- base::pi
 # Set state priors
 level <- m0 <- y[1]
 slope <- mean(diff(y), na.rm=TRUE) 
 # (y[T] - y[1])/(T - 1) true estimate has trend[T] - trend[1]
 m0 <- c(level, slope)
 
 # Set up dynamic linear model in state space form:
 model <- function(X){
 #TREND
 trend <- dlmModPoly(order = 2,
 dV    = exp(X[3]),   # error
 dW    = exp(X[4:5]), # level, slope
 m0    = m0,
 C0    = diag(c(1e7, 1e-4))) # diffuse prior for trend, tight for slope
 if(!is.null(bp)){
 # break-point in the state of slope is done through inflated variance of the level
 trend$JW <- matrix(1, nr = 1, nc = 2) 
 trend$X  <- matrix(exp(X[4:5]), nr = N, nc = 2, byrow = TRUE)
 trend$X[bp,1] <- trend$X[bp,1] * (1 + exp(X[6])) 
 if(drift){trend$X[,2] <- 0}
 }
 if(drift){
 trend$W[2,2] <- 0 # set slope to be a deterministic drift
 }
 #CYCLE
 if(fixed){
 cycle <- dlm(FF = matrix(c(exp(X[8]), 0, 0, 0), nrow = 1), 
 V = 0,
 GG = rbind(cbind(diag(rhoc, 2), diag(1, 2)),
 cbind(diag(0, 2),
 rho*matrix(c(cos(2*PI/cyc), sin(2*PI/cyc),
 -sin(2*PI/cyc), cos(2*PI/cyc)), 
 nrow = 2, byrow = TRUE))
 ),
 W = diag(c(0, 0, 1, 1)), 
 m0 = m0_cycle, #c(0, 0, 0, 0), 
 C0 = 1e-4*diag(4)
 ) 
 }
 else{
 cycle <- dlm(FF = matrix(c(exp(X[8]), 0, 0, 0), nrow = 1), 
 V = 0,
 GG = rbind(cbind(diag(X[7], 2), diag(1, 2)),
 cbind(diag(0, 2),
 X[1]*matrix(c(cos(2*PI/X[2]), sin(2*PI/X[2]),
 -sin(2*PI/X[2]), cos(2*PI/X[2])), 
 nrow = 2, byrow = TRUE))
 ),
 W = diag(c(0, 0, 1, 1)), 
 m0 = m0_cycle, #c(0, 0, 0, 0), 
 C0 = 1e-4*diag(4)
 ) 
 }
 ucm <- trend + cycle
 return(ucm)
 }
 
 if(method == "MLE"){
 # MLE Estimation, lower and upper restrictions are needed (tighter intervals could be more benefitial)
 MLE <- dlmMLE_m(y, parm = c(rho, cyc, -5, -5, -5, -10, rhoc, 0), model, 
 lower = c(0.5, 6,       rep(-50, 4), 0.4, -6),
 upper = c(0.999, 2*cyc, rep(  0, 4), 0.8, 0), method = "NMK"
 )
 #  print(MLE)
 # Estimated parameters
 EstParams <- model(MLE$par)
 #browser()
 if(print_out)print(MLE$par)
 # Smoothed or filtered series
 if(one.sided){
 Smooth_Estimates <- dlmFilter(y, EstParams)
 if(is.null(dim(Smooth_Estimates$m))){
 trend <- dropFirst(Smooth_Estimates$m)
 cycle <- (y - trend)*exp(MLE$par[8])
 }else{
 trend <- dropFirst(Smooth_Estimates$m[, 1])
 cycle <- dropFirst(Smooth_Estimates$m[, 3])*exp(MLE$par[8])
 }
 }else{
 Smooth_Estimates <- dlmSmooth(y, EstParams)
 if(is.null(dim(Smooth_Estimates$s))){
 trend <- dropFirst(Smooth_Estimates$s)
 cycle <- (y - trend)*exp(MLE$par[8])
 }else{
 trend <- dropFirst(Smooth_Estimates$s[, 1])
 cycle <- dropFirst(Smooth_Estimates$s[, 3])*exp(MLE$par[8])
 }
 }
 
 }else{
 # Bayes method starts here, copied from the book
 PARsupport <- function(u){ 
 ## Box type of support for parameters
 (0 < u[1])&&(u[1] < 1) && 
 (PI/60 < u[2])&&(u[2] < 2*PI/3) && 
 (0 < u[3])&&(u[3] < 1) && 
 (0 < u[4])&&(u[4] < 1)
 }
 PARfullCond <- function(u){
 ## log full conditional density for unknown model parameters: c(rho, lambda, ar1, loading)
 mod$GG[3:6, 3:6] <- rbind(cbind(diag(u[3], 2), diag(1, 2)),
 cbind(diag(0, 2),
 u[1]*matrix(c(cos(u[2]), sin(u[2]),
 -sin(u[2]), cos(u[2])), 
 nrow = 2, byrow = TRUE))
 )
 mod$FF[1, 3] <- u[4]
 # Prior information is stored here
 -dlmLL(y, mod) + 
 sum(dnorm(u[1:3], mean = c(rho, 2*PI/cyc, rhoc), sd = rep(0.2, 3), log = TRUE)) +
 dgamma(1/u[4], 5^2/100 , 5/100, log = TRUE) - 2 * log(u[4]) # prior on visibility of the cycle
 }
 # From textbook and dlmGibbsDIG() and gdpGibbs() functions, still could be with errors
 # Initialize the model as in MLE
 mod <- model(c(rho, cyc, -5, -5, -5, -10, rhoc, -1))
 
 every    <- thin + 1
 mcmc     <- n.sample * every
 p        <- ncol(mod$FF) # dim of the state space
 r        <- 2            # number of unknown variances in W
 nobs     <- N
 # signal error term
 shape.y  <- a.y^2 / b.y + 0.5 * nobs
 rate.y   <- a.y   / b.y
 # states errors
 shape.theta <- a.theta^2 / b.theta + 0.5 * nobs
 rate.theta  <- a.theta   / b.theta
 shape.drift <- a.theta^2 / b.theta / drift_res^2 + 0.5 * nobs # smaller prior mean on drift
 rate.drift  <- a.theta   / b.theta / drift_res 
 theta       <- matrix(0, nobs + 1, p)
 gibbsTheta  <- array(0, dim = c(nobs + 1, p, n.sample))
 gibbsV      <- vector("numeric", n.sample) 
 gibbsW      <- matrix(0, nrow = n.sample, ncol = r)
 gibbsPAR    <- matrix(0, nrow = n.sample, ncol = 4)
 it.save     <- 0
 #browser()
 if(print_out) pb <- txtProgressBar(0, mcmc, style = 3)
 for (it in 1:mcmc){ 
 if ( print_out ) setTxtProgressBar(pb, it)
 ## A. Sample the parameters starting from previous value:
 ## generate Cycle matrix rho*H parameters rho and cyc
 a_it   <- mod$GG[5, 5:6]
 rho_it <- sqrt(a_it %*% a_it)
 lambda_it <- acos(a_it[1]/rho_it)  # 2*PI/cyc
 rhoc_it <- mod$GG[3, 3]
 loading_it <- mod$FF[1, 3]
 sampled <- arms(c(rho_it, lambda_it, rhoc_it, loading_it), PARfullCond, PARsupport, 1)
 mod$GG[3:6, 3:6] <- rbind(cbind(diag(sampled[3], 2), diag(1, 2)),
 cbind(diag(0, 2),
 sampled[1]*matrix(c(cos(sampled[2]), sin(sampled[2]),
 -sin(sampled[2]), cos(sampled[2])), 
 nrow = 2, byrow = TRUE))
 )
 mod$FF[1, 3] <- sampled[4]
 
 ## B. generate states - FFBS
 modFilt <- dlmFilter(y, mod, simplify = TRUE)
 theta[] <- dlmBSample(modFilt)
 ## generate V
 
 y.center <- y - tcrossprod(theta[-1, , drop = FALSE], mod$FF)
 SSy <- drop(crossprod(y.center))
 mod$V[] <- 1/rgamma(1, shape = shape.y, rate = rate.y + 0.5 * SSy)
 
 ## generate W
 theta.center <- theta[-1,, drop = FALSE] -
 theta[-(nobs + 1), , drop = FALSE] %*% t(mod$GG)
 SStheta <- drop(sapply(1:p, function(i) crossprod(theta.center[,i])))
 mod$W[1, 1] <-
 1 / rgamma(1, shape = shape.theta,
 rate  = rate.theta + 0.5 * SStheta[1])
 if(drift){mod$W[2, 2] = 0} else {
 mod$W[2, 2] <-
 1 / rgamma(1, shape = shape.drift,
 rate  = rate.drift + 0.5 * SStheta[2])
 }
 ## save current iteration, if appropriate
 if (!(it %% every) )
 {    #browser()
 it.save <- it.save + 1
 #if(print_out & !(it.save %% 100)) print(it.save)
 gibbsTheta[,,it.save] <- theta
 gibbsPAR[it.save,] <- sampled
 gibbsV[it.save] <- diag(mod$V)
 gibbsW[it.save,] <- diag(mod$W)[1:r]
 }
 }
 gibbs_res <- list(dV = gibbsV, dW = gibbsW, theta = gibbsTheta, par = gibbsPAR)
 thetaMean <- ts(dropFirst(apply(gibbs_res$theta[,,-(1:burn)],1:2,mean)),start=1990,freq=4)
 trend <- thetaMean[,1]
 #browser()
 cycle <- ts(apply(t(t(gibbs_res$theta[,3,-(1:burn)])*gibbs_res$par[-(1:burn),4])
 ,1,mean)[-1],start=1990,freq=4) #thetaMean[,3]
 if(print_out){
 cat("\n\n")
 out1 <- gibbs_res$par[-(1:burn),]
 colnames(out1) <- c("rho", "lambda", "ar1", "alpha")
 print(format(mcmcMeans(out1), digits = 1L, scientif = FALSE), quote = FALSE)
 print(format(apply(out1, 2, quantile, probs = c(.05,.95)), digits = 1L, scientif = FALSE), quote = FALSE)
 out2 <- cbind(gibbs_res$dV[-(1:burn)], gibbs_res$dW[-(1:burn),])
 if(drift){out2 <- out2[,-3]
 colnames(out2) <- c("sigma_err", "sigma_level")
 }else{
 colnames(out2) <- c("sigma_err", "sigma_level", "sigma_drift")
 }
 print(format(mcmcMeans(out2), digits = 1L, scientif = FALSE), quote = FALSE)
 print(format(apply(out2, 2, quantile, probs = c(.05,.95)), digits = 1L, scientif = FALSE), quote = FALSE)
 }
 }
 # Error term
 error <- y - trend - cycle
 #browser()
 # Plot the data ---
 if(plot){
 par(mfrow = c(2,1),
 oma = c(1,3,0,0) + 0.1,
 mar = c(1.5,1,1,1) + 0.1)
 plot(y,las = 1,col = "#003755", family = "Danske Text", lwd = 2)
 lines(trend, col = "#6dbace", lwd = 2)
 legend("topleft", legend = c("Observed","Trend"), border = FALSE,
 bty = "n", col = c("#003755", "#6dbace"), lwd = 2)
 title(main = paste(main, " - Trend", sep = ""))
 par(mar = c(1,1,1.5,1) + 0.1)
 plot(cycle*100, las = 1, col = "#e65a6d", lwd = 2)
 title(main = paste(main, " - Cycle", sep = ""))
 abline(h = 0)
 par(mfrow = c(1,1),
 oma = c(0,0,0,0),
 mar = c(5.1,4.1,4.1,2.1))  
 }
 
 # path<-"I:/Dat/38AIALL/Projects/All/Global/Regulatory PD/2.Scripts/1.CC/Development"
 # source(paste0(path,"/functions/test-ucm.ecb.r"),local=TRUE)
 
 # Return the data
 data <- ts.union(y, trend, cycle, error)
 options(warn = 0)
 if(method == "MLE"){
 return(invisible(list(data = data, MLE = MLE)))
 } else{
 if(print_out) close(pb)
 return(invisible(list(data = data, gibbs_res = gibbs_res)))
 }
 
 
 
 }
