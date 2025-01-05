
 # Description:   Unobserved components model ECB version MLE/Bayes
 #                Excluding trend component just cycle

 # Arguments:
 #   y         - time series
 #   one.sided - real time or ex post solution
 #   plot      - plot the figures
 #   cyc       - average cycle length (prior)
 #   rho       - dampening parameter for the stochastic AR-cycle model (based on Mohr (2005))
 #   rhoc      - dampening parameter for the long cycle
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
 main      = "Comparison",
 cyc       = 40,
 rho       = 0.9,
 rhoc      = 0.5,
 fixed     = TRUE,
 method    = c("MLE", "Bayes"),
 a.y       = 0.01, #prior mean
 b.y       = 1000, #prior variance
 a.theta   = 0.01,
 b.theta   = 1000,
 n.sample  = 1100,  # increase n-sample latter, just for testing 1100+
 thin      = 1,     # discard thin iterations for every saved iteration
 burn      = 100,
 print_out = TRUE,
 seed      = 20230930,
 m0_cycle  = c(rep(y[1], 4)),
 drift_res = 1
 ){
 options(warn = -1)
 set.seed(seed)
 require(dlm) # should be already loaded
 if(!method %in% c("MLE", "Bayes")) stop("method should be \"MLE\" or \"Bayes\".")
 if(!"ts" %in% class(y)) stop("series is not a \"ts\" object.")
 N <- length(y)
 PI <- base::pi
 
 # Set up dynamic linear model in state space form:
 model <- function(X){
 #CYCLE
 if(fixed){
 cycle <- dlm(FF = matrix(c(exp(X[5]), 0, 0, 0), nrow = 1), 
 V  = exp(X[3]),   # error
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
 cycle <- dlm(FF = matrix(c(exp(X[5]), 0, 0, 0), nrow = 1), #alpha
 V  = exp(X[3]),                               #signal error
 GG = rbind(cbind(diag(X[4], 2), diag(1, 2)),  #rho_c
 cbind(diag(0, 2),
 X[1]*matrix(c(cos(2*PI/X[2]), sin(2*PI/X[2]),
 -sin(2*PI/X[2]), cos(2*PI/X[2])), 
 nrow = 2, byrow = TRUE))
 ),
 W = diag(c(0, 0, 1, 1)), 
 m0 = m0_cycle, #c(0, 0, 0, 0), 
 C0 = 1e-2*diag(4)
 ) 
 }
 return(cycle)
 }
 
 if(method == "MLE"){
 # MLE Estimation, lower and upper restrictions are needed (tighter intervals could be more beneficial)
 MLE <- dlmMLE_m(y, parm  = c(rho, cyc, -2, rhoc, -1), model, 
 lower = c(0.1,    6, -50, 0.1, -10),
 upper = c(0.999, 48,  -1, 0.9,  10), 
 method = "NMK")
 # Estimated parameters
 EstParams <- model(MLE$par)
 if(print_out)print(MLE$par)
 # Smoothed or filtered series
 if(one.sided){
 Smooth_Estimates <- dlmFilter(y, EstParams)
 cycle <- dropFirst(Smooth_Estimates$m[, 1])*exp(MLE$par[5])
 }else{
 Smooth_Estimates <- dlmSmooth(y, EstParams)
 cycle <- dropFirst(Smooth_Estimates$s[, 1])*exp(MLE$par[5])
 }
 }else{
 # Bayes method
 PARsupport <- function(u){ 
 ## Box type of support for parameters
 (0 < u[1])&&(u[1] < 1) && 
 (PI/60 < u[2])&&(u[2] < 2*PI/3) && 
 (0.3 < u[3])&&(u[3] < 0.7) && 
 (0 < u[4])&&(u[4] < 1)
 }
 PARfullCond <- function(u){
 ## log full conditional density for unknown model parameters: c(rho, lambda, ar1, loading)
 model_init$GG <- rbind(cbind(diag(u[3], 2), diag(1, 2)),
 cbind(diag(0, 2),
 u[1]*matrix(c(cos(u[2]), sin(u[2]),
 -sin(u[2]), cos(u[2])), 
 nrow = 2, byrow = TRUE))
 )
 model_init$FF[1, 1] <- u[4]
 # Prior information is stored here (sd controls tightness)
 -dlmLL(y, model_init) + 
 sum(dnorm(u[1:3], mean = c(rho, 2*PI/cyc, rhoc), sd = rep(0.2, 3), log = TRUE)) +
 dgamma(1/u[4], 5^2/100 , 5/100, log = TRUE) - 2 * log(u[4]) # prior on visibility of the cycle
 }
 # From textbook and dlmGibbsDIG() and gdpGibbs() functions, still could be with errors
 # Initialize the model as in MLE
 model_init <- model(c(rho, cyc, -2, rhoc, -1))
 
 every    <- thin + 1
 mcmc     <- n.sample * every
 p        <- ncol(model_init$FF) # dim of the state space
 r        <- 2            # number of unknown variances in W
 nobs     <- N
 # signal error term
 shape.y  <- a.y^2 / b.y + 0.5 * nobs
 rate.y   <- a.y   / b.y
 theta       <- matrix(0, nobs + 1, p)
 gibbsTheta  <- array(0, dim = c(nobs + 1, p, n.sample))
 gibbsV      <- vector("numeric", n.sample) 
 gibbsPAR    <- matrix(0, nrow = n.sample, ncol = 4)
 it.save     <- 0
 #browser()
 if(print_out) pb <- txtProgressBar(0, mcmc, style = 3)
 for (it in 1:mcmc){ 
 if ( print_out ) setTxtProgressBar(pb, it)
 ## A. Sample the parameters starting from previous value:
 ## generate Cycle matrix rho*H parameters rho and cyc
 a_it   <- model_init$GG[3, 3:4]
 rho_it <- sqrt(a_it %*% a_it)
 lambda_it <- acos(a_it[1]/rho_it)  # 2*PI/cyc
 rhoc_it <- model_init$GG[1, 1]
 loading_it <- model_init$FF[1, 1]
 sampled <- arms(c(rho_it, lambda_it, rhoc_it, loading_it), PARfullCond, PARsupport, 1)
 model_init$GG <- rbind(cbind(diag(sampled[3], 2), diag(1, 2)),
 cbind(diag(0, 2),
 sampled[1]*matrix(c(cos(sampled[2]), sin(sampled[2]),
 -sin(sampled[2]), cos(sampled[2])), 
 nrow = 2, byrow = TRUE))
 )
 model_init$FF[1, 1] <- sampled[4]
 
 ## B. generate states - FFBS
 modFilt <- dlmFilter(y, model_init, simplify = TRUE)
 theta[] <- dlmBSample(modFilt)
 #browser()
 ## generate V
 
 y.center <- y - tcrossprod(theta[-1, , drop = FALSE], model_init$FF)
 SSy <- drop(crossprod(y.center))
 model_init$V[] <- 1/rgamma(1, shape = shape.y, rate = rate.y + 0.5 * SSy)
 
 if (!(it %% every) )
 {    #browser()
 it.save <- it.save + 1
 gibbsTheta[,,it.save] <- theta
 gibbsPAR[it.save,] <- sampled
 gibbsV[it.save] <- diag(model_init$V)
 
 }
 }
 gibbs_res <- list(dV = gibbsV, theta = gibbsTheta, par = gibbsPAR) 
 cycle <- ts(apply(t(t(gibbs_res$theta[,1,-(1:burn)])*gibbs_res$par[-(1:burn),4])
 ,1,mean)[-1],start=c(2006,4),freq=4)
 if(print_out){
 cat("\n\n")
 out1 <- gibbs_res$par[-(1:burn),]
 colnames(out1) <- c("rho", "lambda", "ar1", "alpha")
 print(format(mcmcMeans(out1), digits = 6L, scientif = FALSE), quote = FALSE)
 print(format(apply(out1, 2, quantile, probs = c(.05,.95)), digits = 6L, scientif = FALSE), quote = FALSE)
 out2 <- gibbs_res$dV[-(1:burn)]
 print(format(mean(out2), digits = 6L, scientif = FALSE), quote = FALSE)
 print(format(quantile(out2, probs = c(.05,.95)), digits = 6L, scientif = FALSE), quote = FALSE)
 out <- cbind(out1, out2)
 colnames(out) <- c(colnames(out1), "error")
 }
 }
 # Error term
 error <- y - cycle
 #browser()
 # Plot the data
 if(plot){
 plot(y,las = 1,col = "#003755", lwd = 2)
 lines(cycle, col = "#6dbace", lwd = 2)
 legend("topleft", legend = c("Observed","Cycle"), border = FALSE,
 bty = "n", col = c("#003755", "#6dbace"), lwd = 2)
 #title(main = paste(main, " - Comparison", sep = ""))
 }
 
 # Return the data
 data <- ts.union(y, cycle, error)
 options(warn = 0)
 if(method == "MLE"){
 return(invisible(list(data = data, MLE = MLE, s = Smooth_Estimates$s)))
 } else{
 if(print_out) close(pb)
 return(invisible(list(data = data, gibbs_res = gibbs_res, params = mcmcMeans(out), interval = apply(out, 2, quantile, probs = c(.05,.95)))))
 }
 }
