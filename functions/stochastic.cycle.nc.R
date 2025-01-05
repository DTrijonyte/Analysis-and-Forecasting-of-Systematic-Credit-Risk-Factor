
 # Description:   Stochastic cycle iterations
 
 # Arguments:
 #   cycle  - cycle object (any)
 #   tol    - tolerance level
 #   max.it - max number of iterations
 #   m0     - set the prior value for the cycle's starting point
 #   cyc    - the initial/final value of cycle's length
 #   rho    - the initial value of rho
 #   C0     - prior for variability of the cycle's innovations
 #   nc     - order of a stochastic cycle
 
 function(cycle, tol = 0.001, max.it = 100, m0 = cycle[1], 
 cyc = 32, rho = 0.8, C0 = 1e-4, 
 silent = FALSE,
 nc = 1){
 PI <- base::pi
 # State space description of the cycle object
 model <- function(X){
 #CYCLE 1st order, equivalent to ARMA(2, 1) 
 T_matrix <- X[1]*matrix(c(cos(2*PI/X[2]), sin(2*PI/X[2]),
 -sin(2*PI/X[2]), cos(2*PI/X[2])), 
 nrow = 2, byrow = TRUE)
 lm <- cbind(rbind(matrix(0, 2, 2*(nc - 1)), diag(2*nc - 2)), matrix(0, 2*nc, 2)) # from lmat(2*nc, 2)
 cycle <- dlm(FF = matrix(c(rep(0, (nc - 1)*2), 1, 0), nrow = 1), 
 V = 0,
 GG = kronecker(diag(nc), T_matrix) + lm,
 W = diag(c(exp(X[3]), exp(X[3]), rep(0, (nc - 1)*2))), 
 # second innovation in auxiliary cycle is restricted to 0 as in Mohr
 m0 = rep(m0, 2*nc),
 C0 = C0*diag(2*nc)
 ) 
 return(cycle)
 }
 
 #browser()
 #for(nc in  1:10){
 MLE <- try(dlmMLE_m(cycle, parm = c(rho, cyc, -1), model,
 method = "NMK", 
 lower = c(0.001, 6, -50),
 upper = c(0.999, PU, 0)))
 if(class(MLE) == "try-error") {
 cat("Error occured running nc =", nc, "\n")
 return(list(rho = NA, cyc = NA, nc = nc))
 }
 #if(MLE$message == "Successful convergence") break
 #}
 if(!silent){
 print(MLE$message)
 print(MLE$par)
 }
 #if(MLE$message == "Successful convergence"){cat("Converged at level of cycle =",nc,"\n")}else{
 #  cat("Convergence was not achieved, try changing initial values for rho and/or cyc!")
 #}
 coefs <- MLE$par[1:2]
 Smooth_Estimates <- dlmSmooth(cycle, model(MLE$par))
 cycle_fitted <- dropFirst(Smooth_Estimates$s[, 2*(nc - 1) + 1])
 #browser()
 res <- list(rho = coefs[1], cyc = coefs[2], nc = nc, data = ts.union(cycle, cycle_fitted))
 
 
 return(res)
 }
 
 
