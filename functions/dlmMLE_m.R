
 # Description:   modified version of dlmMLE function to allow for derivative free optimization method

 # Arguments:
 #   y      - function to optimize
 #   parm   - optimization parameters
 #   method - maximum likelihood optimization method
 
 
 function (y, parm, build, method = "L-BFGS-B", ..., debug = FALSE) 
 { 
 logLik <- function(parm, ...) {
 mod <- build(parm, ...)
 return(dlmLL(y = y, mod = mod, debug = debug))
 }
 if(method == "NMK"){ 
 out <- nmkb(parm, logLik, control = list(restarts.max = 10, trace = FALSE), ...)
 }else{
 out <- optim(parm, logLik, method = method, ...)
 }
 
 return(out)
 }
