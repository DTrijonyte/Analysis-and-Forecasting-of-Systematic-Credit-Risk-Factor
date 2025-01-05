 # Description:   Modified accuracy function from forecast

 # Arguments:
 #     f - predicted values
 #     x - actual values
 #     y - insample actual values for scaling errors
 #     D - lag for scaling difference (insample rwf or naive srwf)
 #   rwf - random walk forecasts of the lengths f and x
 # for simplicity function do not have any internal length checks!
 
 function(f, x, y = NULL, D = 1, rwf = NULL)
 {
 res <- x - f
 pe  <- res/x * 100
 spe <- res/(x + f) * 200
 ## Scale dependent measures
 rmse <- sqrt(mean(res^2, na.rm = TRUE))
 mae  <- mean(abs(res), na.rm = TRUE)
 mdae <- median(abs(res), na.rm = TRUE)
 ## Scale independent measures
 mape  <- mean(abs(pe), na.rm = TRUE)
 mdape <- median(abs(pe), na.rm = TRUE)
 rmspe <- sqrt(mean(pe^2, na.rm = TRUE))
 rmdspe <- sqrt(median(pe^2, na.rm = TRUE))
 ## Symmetric measures
 smape  <- mean(abs(spe), na.rm = TRUE)
 smdape <- median(abs(spe), na.rm = TRUE)
 out <- c(rmse,mae,mdae, mape, mdape, rmspe, rmdspe, smape, smdape)
 names(out) <- c("RMSE","MAE","MdAE","MAPE","MdAPE","RMSPE","RMdSPE","sMAPE","sMdAPE")
 if(!is.null(rwf)){
 ## RelMAE
 relmae <- mae/mean(abs(x - rwf), na.rm = TRUE)
 out  <- c(out, mase)
 names(out)[length(out)] <- "RelMAE"
 }
 if(!is.null(y)){
 ## MASE
 q <- res/mean(abs(diff(y, lag = D)), na.rm = TRUE) # scaled errors
 mase <- mean(abs(q), na.rm = TRUE)
 out  <- c(out, mase)
 names(out)[length(out)] <- "MASE"
 }
 return(out)
 }
