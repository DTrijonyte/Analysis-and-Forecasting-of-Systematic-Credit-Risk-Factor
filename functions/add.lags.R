 # Description:   Adds lags and leads to the data.frame

 # Arguments:
 #   data     - data frame or time series matrix
 #   max.lag  - maximum number of lags to be added
 #   max.lead - maximum number of leads to be added
 #   do.melt  - produce stacked long data frames of the variables
 #   by       - each which lag or lead should be added
 
 function (data, max.lag, max.lead = 0, do.melt = FALSE, by = 1) 
 {
 if(do.melt)
 l <- if("Date" %in% names(data)) as.matrix(data[-1]) else as.matrix(data)
 names <- colnames(l)
 if(max.lead > 0){
 res <- ll <- l <- rbind(matrix(NA, nrow = by*max.lead, ncol = ncol(l)), l,
 matrix(NA, nrow = by*max.lag,  ncol = ncol(l)))
 for (i in seq(by*1, by*max.lead, by = by)) {
 ll <- rbind(as.matrix(as.data.frame(ll[-1,])), matrix(NA,nrow = 1, ncol = ncol(ll)))
 colnames(ll) <- paste(names, ".lead", i, sep = "")
 res <- cbind(res, ll)
 }    
 }else  res <- l <- rbind(l, matrix(NA, nrow = by*max.lag, ncol = ncol(l)))
 
 for (i in seq(by*1, by*max.lag, by = by)) {
 l <- rbind(rep(NA,ncol(l)), as.matrix(as.data.frame(l[-nrow(l),])))
 colnames(l) <- paste(names, ".lag", i, sep = "")
 res <- cbind(res, l)
 }
 if(do.melt){
 res <- as.data.frame(melt(res))
 names(res) <- c("Date","List","TD")
 }
 rownames(res) <- NULL
 as.data.frame(res)
 }
