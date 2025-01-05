# Source:        https://otexts.com/fpp2/backcasting.html
 # Description:   Function to reverse a forecast
 
 # Attributes:
 #  object -- a forecast object
 
 
 function(object)
 {
 h <- length(object[["mean"]])
 f <- frequency(object[["mean"]])
 object[["x"]] <- reverse_ts(object[["x"]])
 object[["mean"]] <- ts(rev(object[["mean"]]),
 end=tsp(object[["x"]])[1L]-1/f, frequency=f)
 object[["lower"]] <- object[["lower"]][h:1L,]
 object[["upper"]] <- object[["upper"]][h:1L,]
 return(object)
 }
