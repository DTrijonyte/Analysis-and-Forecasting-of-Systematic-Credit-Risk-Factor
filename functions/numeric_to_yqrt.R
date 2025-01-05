
 # Description:   Takes a numeric input representing a year with a fractional part 
 #                indicating the quarter, and returns a vector with the integer 
 #                part as the year and the fractional part converted to a quarter
 # Modified by:   -
 # Change log:    -
 
 # Attributes:
 #   number -- a numeric value where the integer part represents the year and the fractional
 #             part indicates the position within the year as quarters
 
 function(number){
 year <- floor(number)
 quarter <- (number - floor(number) + 0.25)*4
 return(c(year, quarter))
 }
