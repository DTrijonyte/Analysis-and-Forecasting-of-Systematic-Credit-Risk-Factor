# Description:   Calculates the number of trailing NA values in a vector or time series.

# Attributes:
#  x -- a vector or time series that may contain NA values.

right_na <- function(x) {
  N <- length(x)
  right_na <- N - length(na.trim(x, "right"))
  return(right_na)
}
