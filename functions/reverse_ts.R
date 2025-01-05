# Source:        https://otexts.com/fpp2/backcasting.html
# Description:   Function to reverse time

# Attributes:
#   y -- time series object

function(y) {
  ts(rev(y), start=tsp(y)[1L], frequency=frequency(y))
}
