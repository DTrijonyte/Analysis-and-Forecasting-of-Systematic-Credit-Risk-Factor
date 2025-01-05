
 # Description:   The code identifies which lag has strongest correlation with Z for each variable

 # Arguments:
 # data   -- data frame with all data in one place
 
 max_column <- function(data) {
 col_names <- colnames(data)
 cor_columns <- grep("^Cor", col_names, value = TRUE)
 if (length(cor_columns) == 0) {
 stop("No columns matching the pattern '^Cor'.")
 }
 max_columns <- vector("character", nrow(data))
 for (row in 1:nrow(data)) {
 row_values <- data[row, cor_columns, drop = FALSE]
 max_index <- which.max(abs(as.numeric(row_values)))
 max_col_name <- cor_columns[max_index]
 lag_part <- sub("Cor ", "", max_col_name)
 max_columns[row] <- paste0(rownames(data)[row], ".", lag_part)
 }
 return(max_columns)
 }
 
