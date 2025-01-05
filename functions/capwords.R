
 # Description:   Capitalizing first letters in descriptions and titles (from ?toupper)
 #                Use toupper or tolower to make all letters upper or lower case
 
 # Arguments:
 #  s      - string
 #  strict - only first letter will be capital 
 
 function(s, strict = FALSE) 
 {
 cap    <- function(s) paste(toupper(substring(s, 1, 1)),
 {s <- substring(s, 2); if(strict) tolower(s) else s},
 sep = "", collapse = " ")
 sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
 }
