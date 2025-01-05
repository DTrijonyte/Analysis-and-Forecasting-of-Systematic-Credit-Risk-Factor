
 # Description:   Create new directory, check if the one already exists

 # Arguments: 
 #   dir - path to the directory to be created
 
 function(dir) 
 {
 if(file.exists(dir)) {
 if(!file.info(dir)$isdir) {
 dir.create(dir, recursive = TRUE)
 cat("\n Created  directory ", dir, "\n")
 }
 }
 else {
 dir.create(dir, recursive = TRUE)
 cat("\n Created  directory " ,dir, "\n")
 } 
 }
