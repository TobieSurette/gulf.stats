read.sclen <- function(x, ...){
   # READ.SCLEN - Reads an ASCII gulf length card.

   # Get file list:
   files <- sclen.file.str(...)
   
   # Read data:
   x <- read.csv(file = files, header = TRUE, stringsAsFactors = FALSE)
   
   return(x)
}
