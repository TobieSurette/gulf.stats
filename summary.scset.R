summary.scset <- function(x, year, tow.id, category, ...){
   # SUMMARY.SCSET - Generate a summary of snow crab survey tow data.

   # Load tow data:
   if (missing(x) & !missing(year)){
      x <- read.scset(year = year)
      if (missing(tow.id)){
         tow.id <- unique(x$tow.id)
         flag <- TRUE
      }
   }
      
   if (!missing(category)){
      if (!is.character(category)) stop("'category' must be a vector of character strings.")
      b <- read.scbio(year = unique(x$year))
      b <- b[b$tow.id != "" , ]  
      b <- b[(x$month >= 7), ]
      b <- b[substr(x$tow.id,2,2) != "C", ]
      res <- summary(b, category = category, by = c("year", "tow.id"), ...)
      
      # Remove biological data not in the set data:
      index <- gulf::match(res[c("year", "tow.id")], x[c("year", "tow.id")]) 
      res <- res[!is.na(index), ]
      
      x[category] <- 0
     
      #if (any(is.na(index))) stop("Some biological data was not matched to the tow data.")
      index <- gulf::match(res[c("year", "tow.id")], x[c("year", "tow.id")]) 
      
      x[index, category] <- res[category]
      
      return(x)
   }
   
   if (missing(x) & missing(year)) stop("'x' or 'year' must be specified.")

   # Define flag whether to load all probe data:
   flag <- FALSE
   

      
   # Define 'year':
   if (missing(year)) year <- sort(unique(x$year))
   
   # Subset of tows:
   if (!missing(tow.id) & !missing(x)) x <- x[x$tow.id %in% toupper(tow.id), ]
   if (missing(tow.id)){
      tow.id <- x$tow.id
      flag <- FALSE
   }
   
   if (length(year) != 1) stop("'year' must be a scalar.")
   
   # Read probe files:
   if (flag){
      esonar   <- read.esonar(year = year)
      minilog  <- read.minilog(year = year)
      headline <- read.star.oddi(year = year, type = "depth")
      tilt     <- read.star.oddi(year = year, type = "tilt")
      tow.id   <- unique(c(esonar$tow.id, minilog$tow.id, headline$tow.id, tilt$tow.id))
      tow.id   <- tow.id[!is.na(tow.id)]
   }else{
      esonar   <- read.esonar(year = year, tow.id = tow.id)
      minilog  <- read.minilog(year = year, tow.id = tow.id)
      headline <- read.star.oddi(year = year, type = "depth", tow.id = tow.id)
      tilt     <- read.star.oddi(year = year, type = "tilt", tow.id = tow.id)   
   }
   
   # Tow table for each probe type:
   vars <- c("year", "month", "day", "tow.id")
   res <- unique(rbind(esonar[vars], minilog[vars], headline[vars], tilt[vars]))
   res <- res[!is.na(res$tow.id), ]
   res <- sort(res, by = vars)
   rownames(res) <- NULL
   res[c("esonar", "minilog", "headline", "tilt")] <- FALSE

   #if (length(esonar) > 0)   res$esonar[!is.na(match.data.frame(res[vars], unique(esonar[vars])))] <- TRUE
   #if (length(minilog) > 0)  res$minilog[!is.na(match.data.frame(res[vars], unique(minilog[vars])))] <- TRUE
   #if (length(headline) > 0) res$headline[!is.na(match.data.frame(res[vars], unique(headline[vars])))] <- TRUE
   #if (length(tilt) > 0)     res$tilt[!is.na(match.data.frame(res[vars], unique(tilt[vars])))] <- TRUE

   # Number of observations:
   for (i in 1:nrow(res)){
      res$esonar[i]   <- sum(esonar$tow.id == res$tow.id[i], na.rm = TRUE)
      res$minilog[i]  <- sum(minilog$tow.id == res$tow.id[i], na.rm = TRUE)
      res$headline[i] <- sum(headline$tow.id == res$tow.id[i], na.rm = TRUE)
      res$tilt[i]     <- sum(tilt$tow.id == res$tow.id[i], na.rm = TRUE)
   }

   # Check that time time intervals overlap:
   res$overlap <- NA
   for (i in 1:nrow(res)){
      if ((i %% 10) == 0) cat(paste0("Tows ", i-9, " to ", i, ".\n"))
      int <- list()
      index <- which(esonar$tow.id == res$tow.id[i])
      if (length(index) > 0) int <- c(int, list(range(time(esonar[index, ]))))
      index <- which(minilog$tow.id == res$tow.id[i])
      if (length(index) > 0) int <- c(int, list(range(time(minilog[index, ]))))
      index <- which(headline$tow.id == res$tow.id[i])
      if (length(index) > 0) int <- c(int, list(range(time(headline[index, ]))))
      index <- which(tilt$tow.id == res$tow.id[i])
      if (length(index) > 0) int <- c(int, list(range(time(tilt[index, ]))))   

      d <- min(unlist(lapply(lapply(int, function(x) x[2]), as.numeric))) - 
           max(unlist(lapply(lapply(int, function(x) x[1]), as.numeric)))
      res$overlap[i] <- d
   }

   # Concatenate 'x' and 'res':
   key <- c("year", "tow.id")
   
   index <- match.data.frame(res[key], x[key])
   res$tow.number <- x$tow.number[index]
   res$valid <- x$valid[index]
   index <- match(x[key], res[key])
   if (sum(is.na(index)) > 0){
      tmp <- x[is.na(index), ]
      tmp[setdiff(names(res), names(x))] <- 0
      tmp <- tmp[names(res)]
      res <- rbind(res, tmp)
   }
   
   # Calculate distance between E-Sonar coordinates and logbook coordinates:
   if (length(esonar) > 0){
      tows <- unique(esonar$tow.id)
      res$esonar.distance <- NA
      for (i in 1:length(tows)){
         index <- esonar$tow.id == tows[i] 
         ii <- which(x$tow.id == tows[i])
         if (length(ii) > 0){
            d <- min(distance(longitude(x[ii,]), latitude(x[ii,]), esonar$longitude[index], esonar$latitude[index]))
            res$esonar.distance[match(tows[i], res$tow.id)] <- round(d, 4)
         }
      }
   }
   
   # Re-order columns:
   vars <- c("year", "month", "day", "tow.id", "tow.number", "valid")
   res <- res[c(vars, setdiff(names(res), vars))]
   
   return(res)
}
