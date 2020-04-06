cluster <- function(M){
   # Cluster using 'M' adjacency matrix:

   # Check that adjacency matrix:
   if (!is.matrix(M)) stop("Adjacency matrix must be a matrix.")
   if (dim(M)[1] != dim(M)[2]) stop("Adjacency matrix must be square.")
   if (!is.logical(M)){
      if (!all(unique(as.numeric(M)) %in% 0:1))
         stop("Adjacency matrix must be logical or contain only zeroes or ones.")
      M <- M == 1
   }

   # Convert adjacency matrix to a grouped list:
   grouping <- apply(M, 1, function(x) which(x))

   changed <- TRUE
   while (changed){
      temp <- list()
      temp[[1]] <- grouping[[1]]
      k <- 2
      for (i in 2:length(grouping)){
         flag <- TRUE
         j <- 1
         while (flag & (j <= length(temp))){
            if (length(intersect(grouping[[i]], temp[[j]])) > 0){
               temp[[j]] <- union(temp[[j]], grouping[[i]])
               flag <- FALSE
            }else{
               j <- j + 1
            }
         }
         if (flag){
            temp[[k]] <- grouping[[i]]
            k <- k + 1
         }
      }
      # Check if 'grouping' and 'temp' are identical:
      if (length(grouping) == length(temp)){
         index <- rep(FALSE, length(grouping))
         for (i in 1:length(grouping)){
            if (!(all(grouping[[i]] %in% temp[[i]]) & all(temp[[i]] %in% grouping[[i]]))) index[i] <- TRUE
         }
         if (all(!index)) changed <- FALSE
      }
      if (changed) grouping <- temp
   }

   # Convert 'grouping' list to vector form:
   group <- rep(NA, dim(M)[1])
   for (i in 1:length(grouping)){
      group[grouping[[i]]] <- i
   }

   return(group)
}