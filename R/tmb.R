#' @title Template Model Builder (TMB) Functions
#'
#' @description Functions for use with Template Model Builder (TMB).
#'
#'

#' @describeIn tmb Extract model parameter name(s) from TMB model.
#' @export parameters.cpp
parameters.cpp <- function(x){
   # Return model parameters declared in a TMB C++ code file.

   y <- readLines(x)                       # Read source code.
   y <- y[grep("PARAMETER", y)]            # Find lines with parameter declarations.
   y <- gsub("//.+", "", y)                # Remove comments.
   y <- gsub(";| ", "", y)                 # Remove spacing and separation characters.
   y <- gsub("PARAMETER[_A-Z]*", "", y)    # Remove parameter identifier.
   y <- gsub("[(|)]", "", y)               # Remove parentheses.
   y <- y[y != ""]

   return(y)
}

#' @describeIn tmb Extract data variable name(s) from TMB model.
#' @export data.cpp
data.cpp <- function(x){
   # Returns the data declared in a TMB C++ code file.

   y <- readLines(x)                       # Read source code.
   y <- y[grep("DATA", y)]                 # Find lines with parameter declarations.
   y <- gsub("//.+", "", y)                # Remove comments.
   y <- gsub(";| ", "", y)                 # Remove spacing and separation characters.
   y <- gsub("DATA[_A-Z]*", "", y)         # Remove parameter identifier.
   y <- gsub("[(|)]", "", y)               # Remove parentheses.
   y <- y[y != ""]

   return(y)
}

#' @describeIn tmb Update model parameters using TMB model output.
#' @export update.parameters
update.parameters <- function(x, obj, fixed, random, map){
   if (!missing(obj)){
      if (all(c("par", "fn") %in% names(obj))){
         rep <- sdreport(obj)
         fixed <- summary(rep, "fixed")[, 1]
         random <- summary(rep, "random")[, 1]
      }
   }

   # Update fixed parameters:
   vars <- unique(names(fixed))
   if (length(vars) > 0){
      for (i in 1:length(vars)){
         if (vars[i] %in% names(x)){
            v <- as.numeric(fixed[names(fixed) == vars[i]])
            if (!missing(map)) if (vars[i] %in% names(map)) v <- v[map[[vars[i]]]]
            if (length(v) == 1) v <- rep(v, length(x[[vars[i]]]))
            dim(v) <- dim(x[[vars[i]]])
            x[[vars[i]]][!is.na(v)] <- v[!is.na(v)]
         }else{
            print("'", vars[i], "' not in parameters.")
         }
      }
   }

   # Update random parameters:
   vars <- unique(names(random))
   if (length(vars) > 0){
      for (i in 1:length(vars)){
         if (vars[i] %in% names(x)){
            v <- as.numeric(random[names(random) == vars[i]])
            if (!missing(map)) if (vars[i] %in% names(map)) v <- v[map[[vars[i]]]]
            if (length(v) == 1) v <- rep(v, length(x[[vars[i]]]))
            dim(v) <- dim(x[[vars[i]]])
            x[[vars[i]]][!is.na(v)] <- v[!is.na(v)]
         }else{
            print("'", vars[i], "' not in parameters.")
         }
      }
   }

   return(x)
}

#' @describeIn tmb Update model parameter map.
#' @export update.map
update.map <- function(map, free.variables, fixed.variables){
   if (!missing(free.variables)){
      free.variables <- free.variables[which(free.variables %in% names(map))]
      if (length(free.variables) > 0){
         for (i in 1:length(free.variables)){
            map[[free.variables[i]]] <- factor(1:length(map[[free.variables[i]]]))
         }
      }
   }

   if (!missing(fixed.variables)){
      fixed.variables <- fixed.variables[which(fixed.variables %in% names(map))]
      if (length(fixed.variables) > 0){
         for (i in 1:length(fixed.variables)){
            map[[fixed.variables[i]]] <- factor(rep(NA, length(map[[fixed.variables[i]]])))
         }
      }
   }

   return(map)
}
