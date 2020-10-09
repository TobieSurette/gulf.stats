plot.variogram <- function(x, scale, cex = 1.5, cex.axis = 1.25, xaxt = "s", yaxt = "s", 
                           xlab, ylab, xlim, ylim, show.nugget = TRUE, show.sill = TRUE, show.range = TRUE, ...){
   # PLOT.VARIOGRAM - Graphically display a variogram.
                      
   # Extract model parameters:
   res <- x$empirical
   
   # Determine y-scaling:
   v <- floor(mean(log10(res$semi.variance)))
   if (missing(scale)){
      if (v > 1) scale <- 10^v else scale <- 1
   }
   
   print(scale)
   # Define axis labels:
   if (missing(xlab)) xlab <- "Distance"
   if (!is.null(x$distance.units) & (xlab != "")) xlab <- paste0(xlab, " (", x$distance.units, ")")
   if (missing(ylab)) ylab <- "Semi-variance"
   if ((scale != 1) & (ylab != "")){
      a <- floor(mean(log10(scale)))
      ylab <- parse(text = paste("Semi-variance (10^{", a, "})")) 
   }
   
   h <- seq(0, x$max.distance, len = 1000)

   if (missing(xlim)) xlim <- c(0, x$max.distance)
   if (missing(ylim)) ylim <- c(0, 1.1*max(res$semi.variance / scale))
   
   plot(res$h, res$semi.variance / scale, pch = 21, bg = "grey", 
        xlim = xlim, ylim = ylim, xlab = "", ylab = "", xaxs = "i", yaxs = "i", 
        cex = cex, cex.axis = cex.axis, type = "n", xaxt = "n", yaxt = "n", ...)
   mtext(xlab, 1, 2.5, cex = cex)
   mtext(ylab, 2, 2, cex = cex)
   grid()
   
   lines(par("usr")[1:2], c(x$var, x$var) / scale, col = "chartreuse3", lwd = 2, lty = "dashed")
   if (xaxt != "n") axis(1)
   if (yaxt != "n") axis(2)
      
   # Plot model and parameters:
   if (!is.null(x$sill)){
      nugget <- x$nugget
      sill <- x$sill
      range <- x$range
      lines(h, x$vfun(h, nugget = nugget, sill = sill, range = range) / scale, lwd = 2, col = "red")
      if (show.sill) lines(par("usr")[1:2], c(sill, sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (show.range) lines(c(range, range), c(0, sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (show.nugget) lines(par("usr")[1:2], c(nugget, nugget) / scale, lwd = 1, lty = "dashed", col = "red")    
   }
   
   # Display empirical semi-variances:
   points(res$h, res$semi.variance / scale, pch = 21, bg = "grey", cex = 1.5)
   if (!is.null(x$lag)) text(res$h, res$semi.variance / scale, res$n, pos = 3, cex = 0.75)
   
   # Display parameter texts:
   #if (!is.null(x$sill)){
      #if (show.range){
      #   if (range > mean(par("usr")[1:2])) pos = 2 else pos = 4
      #   text(range, 0.5 * sill / scale, paste0("range = ", round(range,2), " ", x$distance.units), cex = 1.25, pos = pos)                                                  
      #}
      #if (show.sill) text(0, sill / scale + 0.025*diff(par("usr")[3:4]), paste0("sill = ", round((sill) / scale,3)), cex = 1.25, pos = 4)
      #if (show.nugget & (nugget / sill) > 0.01) text(0, nugget / scale - 0.025*diff(par("usr")[3:4]), paste0("nugget = ", round(nugget / scale,3)), cex = 1.25, pos = 4)
   #}   
   #text(0, (x$var / scale) + 0.030*diff(par("usr")[3:4]), paste0("Var(sample) = ", round(x$var / scale,3)), cex = 1.25, pos = 4, col = "chartreuse3") 
   #print(str(x))
   #v <- floor(log10(par("usr")[4]))
   str <- ""
   if (is.finite(range) & !is.na(range))
   str <- c(str, paste0("range = ", round(range,1), x$distance.units)) 
   str <- c(str, paste0("sill = ", round(sill/ 10^v ,3), " x 10^", v),
                 paste0("nugget = ", round(nugget / 10^v,3), " x 10^", v),
                 paste0("Var[z] = ", round(x$var / 10^v,3), " x 10^", v))
            
   #parse(text = )         
            
   dx <- diff(par("usr")[1:2])
   dy <- diff(par("usr")[3:4])
   for (i in 1:length(str)){
      text(par("usr")[1] + 0.55 * dx, 
           par("usr")[3] + 0.30 * dy - 0.06 * (i-1) * dy,
           str[i], pos = 4)
   }
                            
   box() 
}
