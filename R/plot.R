#' Plot Statistical Models
#'
#' @description Functions to graphically display statistical models.
#'
#' @param x \code{variogram} object.
#' @param scale Scalar by which to scale the variogram y axis.
#' @param cex Text size.
#' @param cex.axis Text size for plot axes.
#' @param xaxt,yaxt \code{\link{par}} graphics parameters.
#' @param xlab,ylab Character strings specifying the axes labels.
#' @param xlim,ylim Two-element numerical vectors specifying plot axes limits.
#' @param nugget,sill,range Logical value specifying whether to display reference lines for variogram parameters.
#' @param discrete Logical value specifying whether morphometric data are rounded.

#' @describeIn plot Variogram plot method.
#' @export
plot.variogram <- function(x, scale, cex = 1.5, cex.axis = 1.25, xaxt = "s", yaxt = "s",
                           xlab, ylab, xlim, ylim, nugget = TRUE, sill = TRUE, range = TRUE, ...){

   # Extract model parameters:
   res <- x$empirical

   # Determine y-scaling:
   v <- floor(mean(log10(res$semi.variance)))
   if (missing(scale)) if (v > 1) scale <- (10^v) else scale <- 1

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
      lines(h, x$vfun(h, nugget = x$nugget, sill = x$sill, range = x$range) / scale, lwd = 2, col = "red")
      if (sill) lines(par("usr")[1:2], c(x$sill, x$sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (range) lines(c(x$range, x$range), c(0, x$sill) / scale, lwd = 1, lty = "dashed", col = "red")
      if (nugget) lines(par("usr")[1:2], c(x$nugget, x$nugget) / scale, lwd = 1, lty = "dashed", col = "red")
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

   print(x$range)
   str <- ""
   if (is.finite(x$range) & !is.na(x$range)) str <- c(str, paste0("range = ", round(x$range,1), x$distance.units))
   str <- c(str, paste0("sill = ", round(x$sill / 10^v ,3), " x 10^", v),
                 paste0("nugget = ", round(x$nugget / 10^v,3), " x 10^", v),
                 paste0("Var[z] = ", round(x$var / 10^v,3), " x 10^", v))

   dx <- diff(par("usr")[1:2])
   dy <- diff(par("usr")[3:4])
   for (i in 1:length(str)){
      text(par("usr")[1] + 0.55 * dx,
           par("usr")[3] + 0.30 * dy - 0.06 * (i-1) * dy,
           str[i], pos = 4)
   }

   box()
}

# @describeIn plot Generic morphometric plot method.
#' @export plot.morphometry
plot.morphometry <- function(x, ...) UseMethod("plot.morphometry")

#' @describeIn plot Morphometric plot for snow crab biological data.
#' @rawNamespace S3method(plot.morphometry,scsbio)
plot.morphometry.scsbio <- function(x, y, theta, xlim = c(10, 140), log = TRUE, title, discrete = FALSE, ...){
   v <- morphometry.scsbio(x, y, theta = theta)
   x0 <- seq(xlim[1], xlim[2], len = 1000)
   v0 <- morphometry.scsbio(x0, theta = theta)

   # Prepare plotting area:
   layout(rbind(0, 0, cbind(0, kronecker(c(1,1,2:4), matrix(1, nrow = 5, ncol = 5)), 0), 0, 0))
   par(mar = c(0, 0, 0, 0))

   # Plot data and regression curves:
   plot(x, y, type = "n", xlab = "", ylab = "", xlim = xlim, ylim = c(0, 40), cex.lab = 1.25, xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   if (!discrete){
      points(x[v$p_mature_posterior >= 0.5], y[v$p_mature_posterior >= 0.5], pch = 19, col = "deepskyblue1", cex = 0.1)
      points(x[v$p_mature_posterior < 0.5], y[v$p_mature_posterior < 0.5], pch = 19, col = "chartreuse2", cex = 0.1)
   }else{
      points(jitter(x[v$p_mature_posterior >= 0.5], amount = 0.5),
             jitter(y[v$p_mature_posterior >= 0.5], amount = 0.5), pch = 19, col = "deepskyblue1", cex = 0.1)
      points(jitter(x[v$p_mature_posterior < 0.5], amount = 0.5),
             jitter(y[v$p_mature_posterior < 0.5], amount = 0.5), pch = 19, col = "chartreuse2", cex = 0.1)
   }
   lines(x0, exp(v0$mu_immature), lwd = 2, col = "chartreuse3")
   lines(x0, exp(v0$mu_mature), lwd = 2, col = "deepskyblue3")
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Chela height (mm)", 2, 2.5, cex = 1.0)
   if (!missing(title)) mtext(title, 3, 1.0, cex = 1.25)
   legend("topleft",
          legend = c("Mature", "Immature"),
          col = c("deepskyblue3", "chartreuse3"),
          lwd = 2)
   box()

   # Plot maturity proportions:
   cw <- round(x)
   res <- aggregate(list(k = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), sum)
   res$n <- aggregate(list(n = (v$p_mature_posterior  >= 0.5)), by = list(cw = cw), length)$n
   res$p_mature <-  res$k / res$n
   plot(xlim, c(0, 1.1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")
   grid()
   gbarplot(res$p_mature, res$cw, col = "grey90", width = 1, add = TRUE)
   lines(x0, v0$p_mature, col = "red", lwd = 2)
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   legend("topleft",
          legend = c("model", "classification"),
          col = c("red", "grey30"),
          pch = c(NA, 22),
          pt.bg = c(NA, "grey60"),
          pt.cex = c(1,2.5),
          lwd = c(2, 0))
   box()
   mtext("Mature Proportion", 2, 2.5, cex = 1.0)

   # Residual plots
   if (!discrete){
      r_imm <- (log(y[v$p_mature_posterior < 0.5]) - v$mu_immature[v$p_mature_posterior < 0.5]) / exp(theta[["log_sigma"]])
      r_mat <- (log(y[v$p_mature_posterior >= 0.5]) - v$mu_mature[v$p_mature_posterior >= 0.5]) / exp(theta[["log_sigma"]])
   }else{
      r_imm <- (log(jitter(y[v$p_mature_posterior < 0.5], amount = 0.5)) - v$mu_immature[v$p_mature_posterior < 0.5]) / exp(theta[["log_sigma"]])
      r_mat <- (log(jitter(y[v$p_mature_posterior >= 0.5], amount = 0.5)) - v$mu_mature[v$p_mature_posterior >= 0.5]) / exp(theta[["log_sigma"]])
   }
   res$p_model <- morphometry.scsbio(res$cw, theta = theta)$p_mature
   res$logit_p_model <- log(res$p_model / (1-res$p_model))
   res$logit_p <- log(res$p_mature / (1-res$p_mature))
   res$residuals_logit_p <- res$logit_p - res$logit_p_model

   # Mature residual plot:
   if (!discrete){
      plot(x[v$p_mature_posterior >= 0.5], r_mat, pch = 19, xlab = "", ylab = "", xlim = xlim, col = "deepskyblue1", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }else{
      plot(jitter(x[v$p_mature_posterior >= 0.5], amount = 0.5), r_mat, pch = 19, xlab = "", ylab = "", xlim = xlim, col = "deepskyblue1", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }
   grid()
   abline(0, 0, col = "deepskyblue1", lwd = 2)
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Mature residuals", 2, 2.5, cex = 1.0)
   box()

   # Immature residual plot:
   if (!discrete){
      plot(x[v$p_mature_posterior < 0.5], r_imm,  xlab = "",  ylab = "", xlim = xlim, pch = 19, col = "chartreuse2", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }else{
      plot(jitter(x[v$p_mature_posterior < 0.5], amount = 0.5), r_imm,  xlab = "",  ylab = "", xlim = xlim, pch = 19, col = "chartreuse2", bg = "grey", cex = 0.1, xaxt = "n", xaxs = "i", ylim = c(-3,3))
   }
   grid()
   abline(0, 0, col = "chartreuse3", lwd = 2)
   gulf.graphics::vline(95, lty = "dashed", col = "red")
   mtext("Immature residuals", 2, 2.5, cex = 1.0)
   box()

   axis(1)
   mtext("Carapace width (mm)", 1, 2.5, cex = 1.15)
}


