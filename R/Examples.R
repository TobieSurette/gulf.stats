path <- "U:/Statistical Models/R Object Classes/splm class/"
#U:\Statistical Models\R Object Classes\splm class
source(paste0(path, "source splm.R"))

# Generate sample data set:
x <- 10 * sort(runif(100))
y <- x * NA
y[x < 5] <- x[x<5] + 0.2 * rnorm(sum(x<5))
y[x > 5] <- 0.2 * rnorm(sum(x>=5)) + 5

# Things that should work:

# 'splm' examples:
splm()
splm(order = 3)
splm(slope = 1:5)
splm(intercept = 1, slope = c(-1, 0, 1, 0, 1))
splm(transition = 1:5)
splm(transition = 1:5, transition.basis = "step")
splm(intercept = 0, slope = c(0, 1), transition = c(0), window = 0.1, transition.basis = "logistic")
splm(intercept = 0, slope = c(0, -1, 1), transition = c(0, 1), window = 0.1, transition.basis = "linear")
splm(intercept = 0, slope = c(0, -1, 1), transition = c(0, 1), window = c(0.1, 0.2), transition.basis = "linear")
splm(intercept = 0, slope = c(0, -1, 1), transition = c(0, 1), window = c(0.1, 0.2), sigma = 1, transition.basis = "linear")
splm(intercept = 0, slope = c(0,1), transition = 0, window = 0.1, transition.basis = "logistic", error.model="expplus", sigma=c(0,0,0))
splm(intercept = 0, slope = c(0,1), transition = 0, window = 0.1, transition.basis = "logistic", error.model = "splm", sigma = c(0,0))

# 'theta' extraction examples:
x <- splm(intercept = 0, slope = c(0, -1, 1), transition = c(0, 1), window = c(0.1, 0.2), sigma = 1, transition.basis = "linear")
fix(x) <- "transition"  # Fix 'transition' parameters.
theta(x) # Free parameters.
theta(x, full = TRUE) # Full parameter vector.
theta(x, transform = TRUE) # Back-transformed transition parameters.
theta(x, "transition")     # Block of transition parameters.
theta(x, "transition", transform = TRUE)   # Block of back-transformed transition parameters.
theta(x, "window2", transform = TRUE) # Particular parameter.
theta(x, "intercept")
theta(x, 2:5)  # Second to fifth parameters from full parameter.

# 'theta' assignment examples:
theta(x, "window") <- 0.5  # Assign 0.5 to all 'window' parameters.
theta(x, "window2") <- 0.1 # Change second parameter.
theta(x, "intercept") <- 3 
theta(x) <- theta(x)  # Full vector assignment.
p <- theta(x)
p[1] = 5
theta(x) <- p  # Assign parameter vector.
theta(x) <- NULL  # Erase free parameter values.
theta(x) <- p  # Reassign parameter vector.
theta(x) <- NA  # Another way of erasing parameters.
theta(x) <- p  # Restore.
theta(x, full = TRUE) <- NA # Erase all parameters.

# Hockey stick plot:
windows(height = 7, width = 9)
t <- seq(0, 2, len = 1000)
v <- splm(intercept = 0, slope = c(1, 0), transition = c(1), window = 0.00001, transition.basis = "logistic")
plot(t, predict(v, t), type = "l", xlab = "x", ylab = "y", 
     cex.lab = 1.5, lwd = 2, ylim = c(0, 1.2), xaxs = "i", yaxs = "i", cex.axis = 1.3)
xp <- theta(v, "transition")
lines(c(xp, xp), c(0, predict(v, xp)), lwd = 2, lty = "dashed", col = "red") 
text(0.5, predict(v, 0.5), paste0("slope = ", theta(v, "slope0")), srt = 47, pos = 3, cex = 1.5)
text(1.5, predict(v, 1.5), paste0("slope = ", theta(v, "slope1")), srt = 0, pos = 3, cex = 1.5)
arrows(xp+0.1, 0.5, xp+0.01, 0.5, len = 0.15, lwd = 2)
text(xp+0.1, 0.5, expression(x[p]), srt = 0, pos = 4, cex = 1.5)
arrows(0.3, 0.1, 0.04, 0.0133, len = 0.15, lwd = 2)
text(0.3, 0.1, paste0("intercept = ", theta(v, "intercept")), srt = 0, pos = 4, cex = 1.5)


# Example with a single change-point:
t <- seq(-1.75, 1.75, len = 1000)
v <- splm(intercept = 0, slope = c(0, 1), transition = c(0), window = 0.1, transition.basis = "logistic")

windows(width = 10, height = 7)
plot(t, predict(v, t), type = "n", lwd = 2, xlab = "x", ylab = "y", cex.lab = 1.5)  
grid()
w <- seq(1.5, 0, by = -0.3)
cols <- rainbow(length(w))
cols <- rev(grey(seq(.9, 0, len = length(w))))
cols <- colorRampPalette(c("black", "grey95", "black"))(length(w))
lty = c("dotted", "solid", "dashed")
lty <- lty[1:length(w) %% length(lty) + 1]
cols[length(cols)] <- "red"
lty[length(lty)] <- "solid"
for (i in  1:length(w)){
   if (w[i] == 0){
      lines(t, predict(v, t, smooth = FALSE), lwd = 3, col = cols[i], lty = lty[i])
   }else{
      theta(v, "window") <- log(w[i])
      lines(t, predict(v, t), lwd = 3, col = cols[i], lty = lty[i])
   }
}
theta(v, "window") <- log(w[1])
lines(t, predict(v, t), lwd = 3, col = cols[1], lty = lty[1])
legend("topleft", legend = rev(paste("w = ", w)), lwd = 3, col = rev(cols), lty = rev(lty), bg = "white", cex = 1.3)

 
# Example with two change-points:
t <- seq(-1.5, 2.75, len = 1000)
v <- splm(intercept = 0, slope = c(0, -1, 1), transition = c(0, 1), window = 0.1, transition.basis = "linear")

windows(width = 10, height = 7)
plot(t, predict(v, t), type = "n", lwd = 2, xlab = "x", ylab = "y", cex.lab = 1.5)  
grid()
w <- seq(1.5, 0, by = -0.3)
cols <- rainbow(length(w))
cols <- rev(grey(seq(.9, 0, len = length(w))))
cols <- colorRampPalette(c("black", "grey95", "black"))(length(w))
lty = c("dotted", "solid", "dashed")
lty <- lty[1:length(w) %% length(lty) + 1]
cols[length(cols)] <- "red"
lty[length(lty)] <- "solid"
for (i in  1:length(w)){
   if (w[i] == 0){
      lines(t, predict(v, t, smooth = FALSE), lwd = 3, col = cols[i], lty = lty[i])
   }else{
      theta(v, "window") <- log(w[i])
      lines(t, predict(v, t), lwd = 3, col = cols[i], lty = lty[i])
   }
}
theta(v, "window") <- log(w[1])
lines(t, predict(v, t), lwd = 3, col = cols[1], lty = lty[1])
legend("topleft", legend = rev(paste("w = ", w)), lwd = 3, col = rev(cols), lty = rev(lty), bg = "white", cex = 1.3)

# Example with a multiple change-points:
t <- seq(-2, 3, len = 1000)
v <- splm(intercept = 0, slope = c(-1, 1, -1, 1, -1), transition = c(-1, 0, 1, 2), window = 0.1, transition.basis = "logistic")

windows(width = 10, height = 7)
plot(t, predict(v, t), type = "l", lwd = 2, xlab = "x", ylab = "y", cex.lab = 1.5)  
grid()
w <- seq(1.5, 0, by = -0.3)
cols <- rainbow(length(w))
cols <- rev(grey(seq(.9, 0, len = length(w))))
cols <- colorRampPalette(c("black", "grey95", "black"))(length(w))
lty = c("dotted", "solid", "dashed")
lty <- lty[1:length(w) %% length(lty) + 1]
cols[length(cols)] <- "red"
lty[length(lty)] <- "solid"
for (i in  1:length(w)){
   if (w[i] == 0){
      lines(t, predict(v, t, smooth = FALSE), lwd = 3, col = cols[i], lty = lty[i])
   }else{
      theta(v, "window") <- log(w[i])
      lines(t, predict(v, t), lwd = 3, col = cols[i], lty = lty[i])
   }
}
theta(v, "window") <- log(w[1])
lines(t, predict(v, t, smooth = FALSE), lwd = 3, col = cols[1], lty = lty[1])
legend("topleft", legend = rev(paste("w = ", w)), lwd = 3, col = rev(cols), lty = rev(lty), bg = "white", cex = 1.3)


t <- seq(-2, 3, len = 1000)
v <- splm(intercept = 0, slope = c(-1, 2, -1, 2, -1), transition = c(-1, 0, 1, 2), window = 0.1, transition.basis = "logistic")

plot(v)
theta(v, "window") <- log(0.4)
plot(v, col = "red", add = TRUE, lwd = 2)
theta(v, "window") <- log(0.8)
plot(v, col = "blue", add = TRUE, lwd = 2)

theta(v, "window") <- log(0.000000001)
lines(t, predict(v, t), col = "green", lwd = 2)  

theta(v, "window") <- log(1)

predict(v, t, smooth = FALSE)


