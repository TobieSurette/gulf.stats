library(gulf)

# Define basis functions:
istep <- function(x){
   v <- rep(0, length(x))
   v[x > 0] <- x[x > 0]
   return(v)
}


quadratic <- function(x){
   v <- 0.5 * (x * (x + 1) + 0.25)
   v[x < -0.5] <- 0
   v[x > 0.5] <- x[x > 0.5]
   return(v)
}

ilogistic <- function(x) return(0.25*log(1 + exp(4*x)))   
   

xp <- 0
w <- 1
t <- seq(-1.25, 1.25, len = 1000)

# Plot integrated bases:
clg()
windows()
lty <- c("dotted", "dashed", "solid")
cols <- c(grey(0), grey(0.3), grey(0.6))
plot(range(t), c(0, max(t)), type = "n", xlab = "x", ylab = "y", cex.lab = 1.5)
lines(t, w*istep((t-xp)/w), col = cols[1], lwd = 3, lty = lty[1])
lines(t, w*quadratic((t-xp)/w), col = cols[2], lwd = 3, lty = lty[2]) 
w <- 1*w
lines(t, w*ilogistic((t-xp)/w), col = cols[3], lwd = 3, lty = lty[3]) 
# Window width:
lines(c(xp-w/2, xp-w/2), par("usr")[3:4], lwd = 2, lty = "dashed")
lines(c(xp+w/2, xp+w/2), par("usr")[3:4], lwd = 2, lty = "dashed")

legend("topleft", 
       legend = c("i-step", "quadratic", "i-logistic"), 
       lwd = 3, lty = lty, 
       col = cols, 
       cex = 1.5,
       bg = "white",
       border = "white")
       
h <- 0.65
arrows(xp, h, xp-w/2, h, length = 0.15, lwd = 2)
arrows(xp, h, xp+w/2, h, length = 0.15, lwd = 2)
text(par("usr")[1] + diff(par("usr")[1:2]) / 2, h - 0.025 * diff(par("usr")[3:4]), "w", font = 3, cex = 1.40)


# Define basis functions:
step <- function(x){
   v <- rep(0, length(x))
   v[x > 0] <- 1
   return(v)
}


dquadratic <- function(x){
   v <- x + 0.5
   v[x < -0.5] <- 0
   v[x > 0.5] <- 1
   return(v)
}

logistic <- function(x) return(1/(1 + exp(-4*x)))  

# Plot bases:
clg()
windows()
lty <- c("dotted", "dashed", "solid")
cols <- c(grey(0), grey(0.3), grey(0.6))
plot(range(t), c(0, max(t)), type = "n", xlab = "x", ylab = "y", cex.lab = 1.5)
lines(t, w*step((t-xp)/w), col = cols[1], lwd = 3, lty = lty[1])
lines(t, w*dquadratic((t-xp)/w), col = cols[2], lwd = 3, lty = lty[2]) 
w <- 1*w
lines(t, w*logistic((t-xp)/w), col = cols[3], lwd = 3, lty = lty[3]) 
# Window width:
lines(c(xp-w/2, xp-w/2), par("usr")[3:4], lwd = 2, lty = "dashed", col = "red")
lines(c(xp+w/2, xp+w/2), par("usr")[3:4], lwd = 2, lty = "dashed", col = "red")

legend("topleft", 
       legend = c("step", "p-linear", "logistic"), 
       lwd = 3, lty = lty, 
       col = cols, 
       cex = 1.5,
       bg = "white",
       border = "white")
       
h <- 0.5
arrows(xp, h, xp-w/2, h, length = 0.15, lwd = 2)
arrows(xp, h, xp+w/2, h, length = 0.15, lwd = 2)
text(par("usr")[1] + diff(par("usr")[1:2]) / 2, h - 0.025 * diff(par("usr")[3:4]), "w", font = 3, cex = 1.40)



