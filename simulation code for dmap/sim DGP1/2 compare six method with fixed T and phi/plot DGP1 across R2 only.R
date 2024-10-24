
R2 = seq(0.1, 0.9, by=0.1)

MSPE_R2a = read.table("MSPE across R2(N=2^13, J=2^6, p=15, T=3)(Case 1).txt", header=TRUE)
MSPE_R2b = read.table("MSPE across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt", header=TRUE)
RISK_R2a = read.table("RISK across R2(N=2^13, J=2^6, p=15, T=3)(Case 1).txt", header=TRUE)
RISK_R2b = read.table("RISK across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt", header=TRUE)

name_method = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "gMAP", "DP")


windows(height=6, width=6, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(2, 2), mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   

## Case 1, MSPE_R2 ##
plot(R2, log(MSPE_R2a[,1], 10), ylim=c(range(log(MSPE_R2a, 10))[1], 0.4), type="n", #
      xlab=expression(R^2), ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2a[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2a[, i], 10), col=i+1, pch=i+1)
}
# legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1), J=2^6, p=15", side=3, adj=0, font=3)
rect(0.48, -0.01, 0.52, 0.08, lty=1, border ="gray")  # same as xlim and ylim


## Case 2, MSPE_R2 ##
plot(R2, log(MSPE_R2b[,1], 10), ylim=c(range(log(MSPE_R2b, 10))[1], 1), type="n",
      xlab=expression(R^2), ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2b[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2b[, i], 10), col=i+1, pch=i+1)
}
mtext("(Case 2), J=2^6, p=15", side=3, adj=0, font=3)



## Case 1, Risk_R2 ##
plot(R2, RISK_R2a[,1], ylim=c(0, 0.6), type="n", #range(RISK_R2a)[1]
      xlab=expression(R^2), ylab="RISK")
for(i in 1:length(name_method)){
 lines(R2, RISK_R2a[, i], col=i+1, lty=i+1, lwd=1)
 points(R2, RISK_R2a[, i], col=i+1, pch=i+1)
}
legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1), J=2^6, p=15", side=3, adj=0, font=3)

## Case 2, Risk_R2 ##
plot(R2, RISK_R2b[,1], ylim=c(0, 1), type="n",  #range(RISK_R2b)[1]
      xlab=expression(R^2), ylab="RISK")
for(i in 1:length(name_method)){
 lines(R2, RISK_R2b[, i], col=i+1, lty=i+1, lwd=1)
 points(R2, RISK_R2b[, i], col=i+1, pch=i+1)
}
# legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2), J=2^6, p=15", side=3, adj=0, font=3)

## add a local image in the first plot ##
par(mar=c(3,3,1,1), font=3, fig = c(0.07, 0.31, 0.58, 0.82), new = T)  
plot(R2, log(MSPE_R2[,1], 10), ylim=c(-0.01, 0.08), type="n", #range(log(MSPE_R2, 10))[1]
      xlim=c(0.48, 0.52), cex=0.5, xlab="", ylab="")
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2[, i], 10), col=i+1, pch=i+1)
}
















