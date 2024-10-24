####################################################################################
##### plot MSPE across N, J, p together ###################################



MSPE_N = read.table("MSPE across N(J=2^6, p=15, R2=0.5, T=5)(Case 1).txt", header=TRUE)
MSPE_Nb = read.table("MSPE across N(J=2^6, p=15, R2=0.5, T=5)(Case 2).txt", header=TRUE)
N_vec = (1:5)*2^12

MSPE_J = read.table("MSPE across J(N=2^13, p=15, R2=0.5, T=5)(Case 1).txt", header=TRUE)
MSPE_Jb = read.table("MSPE across J(N=2^13, p=15, R2=0.5, T=5)(Case 2).txt", header=TRUE)
J_vec = 2^(1:6)

MSPE_p = read.table("MSPE across p(N=2^13, J=2^6, R2=0.5, T=5)(Case 1).txt", header=TRUE)
MSPE_pb = read.table("MSPE across p(N=2^13, J=2^6, R2=0.5, T=5)(Case 2).txt", header=TRUE)
p_vec = c(5, 10, 15, 20, 30)

name_method = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")


windows(height=9, width=6, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(3, 2), mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   

## Case 1, MSPE_N ##
plot(N_vec, log(MSPE_N[,1], 10), ylim=c(0, 0.5), type="n", xaxt="n",  #range(log(MSPE_N, 10))[1]
      xlab="Sample size N", ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(N_vec, log(MSPE_N[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(N_vec, log(MSPE_N[, i], 10), col=i+1, pch=i+1)
}
#xlabepr = c(expression(2^12), expression(2^13), expression(2^14), expression(2^15), expression(2^16))
axis(1, at=N_vec, labels=N_vec)
legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1), R^2=0.5, J=2^6, p=15", side=3, adj=0, font=3)
# rect(0.47, 0.02, 0.53, 0.12, lty=1, border ="gray")  # same as xlim and ylim


## Case 2, MSPE_N ##
plot(N_vec, log(MSPE_Nb[,1], 10), ylim=c(0, 2.5), type="n", xaxt="n",  #range(log(MSPE_Nb, 10))[1]
      xlab="Sample size N", ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(N_vec, log(MSPE_Nb[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(N_vec, log(MSPE_Nb[, i], 10), col=i+1, pch=i+1)
}
# xlabepr = c(expression(2^12), expression(2^13), expression(2^14), expression(2^15), expression(2^16))
axis(1, at=N_vec, labels=N_vec)
# legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2), R^2=0.5, J=2^6, p=15", side=3, adj=0, font=3)


## Case 1, MSPE_J ##
plot(J_vec, MSPE_J[,1], ylim=c(1.1, 1.7), type="n", #range(MSPE_J)[1]
            xlab="Number of machines J", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(J_vec, MSPE_J[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, MSPE_J[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), N=2^13, p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 2, MSPE_J ##
plot(J_vec, MSPE_Jb[,1], ylim=c(1.1, 4), type="n",  #range(MSPE_Jb)[1]
            xlab="Number of machines J", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(J_vec, MSPE_Jb[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, MSPE_Jb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), N=2^13, p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 1, MSPE_p ##
plot(p_vec, MSPE_p[,1], ylim=c(1.0, 2.5), type="n",
            xlab="Number of covariates p", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(p_vec, MSPE_p[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, MSPE_p[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), N=2^13, J=2^6, R^2=0.5", side=3, adj=0, font=3)

## Case 2, MSPE_p ##
plot(p_vec, MSPE_pb[,1], ylim=c(1.0, 3), type="n",
            xlab="Number of covariates p", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(p_vec, MSPE_pb[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, MSPE_pb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), N=2^13, J=2^6, R^2=0.5", side=3, adj=0, font=3)









####################################################################################
##### plot Risk across N, J, p together ###################################

RISK_N = read.table("RISK across N(J=2^6, p=15, R2=0.5, T=5)(Case 1).txt", header=TRUE)
RISK_Nb = read.table("RISK across N(J=2^6, p=15, R2=0.5, T=5)(Case 2).txt", header=TRUE)
N_vec = (1:5)*2^12

RISK_J = read.table("RISK across J(N=2^13, p=15, R2=0.5, T=5)(Case 1).txt", header=TRUE)
RISK_Jb = read.table("RISK across J(N=2^13, p=15, R2=0.5, T=5)(Case 2).txt", header=TRUE)
J_vec = 2^(1:6)

RISK_p = read.table("RISK across p(N=2^13, J=2^6, R2=0.5, T=5)(Case 1).txt", header=TRUE)
RISK_pb = read.table("RISK across p(N=2^13, J=2^6, R2=0.5, T=5)(Case 2).txt", header=TRUE)
p_vec = c(5, 10, 15, 20, 30)

name_method = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")


windows(height=9, width=6, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(3, 2), mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   

## Case 1, RISK_N ##
plot(N_vec, RISK_N[,1], ylim=c(0, 1), type="n", xaxt="n",  #range(log(RISK_N, 10))[1]
      xlab="Sample size N", ylab="RISK")
for(i in 1:length(name_method)){
 lines(N_vec, RISK_N[, i], col=i+1, lty=i+1, lwd=1)
 points(N_vec, RISK_N[, i], col=i+1, pch=i+1)
}
# xlabepr = c(expression(2^12), expression(2^13), expression(2^14), expression(2^15), expression(2^16))
axis(1, at=N_vec, labels=N_vec)
legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1), R^2=0.5, J=2^6, p=15", side=3, adj=0, font=3)
# rect(0.47, 0.02, 0.53, 0.12, lty=1, border ="gray")  # same as xlim and ylim


## Case 2, RISK_N ##
plot(N_vec, RISK_Nb[,1], ylim=c(0, 4), type="n", xaxt="n",  #range(log(RISK_N, 10))[1]
      xlab="Sample size N", ylab="RISK")
for(i in 1:length(name_method)){
 lines(N_vec, RISK_Nb[, i], col=i+1, lty=i+1, lwd=1)
 points(N_vec, RISK_Nb[, i], col=i+1, pch=i+1)
}
# xlabepr = c(expression(2^12), expression(2^13), expression(2^14), expression(2^15), expression(2^16))
axis(1, at=N_vec, labels=N_vec)
# legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2), R^2=0.5, J=2^6, p=15", side=3, adj=0, font=3)


## Case 1, RISK_J ##
plot(J_vec, RISK_J[,1], ylim=c(0.1, 0.7), type="n", #range(RISK_J)[1]
            xlab="Number of machines J", ylab="RISK")
for(i in 1:length(name_method)){
  lines(J_vec, RISK_J[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, RISK_J[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), N=2^13, p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 2, RISK_J ##
plot(J_vec, RISK_Jb[,1], ylim=c(0, 3), type="n",  #range(RISK_Jb)[1]
            xlab="Number of machines J", ylab="RISK")
for(i in 1:length(name_method)){
  lines(J_vec, RISK_Jb[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, RISK_Jb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), N=2^13, p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 1, RISK_p ##
plot(p_vec, RISK_p[,1], ylim=c(0, 1.2), type="n",
            xlab="Number of covariates p", ylab="RISK")
for(i in 1:length(name_method)){
  lines(p_vec, RISK_p[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, RISK_p[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), N=2^13, J=2^6, R^2=0.5", side=3, adj=0, font=3)

## Case 2, RISK_p ##
plot(p_vec, RISK_pb[,1], ylim=c(0, 3), type="n",
            xlab="Number of covariates p", ylab="RISK")
for(i in 1:length(name_method)){
  lines(p_vec, RISK_pb[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, RISK_pb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), N=2^13, J=2^6, R^2=0.5", side=3, adj=0, font=3)






