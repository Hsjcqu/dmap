################################################################################
################ plot MSPE or RISK across R2, J, p, respectively, ###########
#############  N=2^13 for DGP1, fix T=3  ##########

####### plot MSPE across R2, J, p together #############################

MSPE_R2 = read.table("MSPE across R2(N=2^13, J=2^6, p=15, T=3)(Case 1).txt", header=TRUE)
MSPE_R2b = read.table("MSPE across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt", header=TRUE)
R2 = seq(0.1, 0.9, by=0.1)

MSPE_J = read.table("MSPE across J(N=2^13, p=15, R2=0.5, T=3)(Case 1).txt", header=TRUE)
MSPE_Jb = read.table("MSPE across J(N=2^13, p=15, R2=0.5, T=3)(Case 2).txt", header=TRUE)
J_vec = 2^(1:6)

MSPE_p = read.table("MSPE across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 1).txt", header=TRUE)
MSPE_pb = read.table("MSPE across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 2).txt", header=TRUE)
p_vec = c(5, 10, 15, 20, 30)

name_method = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "gMAP", "DP")


windows(height=10, width=7, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(3, 2), mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   

## Case 1, MSPE_R2 ##
plot(R2, log(MSPE_R2[,1], 10), ylim=c(range(log(MSPE_R2, 10))[1], 0.4), type="n", #
      xlab=expression(R^2), ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2[, i], 10), col=i+1, pch=i+1)
}
legend("topright", name_method, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1), J=2^6, p=15", side=3, adj=0, font=3)
rect(0.47, 0.02, 0.53, 0.12, lty=1, border ="gray")  # same as xlim and ylim


## Case 2, MSPE_R2 ##
plot(R2, log(MSPE_R2b[,1], 10), ylim=c(range(log(MSPE_R2b, 10))[1], 1), type="n",
      xlab=expression(R^2), ylab=expression(log[10](MSPE)))
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2b[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2b[, i], 10), col=i+1, pch=i+1)
}
mtext("(Case 2), J=2^6, p=15", side=3, adj=0, font=3)

## Case 1, MSPE_J ##
plot(J_vec, MSPE_J[,1], ylim=c(1, 1.25), type="n", #range(MSPE_J)[1]
            xlab="Number of machines J", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(J_vec, MSPE_J[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, MSPE_J[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 2, MSPE_J ##
plot(J_vec, MSPE_Jb[,1], ylim=c(1, 2.2), type="n",  #range(MSPE_Jb)[1]
            xlab="Number of machines J", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(J_vec, MSPE_Jb[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, MSPE_Jb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 1, MSPE_p ##
plot(p_vec, MSPE_p[,1], ylim=c(1.0, 1.4), type="n",
            xlab="Number of covariates p", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(p_vec, MSPE_p[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, MSPE_p[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), J=2^6, R^2=0.5", side=3, adj=0, font=3)

## Case 2, MSPE_p ##
plot(p_vec, MSPE_pb[,1], ylim=c(1.0, 2.5), type="n",
            xlab="Number of covariates p", ylab="MSPE")
for(i in 1:length(name_method)){
  lines(p_vec, MSPE_pb[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, MSPE_pb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), J=2^6, R^2=0.5", side=3, adj=0, font=3)

## add a local image in the first plot ##
par(mar=c(3,3,1,1), font=3, fig = c(0.06, 0.28, 0.71, 0.88), new = T)  
plot(R2, log(MSPE_R2[,1], 10), ylim=c(0.02, 0.12), type="n", #range(log(MSPE_R2, 10))[1]
      xlim=c(0.47, 0.53), cex=0.5, xlab="", ylab="")
for(i in 1:length(name_method)){
 lines(R2, log(MSPE_R2[, i], 10), col=i+1, lty=i+1, lwd=1)
 points(R2, log(MSPE_R2[, i], 10), col=i+1, pch=i+1)
}






####### plot Risk across R2, J, p together ##############################

RISK_R2 = read.table("RISK across R2(N=2^13, J=2^6, p=15, T=3)(Case 1).txt", header=TRUE)
RISK_R2b = read.table("RISK across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt", header=TRUE)
R2 = seq(0.1, 0.9, by=0.1)

RISK_J = read.table("RISK across J(N=2^13, p=15, R2=0.5, T=3)(Case 1).txt", header=TRUE)
RISK_Jb = read.table("RISK across J(N=2^13, p=15, R2=0.5, T=3)(Case 2).txt", header=TRUE)
J_vec = 2^(1:6)

RISK_p = read.table("RISK across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 1).txt", header=TRUE)
RISK_pb = read.table("RISK across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 2).txt", header=TRUE)
p_vec = c(5, 10, 15, 20, 30)

name_method = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "gMAP", "DP")

windows(height=10, width=7, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(3, 2), mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   

## Case 1, Risk_R2 ##
plot(R2, RISK_R2[,1], ylim=c(0, 0.6), type="n", #range(RISK_R2)[1]
      xlab=expression(R^2), ylab="RISK")
for(i in 1:length(name_method)){
 lines(R2, RISK_R2[, i], col=i+1, lty=i+1, lwd=1)
 points(R2, RISK_R2[, i], col=i+1, pch=i+1)
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

## Case 1, Risk_J ##
plot(J_vec, RISK_J[,1], ylim=c(0, 0.3), type="n",  #range(RISK_J)[1]
            xlab="Number of machines J", ylab="RISK")
for(i in 1:length(name_method)){
  lines(J_vec, RISK_J[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, RISK_J[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 2, Risk_J ##
plot(J_vec, RISK_Jb[,1], ylim=c(0, 1), type="n",  #range(RISK_Jb)[1]
            xlab="Number of machines J", ylab="RISK")
for(i in 1:length(name_method)){
  lines(J_vec, RISK_Jb[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, RISK_Jb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), p=15, R^2=0.5", side=3, adj=0, font=3)

## Case 1, Risk_p ##
plot(p_vec, RISK_p[,1], ylim=c(0, 0.4), type="n",  # range(RISK_p)[1]
            xlab="Number of covariates p", ylab="RISK")
for(i in 1:length(name_method)){
  lines(p_vec, RISK_p[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, RISK_p[, i], col=i+1, pch=i+1)
}
mtext("(Case 1), J=2^6, R^2=0.5", side=3, adj=0, font=3)

## Case 2, Risk_p ##
plot(p_vec, RISK_pb[,1], ylim=c(0, 1), type="n",   #range(RISK_pb)[1]
            xlab="Number of covariates p", ylab="RISK")
for(i in 1:length(name_method)){
  lines(p_vec, RISK_pb[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, RISK_pb[, i], col=i+1, pch=i+1)
}
mtext("(Case 2), J=2^6, R^2=0.5", side=3, adj=0, font=3)


