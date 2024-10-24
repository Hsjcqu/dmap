##############################################################################
#### plot comparison of phi=2 and phi=log(n) for DMAP-SA and DMAP-SA 
###  for Example 1 (i.e., DGPfun1) #########

#### plot Case 1a (nested, R2=0.5) ######
MSPE_T_compare = read.table("MSPE across T(phi=2 and phi=log(n))(Case 1a).txt", header=TRUE)
RISK_T_compare = read.table("RISK across T(phi=2 and phi=log(n))(Case 1a).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=4, width=4.5, points=9)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, MSPE_T_compare[,1], ylim=c(range(MSPE_T_compare)[1], 1.35), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)


windows(height=4, width=4.5, points=9)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, RISK_T_compare[,1], ylim=c(range(RISK_T_compare)[1], 0.35), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)





#### plot Case 1b (nested, R2=0.9) ######
MSPE_T_compare2 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 1b).txt", header=TRUE)
RISK_T_compare2 = read.table("RISK across T(phi=2 and phi=log(n))(Case 1b).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=4, width=4.5, points=9)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, MSPE_T_compare2[,1], ylim=c(range(MSPE_T_compare2)[1], 0.26), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare2[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare2[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)


windows(height=4, width=4.5, points=9)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, RISK_T_compare2[,1], ylim=c(range(RISK_T_compare2)[1], 0.15), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare2[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare2[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)





#### plot Case 2a (group, R2=0.5) ######
MSPE_T_compare3 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 2a).txt", header=TRUE)
RISK_T_compare3 = read.table("RISK across T(phi=2 and phi=log(n))(Case 2a).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=4, width=4.5, points=9)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, MSPE_T_compare3[,1], ylim=c(range(MSPE_T_compare3)[1], 1.425), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare3[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare3[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)


windows(height=4, width=4.5, points=9)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, RISK_T_compare3[,1], ylim=c(range(RISK_T_compare3)[1], 0.425), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare3[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare3[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)




#### plot Case 2b (group, R2=0.9) ######
MSPE_T_compare4 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 2b).txt", header=TRUE)
RISK_T_compare4 = read.table("RISK across T(phi=2 and phi=log(n))(Case 2b).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=4, width=4.5, points=9)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, MSPE_T_compare4[,1], ylim=c(range(MSPE_T_compare4)[1], 0.515), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare4[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare4[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)


windows(height=4, width=4.5, points=9)
par(mar=c(4,4,2,2), font=3, bg="white", font.lab=3, font.axis=5)   
plot(T_vec, RISK_T_compare4[,1], ylim=c(range(RISK_T_compare4)[1], 0.41), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare4[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare4[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)



############################################################################################
### plot MSPE for Case 1a-2b together #####################################################

MSPE_T_compare = read.table("MSPE across T(phi=2 and phi=log(n))(Case 1a).txt", header=TRUE)
MSPE_T_compare2 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 1b).txt", header=TRUE)
MSPE_T_compare3 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 2a).txt", header=TRUE)
MSPE_T_compare4 = read.table("MSPE across T(phi=2 and phi=log(n))(Case 2b).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=8, width=8, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(2,2), mar=c(4,4,2,1.5), font=3, bg="white", font.lab=3, font.axis=5) 
  
plot(T_vec, MSPE_T_compare[,1], ylim=c(range(MSPE_T_compare)[1], 1.3), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1) and R^2=0.5", side=3, adj=0, font=3)

plot(T_vec, MSPE_T_compare2[,1], ylim=c(range(MSPE_T_compare2)[1], 0.2), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare2[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare2[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1) and R^2=0.9", side=3, adj=0, font=3)

plot(T_vec, MSPE_T_compare3[,1], ylim=c(range(MSPE_T_compare3)[1], 1.2), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare3[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare3[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2) and R^2=0.5", side=3, adj=0, font=3)

plot(T_vec, MSPE_T_compare4[,1], ylim=c(range(MSPE_T_compare4)[1], 0.303), type="n",
      xlab="Number of iterations T", ylab="MSPE")
for(i in 1:length(name_method2)){
 lines(T_vec, MSPE_T_compare4[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, MSPE_T_compare4[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2) and R^2=0.9", side=3, adj=0, font=3)



### plot Risk for Case 1a-2b together #####################################################

RISK_T_compare = read.table("RISK across T(phi=2 and phi=log(n))(Case 1a).txt", header=TRUE)
RISK_T_compare2 = read.table("RISK across T(phi=2 and phi=log(n))(Case 1b).txt", header=TRUE)
RISK_T_compare3 = read.table("RISK across T(phi=2 and phi=log(n))(Case 2a).txt", header=TRUE)
RISK_T_compare4 = read.table("RISK across T(phi=2 and phi=log(n))(Case 2b).txt", header=TRUE)

T_vec = 1:8
name_method2 = c(expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )')),
                  expression(paste('gMAP( ', phi,'=log(N) )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )')))

windows(height=8, width=8, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(2,2), mar=c(4,4,2,1.5), font=3, bg="white", font.lab=3, font.axis=5)   

plot(T_vec, RISK_T_compare[,1], ylim=c(range(RISK_T_compare)[1], 0.3), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare[, i], col=i+1, pch=i+1)
}
legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1) and R^2=0.5", side=3, adj=0, font=3)

plot(T_vec, RISK_T_compare2[,1], ylim=c(range(RISK_T_compare2)[1], 0.085), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare2[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare2[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 1) and R^2=0.9", side=3, adj=0, font=3)

plot(T_vec, RISK_T_compare3[,1], ylim=c(range(RISK_T_compare3)[1], 0.23), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare3[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare3[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2) and R^2=0.5", side=3, adj=0, font=3)

plot(T_vec, RISK_T_compare4[,1], ylim=c(range(RISK_T_compare4)[1], 0.19), type="n",
      xlab="Number of iterations T", ylab="Risk")
for(i in 1:length(name_method2)){
 lines(T_vec, RISK_T_compare4[, i], col=i+1, lty=i+1, lwd=1)
 points(T_vec, RISK_T_compare4[, i], col=i+1, pch=i+1)
}
#legend("topleft", name_method2, col=2:7, lty=2:7, pch=2:7, bty="n", lwd=1)
mtext("(Case 2) and R^2=0.9", side=3, adj=0, font=3)


### it seems T=2 or 3 is good choice for DMAP-SL ###







