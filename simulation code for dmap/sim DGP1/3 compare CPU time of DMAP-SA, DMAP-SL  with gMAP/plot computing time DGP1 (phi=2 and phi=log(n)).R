## we fix T=3  ######

########## plot CPU time across N, J, p together #######################################

Time_N = read.table("CPU time across N(J=2^6, R2=0.5, p=15)DGP1.txt", header=TRUE)
N_vec = 2^(12:16)
Time_J = read.table("CPU time across J(N=2^13, R2=0.5, p=15)DGP1.txt", header=TRUE)
J_vec = 2^(2:7)
Time_p = read.table("CPU time across p(N=2^13, J=2^6, R2=0.5)DGP1.txt", header=TRUE)
p_vec=c(5, 10, 15, 20, 30)
Time_R2 = read.table("CPU time across R2(N=2^13, J=2^6, p=15)DGP1.txt", header=TRUE)
R2_vec=c(0.1, 0.3, 0.5, 0.7, 0.9)


name_method = c(expression(paste('DMAP-SA( ', italic(phi),'=2 )')), 
                 expression(paste('DMAP-SL( ', phi,'=2 )(T=3)')),
                 expression(paste('gMAP( ', phi,'=2 )')),
                 expression(paste('DMAP-SA( ', phi,'=log(n) )')),
                 expression(paste('DMAP-SL( ', phi,'=log(n) )(T=3)')),
                 expression(paste('gMAP( ', phi,'=log(N) )'))
)

windows(height=6, width=6, points=10)
# png(file="xifuweight.png", width=4, height=4, pointsize =9, units="in",res=400)
par(mfrow=c(2,2), mar=c(4,4,2,1.5), font=3, bg="white", font.lab=3, font.axis=5)   

plot(N_vec, Time_N[,1], ylim=c(range(Time_N)[1], 0.9), type="n", xaxt="n",  # axes=F, 
            xlab="Total sample size N", ylab="CPU time(s)")
for(i in 1:length(name_method)){
  lines(N_vec, Time_N[, i], col=i+1, lty=i+1, lwd=1)
  points(N_vec, Time_N[, i], col=i+1, pch=i, cex=1.1)
}
#legend("topleft", name_method, col=2:7, lty=2:7, pch=1:6, bty="n", lwd=1, text.font=4)
xlabepr = c(expression(2^12), expression(2^13), expression(2^14), expression(2^15), expression(2^16))
axis(1, at=N_vec, labels=xlabepr)
mtext("(J, R^2, p)=(2^6, 0.5, 15)", side = 3, adj=0, font=3)

plot(J_vec, Time_J[,1], ylim=c(range(Time_J)[1], 0.18), type="n", #xaxt="n",  # axes=F, 
            xlab="Numer of machines J", ylab="CPU time(s)")
for(i in 1:length(name_method)){
  lines(J_vec, Time_J[, i], col=i+1, lty=i+1, lwd=1)
  points(J_vec, Time_J[, i], col=i+1, pch=i, cex=1.1)
}
#legend("topleft", name_method, col=2:7, lty=2:7, pch=1:6, bty="n", lwd=1, text.font=4)
#xlabepr = c(expression(2^2), expression(2^3), expression(2^4), expression(2^5), expression(2^6), expression(2^7))
#axis(1, at=J_vec, labels=xlabepr)
mtext("(N, R^2, p)=(2^13, 0.5, 15)", side = 3, adj=0, font=3)

plot(p_vec, Time_p[,1], ylim=c(range(Time_p)[1], 0.5), type="n", xaxt="n",  # axes=F, 
            xlab="Numer of covariates p", ylab="CPU time(s)")
for(i in 1:length(name_method)){
  lines(p_vec, Time_p[, i], col=i+1, lty=i+1, lwd=1)
  points(p_vec, Time_p[, i], col=i+1, pch=i, cex=1.1)
}
#legend("topleft", name_method, col=2:7, lty=2:7, pch=1:6, bty="n", lwd=1, text.font=4)
axis(1, at=p_vec, labels=p_vec)
mtext("(N, J, R^2)=(2^13, 2^6, 0.5)", side = 3, adj=0, font=3)


plot(R2_vec, Time_R2[,1], ylim=c(range(Time_R2)[1], 0.5), type="n", xaxt="n",  # axes=F, 
            xlab=expression(R^2), ylab="CPU time(s)")
for(i in 1:length(name_method)){
  lines(R2_vec, Time_R2[, i], col=i+1, lty=i+1, lwd=1)
  points(R2_vec, Time_R2[, i], col=i+1, pch=i, cex=1.1)
}
legend("topleft", name_method, col=2:7, lty=2:7, pch=1:6, bty="n", lwd=1, text.font=4)
axis(1, at=R2_vec, labels=R2_vec)
mtext("(N, J, p)=(2^13, 2^6, 15)", side = 3, adj=0, font=3)
















