########################################################################
### compare six methods for DGP1 and fix N=2^13 , T=3 ############################
### consider two cases for candidate models: Case 2(group) ###########

source("repeat functions.R")


##### Case 2 (group) ###########
### vary R^2 via sigma0= sqrt((1-R2)/R2) #######
R2_vec = seq(0.1, 0.9, by=0.1)
MSPE_R2 = matrix(NA, length(R2_vec), 6) 
RISK_R2 = matrix(NA, length(R2_vec), 6) 
for(i in 1:length(R2_vec)){
  result2 = repfun1(num_rep=200, N=2^13, n_test=100,  J=2^6, p=15, p0=200, rho=0.5, 
                 sigma0=sqrt((1-R2_vec[i])/R2_vec[i]), T=3, cm_type="group", n_group=5)
  MSPE_R2[i, ] = result2$mspe_mean[1, ]
  RISK_R2[i, ] = result2$risk_mean[1, ] 
  cat("R2 = ", R2_vec[i], "\n")
}
colnames(MSPE_R2) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "gMAP", "DP")
colnames(RISK_R2) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "gMAP", "DP")
write.table(MSPE_R2, file="MSPE across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt")
write.table(RISK_R2, file="RISK across R2(N=2^13, J=2^6, p=15, T=3)(Case 2).txt")




########################################################################
### vary N  ##################################################
N_vec = (1:5)*2^12
R2=0.5

result1 = repfun1_vary_N(num_rep=200, N=N_vec, n_test=100, J=2^6, p=15, p0=200,  
              rho=0.5, sigma0=sqrt((1-R2)/R2), T=3, cm_type="group", n_group=5)
MSPE_N = result1$mspe_mean
RISK_N = result1$risk_mean

colnames(MSPE_N) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
colnames(RISK_N) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
write.table(MSPE_N, file="MSPE across N(J=2^6, p=15, R2=0.5, T=3)(Case 2).txt")
write.table(RISK_N, file="RISK across N(J=2^6, p=15, R2=0.5, T=3)(Case 2).txt")



### vary J ##############
J_vec = 2^(1:6)
R2=0.5

result3 = repfun1_vary_J(num_rep=200, N=2^13, n_test=100, 
              J=J_vec, p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=3, cm_type="group", n_group=5)
MSPE_J = result3$mspe_mean
RISK_J = result3$risk_mean

colnames(MSPE_J) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
colnames(RISK_J) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
write.table(MSPE_J, file="MSPE across J(N=2^13, p=15, R2=0.5, T=3)(Case 2).txt")
write.table(RISK_J, file="RISK across J(N=2^13, p=15, R2=0.5, T=3)(Case 2).txt")


##  vary p (dimension) ###########
p_vec = c(5, 10, 15, 20, 30)
R2=0.5

result4 = repfun1_vary_p(num_rep=200, N=2^13, n_test=100, 
              J=2^6, p=p_vec, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=3, cm_type="group", n_group=5)
MSPE_p = result4$mspe_mean
RISK_p = result4$risk_mean

colnames(MSPE_p) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
colnames(RISK_p) = c("gMAP", "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP")
write.table(MSPE_p, file="MSPE across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 2).txt")
write.table(RISK_p, file="RISK across p(N=2^13, J=2^6, R2=0.5, T=3)(Case 2).txt")















