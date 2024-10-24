############################################################################
### compare phi=2 and phi=log(n) for DGP1 ########################################
### fix N=2^13, J=2^6,  p=15, MSPE and Risk across T #########################

source("repeat functions.R")



### Case 1a: nested candidate models, R2 = 0.5  ###
T_vec = 1:8
n_T = length(T_vec)
result_compare1 = repfun1_compare_acrossT(num_rep=200, N=2^13, n_test=100, J=2^6, p=15,
                             p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T_vec=T_vec, cm_type="nested", n_group=NULL)

MSPE_T_compare1 = matrix(NA, length(T_vec), 6) 
RISK_T_compare1 = matrix(NA, length(T_vec), 6) 
for(i in 1:n_T){
  MSPE_T_compare1[i, ] = result_compare1$mspe_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  RISK_T_compare1[i, ] =  result_compare1$risk_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  cat("T = ", T_vec[i], "\n")
}
colnames(MSPE_T_compare1) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
colnames(RISK_T_compare1) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
write.table(MSPE_T_compare1, file="MSPE across T(phi=2 and phi=log(n))(Case 1a).txt")
write.table(RISK_T_compare1, file="RISK across T(phi=2 and phi=log(n))(Case 1a).txt")



### Case 1b: nested candidate models, R2 = 0.9  ###
T_vec = 1:8
n_T = length(T_vec)
result_compare2 = repfun1_compare_acrossT(num_rep=200, N=2^13, n_test=100, J=2^6, p=15,
                             p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T_vec=T_vec, cm_type="nested", n_group=NULL)

MSPE_T_compare2 = matrix(NA, length(T_vec), 6) 
RISK_T_compare2 = matrix(NA, length(T_vec), 6) 
for(i in 1:n_T){
  MSPE_T_compare2[i, ] = result_compare2$mspe_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  RISK_T_compare2[i, ] =  result_compare2$risk_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  cat("T = ", T_vec[i], "\n")
}
colnames(MSPE_T_compare2) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
colnames(RISK_T_compare2) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
write.table(MSPE_T_compare2, file="MSPE across T(phi=2 and phi=log(n))(Case 1b).txt")
write.table(RISK_T_compare2, file="RISK across T(phi=2 and phi=log(n))(Case 1b).txt")





### Case 2a: grouped candidate models, R2 = 0.5  ###
T_vec = 1:8
n_T = length(T_vec)
result_compare3 = repfun1_compare_acrossT(num_rep=200, N=2^13, n_test=100, J=2^6, p=15,
                             p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T_vec=T_vec, cm_type="group", n_group=5)

MSPE_T_compare3 = matrix(NA, length(T_vec), 6) 
RISK_T_compare3 = matrix(NA, length(T_vec), 6) 
for(i in 1:n_T){
  MSPE_T_compare3[i, ] = result_compare3$mspe_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  RISK_T_compare3[i, ] =  result_compare3$risk_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  cat("T = ", T_vec[i], "\n")
}
colnames(MSPE_T_compare3) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
colnames(RISK_T_compare3) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
write.table(MSPE_T_compare3, file="MSPE across T(phi=2 and phi=log(n))(Case 2a).txt")
write.table(RISK_T_compare3, file="RISK across T(phi=2 and phi=log(n))(Case 2a).txt")



### Case 2b: grouped candidate models, R2 = 0.9  ###
T_vec = 1:8
n_T = length(T_vec)
result_compare4 = repfun1_compare_acrossT(num_rep=200, N=2^13, n_test=100, J=2^6, p=15,
                             p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T_vec=T_vec, cm_type="group", n_group=5)

MSPE_T_compare4 = matrix(NA, length(T_vec), 6) 
RISK_T_compare4 = matrix(NA, length(T_vec), 6) 
for(i in 1:n_T){
  MSPE_T_compare4[i, ] = result_compare4$mspe_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  RISK_T_compare4[i, ] =  result_compare4$risk_mean[1, c(1, 2, 2+i, 3+n_T, 4+n_T, 4+length(T_vec)+i)]
  cat("T = ", T_vec[i], "\n")
}
colnames(MSPE_T_compare4) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
colnames(RISK_T_compare4) = c("gMAP(phi=2)", "DMAP-SA(phi=2)", "DMAP-SL(phi=2)", 
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))")
write.table(MSPE_T_compare4, file="MSPE across T(phi=2 and phi=log(n))(Case 2b).txt")
write.table(RISK_T_compare4, file="RISK across T(phi=2 and phi=log(n))(Case 2b).txt")


