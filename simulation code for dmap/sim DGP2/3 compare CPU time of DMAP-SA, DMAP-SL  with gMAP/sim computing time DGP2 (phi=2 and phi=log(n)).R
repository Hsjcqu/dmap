## we fix T=5  ######

source("repeat functions.R")


##### (1) CPU time across N for DMAP-SA, DMAP-SL, gMAP  ###
R2=0.5
N_vec = (1:5)*2^12  ##
t_cost1 = matrix(NA, length(N_vec), 6)
for(i in 1:length(N_vec)){
  result1 = rep_compare_time_fun2(num_rep=20, N=N_vec[i], n_test=100, 
             J=2^6, p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
  t_cost1[i, ] = result1$time_mean[1, ] # result: 
  cat("N = ", N_vec[i], "\n")
}

colnames(t_cost1) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost1, file="CPU time across N(J=2^6, R2=0.5, p=15)DGP2.txt")


##### (2) CPU time across J for DMAP-SA, DMAP-SL, gMAP  ###
R2=0.5
J_vec=2^(2:7)
t_cost2 = matrix(NA, length(J_vec), 6)
for(i in 1:length(J_vec)){
  result2 = rep_compare_time_fun2(num_rep=20, N=2^13, n_test=100, 
             J=J_vec[i], p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
  t_cost2[i, ] = result2$time_mean[1, ] # result: 
  cat("J = ", J_vec[i], "\n")
}

colnames(t_cost2) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost2, file="CPU time across J(N=2^13, R2=0.5, p=15)DGP2.txt")



##### (3) CPU time  across p for DMAP-SA, DMAP-SL, gMAP ###
R2=0.5
p_vec=c(5, 10, 15, 20, 30)
t_cost3 = matrix(NA, length(p_vec), 6)
for(i in 1:length(p_vec)){
  result3 = rep_compare_time_fun2(num_rep=20, N=2^13, n_test=100, 
             J=2^6, p=p_vec[i], p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
  t_cost3[i, ] = result3$time_mean[1, ] # result: 
  cat("p = ", p_vec[i], "\n")
}

colnames(t_cost3) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost3, file="CPU time across p(N=2^13, J=2^6, R2=0.5)DGP2.txt")






###############################################################################################
###### new  #################################################################################
source("repeat functions.R")

##### (1) CPU time across N for DMAP-SA, DMAP-SL, gMAP  ###
R2=0.5
N_vec = 2^(12:16)  ##

result1 = repfun2_compare_time_N(num_rep=200, N=N_vec, n_test=100, 
             J=2^6, p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
t_cost1 = result1$time_mean # result: 


colnames(t_cost1) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost1, file="CPU time across N(J=2^6, R2=0.5, p=15)DGP2.txt")


##### (2) CPU time across J for DMAP-SA, DMAP-SL, gMAP  ###
R2=0.5
J_vec=2^(2:7)

result2 = repfun2_compare_time_J(num_rep=200, N=2^13, n_test=100, 
             J=J_vec, p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
t_cost2 = result2$time_mean # result: 


colnames(t_cost2) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost2, file="CPU time across J(N=2^13, R2=0.5, p=15)DGP2.txt")



##### (3) CPU time  across p for DMAP-SA, DMAP-SL, gMAP ###
R2=0.5
p_vec=c(5, 10, 15, 20, 30)

result3 = repfun2_compare_time_p(num_rep=200, N=2^13, n_test=100, 
             J=2^6, p=p_vec, p0=200, rho=0.5, sigma0=sqrt((1-R2)/R2), T=5)
t_cost3 = result3$time_mean # result: 


colnames(t_cost3) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost3, file="CPU time across p(N=2^13, J=2^6, R2=0.5)DGP2.txt")


##### (4) CPU time across R^2 for DMAP-SA, DMAP-SL, gMAP  ###
R2_vec=c(0.1, 0.3, 0.5, 0.7, 0.9)

t_cost4 = matrix(NA, length(R2_vec), 6)
for(i in 1:length(R2_vec)){
  result4 = rep_compare_time_fun2(num_rep=200, N=2^13, n_test=100, 
             J=2^6, p=15, p0=200, rho=0.5, sigma0=sqrt((1-R2_vec[i])/R2_vec[i]), T=5)
  t_cost4[i, ] = result4$time_mean[1, ] # result: 
  cat("R^2 = ", R2_vec[i], "\n")
}

colnames(t_cost4) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
write.table(t_cost4, file="CPU time across R2(N=2^13, J=2^6, p=15)DGP2.txt")








