## we fix T=3  #######################################################################################
source("repeat functions.R")

## N=2^12, R^2=0.5, J=2^5, p=15 ### 
result1 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^12, n_test=100, J=2^5, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=3, cm_type="nested", n_group=NULL)
re1 = t(cbind(result1$time_mean[1, ], result1$mspe_mean[1,]))
rownames(re1) = c("Time(s)", "MSPE"); re1

## N=2^13, R^2=0.5, J=2^5, p=15 ### 
result2 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^5, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=3, cm_type="nested", n_group=NULL)
re2 = t(cbind(result2$time_mean[1, ], result2$mspe_mean[1,]))
rownames(re2) = c("Time(s)", "MSPE"); re2

## N=2^13, R^2=0.5, J=2^6, p=15 ### 
result3 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=3, cm_type="nested", n_group=NULL)
re3 = t(cbind(result3$time_mean[1, ], result3$mspe_mean[1,]))
rownames(re3) = c("Time(s)", "MSPE"); re3

## N=2^13, R^2=0.9, J=2^6, p=15 ### 
result4 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T=3, cm_type="nested", n_group=NULL)
re4 = t(cbind(result4$time_mean[1, ], result4$mspe_mean[1,]))
rownames(re4) = c("Time(s)", "MSPE"); re4


## N=2^13, R^2=0.9, J=2^6, p=30 ### 
result5 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=30, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T=3, cm_type="nested", n_group=NULL)
re5 = t(cbind(result5$time_mean[1, ], result5$mspe_mean[1,]))
rownames(re5) = c("Time(s)", "MSPE"); re5

re = rbind(re1, re2, re3, re4, re5); re 
write.table(re, file="MSPE and time 5 settings(Case 1)(T=3).txt")


S1 = "(N, R^2, J, p)=(2^12, 0.5, 2^5, 15)"
S2 = "(N, R^2, J, p)=(2^13, 0.5, 2^5, 15)"
S3 = "(N, R^2, J, p)=(2^13, 0.5, 2^6, 15)"
S4 = "(N, R^2, J, p)=(2^13, 0.9, 2^6, 15）"
S5 = "(N, R^2, J, p)=(2^13, 0.9, 2^6, 30）"
setting = c(S1, S1, S2, S2, S3, S3, S4, S4, S5, S5)
cbind(setting, re)

# reA = read.table("MSPE and time 5 settings(Case 1)(T=3).txt", header=T)





## we fix T=1  #######################################################################################
source("repeat functions.R")

## N=2^12, R^2=0.5, J=2^5, p=15 ### 
result1 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^12, n_test=100, J=2^5, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=1, cm_type="nested", n_group=NULL)
re1 = t(cbind(result1$time_mean[1, ], result1$mspe_mean[1,]))
rownames(re1) = c("Time(s)", "MSPE"); re1

## N=2^13, R^2=0.5, J=2^5, p=15 ### 
result2 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^5, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=1, cm_type="nested", n_group=NULL)
re2 = t(cbind(result2$time_mean[1, ], result2$mspe_mean[1,]))
rownames(re2) = c("Time(s)", "MSPE"); re2

## N=2^13, R^2=0.5, J=2^6, p=15 ### 
result3 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.5)/0.5), T=1, cm_type="nested", n_group=NULL)
re3 = t(cbind(result3$time_mean[1, ], result3$mspe_mean[1,]))
rownames(re3) = c("Time(s)", "MSPE"); re3

## N=2^13, R^2=0.9, J=2^6, p=15 ### 
result4 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=15, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T=1, cm_type="nested", n_group=NULL)
re4 = t(cbind(result4$time_mean[1, ], result4$mspe_mean[1,]))
rownames(re4) = c("Time(s)", "MSPE"); re4


## N=2^13, R^2=0.9, J=2^6, p=30 ### 
result5 = rep_compare_MSPEandTime_fun1(num_rep=200, N=2^13, n_test=100, J=2^6, p=30, 
              p0=200, rho=0.5, sigma0=sqrt((1-0.9)/0.9), T=1, cm_type="nested", n_group=NULL)
re5 = t(cbind(result5$time_mean[1, ], result5$mspe_mean[1,]))
rownames(re5) = c("Time(s)", "MSPE"); re5

re = rbind(re1, re2, re3, re4, re5); re 
write.table(re, file="MSPE and time 5 settings(Case 1)(T=1).txt")


S1 = "(N, R^2, J, p)=(2^12, 0.5, 2^5, 15)"
S2 = "(N, R^2, J, p)=(2^13, 0.5, 2^5, 15)"
S3 = "(N, R^2, J, p)=(2^13, 0.5, 2^6, 15)"
S4 = "(N, R^2, J, p)=(2^13, 0.9, 2^6, 15）"
S5 = "(N, R^2, J, p)=(2^13, 0.9, 2^6, 30）"
setting = c(S1, S1, S2, S2, S3, S3, S4, S4, S5, S5)
cbind(setting, re)






