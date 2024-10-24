
source("example functions.R")
source("main functions.R")





repfun1 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = ris = matrix(NA, num_rep, 6)
  colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re4 = dp_single(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T)
    res[s, 1] = re1$mspe
    res[s, 2] = re2$mspe
    res[s, 3] = re1$mspe_ew
    res[s, 4] = re2$mspe_ew  
    res[s, 5] = re3$mspe
    res[s, 6] = re4$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    ris[s, 3] = mean((mu_train-re1$mu_w_est_ew)^2)
    ris[s, 4] = mean((mu_train-re2$mu_w_est_ew)^2)  
    ris[s, 5] = mean((mu_train-re3$mu_w_est)^2)
    ris[s, 6] = mean((mu_train-re4$mu_est)^2)

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  }  
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}




repfun2 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = ris = matrix(NA, num_rep, 6)
  colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re4 = dp_single(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T)
    res[s, 1] = re1$mspe
    res[s, 2] = re2$mspe
    res[s, 3] = re1$mspe_ew
    res[s, 4] = re2$mspe_ew  
    res[s, 5] = re3$mspe
    res[s, 6] = re4$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    ris[s, 3] = mean((mu_train-re1$mu_w_est_ew)^2)
    ris[s, 4] = mean((mu_train-re2$mu_w_est_ew)^2)  
    ris[s, 5] = mean((mu_train-re3$mu_w_est)^2)
    ris[s, 6] = mean((mu_train-re4$mu_est)^2)

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}




#### compare  phi=2 and phi=log(n) for DMAP-SL, DMAP-SA, gMAP methods #############

rep_compare_fun1 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = ris = matrix(NA, num_rep, 6)
  name_method = c("DMAP-SA", "DMAP-SL", "MA-global", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(res) = name_method
  colnames(ris) = name_method
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    res[s, 2] = re2$mspe
    res[s, 3] = re3$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    ris[s, 3] = mean((mu_train-re3$mu_w_est)^2)

    re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4] = re1b$mspe
    res[s, 5] = re2b$mspe
    res[s, 6] = re3b$mspe
    ris[s, 4] = mean((mu_train-re1b$mu_w_est)^2)
    ris[s, 5] = mean((mu_train-re2b$mu_w_est)^2)
    ris[s, 6] = mean((mu_train-re3b$mu_w_est)^2)
   
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}



rep_compare_fun2 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = ris = matrix(NA, num_rep, 6)
  name_method = c("DMAP-SA", "DMAP-SL", "MA-global", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(res) = name_method
  colnames(ris) = name_method
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    res[s, 2] = re2$mspe
    res[s, 3] = re3$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    ris[s, 3] = mean((mu_train-re3$mu_w_est)^2)

    re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4] = re1b$mspe
    res[s, 5] = re2b$mspe
    res[s, 6] = re3b$mspe
    ris[s, 4] = mean((mu_train-re1b$mu_w_est)^2)
    ris[s, 5] = mean((mu_train-re2b$mu_w_est)^2)
    ris[s, 6] = mean((mu_train-re3b$mu_w_est)^2)

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}






##### compare computing time  ##########################################################

rep_compare_time_fun1 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = matrix(NA, num_rep, 6)
  name_method = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(res) = name_method
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$t_cost
    res[s, 2] = re2$t_cost
    res[s, 3] = re3$t_cost

    re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4] = re1b$t_cost
    res[s, 5] = re2b$t_cost
    res[s, 6] = re3b$t_cost
   
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 

  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  list(time_mean=time_mean, time_median=time_median)
}




rep_compare_time_fun2 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res = matrix(NA, num_rep, 6)
  name_method = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
               "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(res) = name_method
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$t_cost
    res[s, 2] = re2$t_cost
    res[s, 3] = re3$t_cost

    re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4] = re1b$t_cost
    res[s, 5] = re2b$t_cost
    res[s, 6] = re3b$t_cost
   
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 

  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  list(time_mean=time_mean, time_median=time_median)
}




##### compare CPUtime and MSPE DMAP-SL and DMAP-SL2 (two initial estiamtors) #################################

rep_compare_MSPEandTime_fun1 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res_time = matrix(NA, num_rep, 4)
  res_mspe = matrix(NA, num_rep, 4)
  name_method = c("DMAP-SL1(phi=2)",  "DMAP-SL2(phi=2)",  "DMAP-SL2(phi=log(n))", "DMAP-SL2(phi=log(n))")
  colnames(res_time) = colnames(res_mspe) = name_method
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train

    re1 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    re2 = dmap_sl2(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res_time[s, 1] = re1$t_cost
    res_time[s, 2] = re2$t_cost
    res_mspe[s, 1] = re1$mspe
    res_mspe[s, 2] = re2$mspe

    re3 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    re4 = dmap_sl2(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res_time[s, 3] = re3$t_cost
    res_time[s, 4] = re4$t_cost
    res_mspe[s, 3] = re3$mspe
    res_mspe[s, 4] = re4$mspe
   
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 

  time_mean = rbind(colMeans(res_time), apply(res_time, 2, sd))
  time_median = rbind(apply(res_time, 2, median), apply(res_time, 2, sd))

  mspe_mean = rbind(colMeans(res_mspe), apply(res_mspe, 2, sd))
  mspe_median = rbind(apply(res_mspe, 2, median), apply(res_mspe, 2, sd))

  list(time_mean=time_mean, time_median=time_median, mspe_mean=mspe_mean, mspe_median=mspe_median)
}











