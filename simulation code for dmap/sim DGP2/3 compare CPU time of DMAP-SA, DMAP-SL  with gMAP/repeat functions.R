
source("example functions.R")
source("main functions.R")

#####################################################################################
### compare six methods across N R^2, J, p ############################################
## Example 1 ###
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



repfun1_vary_N = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_N = length(N)
  res = ris = matrix(NA, num_rep, 6*n_N)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=max(N), n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    Y = dat1$y.train
    X = dat1$x.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(N)){
       mu_train = dat1$mu_train[1:N[j]]
       Xj = as.matrix(X[1:N[j], ]) 
       Yj = as.vector(Y[1:N[j]])
       re1 = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T)
       res[s, j] = re1$mspe
       res[s, n_N+j ] = re2$mspe
       res[s, 2*n_N+j] = re3$mspe
       res[s, 3*n_N+j] = re2$mspe_ew  
       res[s, 4*n_N+j] = re3$mspe_ew
       res[s, 5*n_N+j] = re4$mspe

       ris[s, j] = mean((mu_train-re1$mu_w_est)^2)
       ris[s, n_N+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 2*n_N+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 3*n_N+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 4*n_N+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 5*n_N+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
   
  mspe_mean_N = matrix(mspe_mean[1, ],  nrow=n_N)
  risk_mean_N =  matrix(risk_mean[1, ],  nrow=n_N)
  mspe_median_N = matrix(mspe_median[1, ],  nrow=n_N)
  risk_median_N = matrix(risk_median[1, ],  nrow=n_N)

  colnames(mspe_mean_N) =  colnames(risk_mean_N) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_N) =  colnames(risk_median_N) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_N, mspe_median=mspe_median_N,
       risk_mean=risk_mean_N, risk_median=risk_median_N)
}




repfun1_vary_J = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_J = length(J)
  res = ris = matrix(NA, num_rep, 1+5*n_J)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train
    re1 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)

    for(j in 1:length(J)){
       re2 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T)
       res[s, 1+j ] = re2$mspe
       res[s, 1+n_J+j] = re3$mspe
       res[s, 1+2*n_J+j] = re2$mspe_ew  
       res[s, 1+3*n_J+j] = re3$mspe_ew
       res[s, 1+4*n_J+j] = re4$mspe

       ris[s, 1+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 1+n_J+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 1+2*n_J+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 1+3*n_J+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 1+4*n_J+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd)) 
  
  mspe_mean_J = cbind(mspe_mean[1, 1],  matrix(mspe_mean[1, 2:(1+5*n_J)],  nrow=n_J))
  risk_mean_J = cbind(risk_mean[1, 1],  matrix(risk_mean[1, 2:(1+5*n_J)],  nrow=n_J))
  mspe_median_J =  cbind(mspe_median[1, 1],  matrix(mspe_median[1, 2:(1+5*n_J)],  nrow=n_J))
  risk_median_J = cbind(risk_median[1, 1],  matrix(risk_median[1, 2:(1+5*n_J)],  nrow=n_J))

  colnames(mspe_mean_J) =  colnames(risk_mean_J) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_J) =  colnames(risk_median_J) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_J, mspe_median=mspe_median_J,
       risk_mean=risk_mean_J, risk_median=risk_median_J)
}


repfun1_vary_p = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_p = length(p)
  res = ris = matrix(NA, num_rep, 6*n_p)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    Y = dat1$y.train
    X = dat1$x.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train
    for(j in 1:length(p)){
       Xj = as.matrix(X[, 1:p[j]]) 
       X_newj = as.matrix(X_new[, 1:p[j]]) 
       re1 = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T)
       res[s, j] = re1$mspe
       res[s, n_p+j ] = re2$mspe
       res[s, 2*n_p+j] = re3$mspe
       res[s, 3*n_p+j] = re2$mspe_ew  
       res[s, 4*n_p+j] = re3$mspe_ew
       res[s, 5*n_p+j] = re4$mspe

       ris[s, j] = mean((mu_train-re1$mu_w_est)^2)
       ris[s, n_p+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 2*n_p+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 3*n_p+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 4*n_p+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 5*n_p+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
   
  mspe_mean_p = matrix(mspe_mean[1, ],  nrow=n_p)
  risk_mean_p =  matrix(risk_mean[1, ],  nrow=n_p)
  mspe_median_p = matrix(mspe_median[1, ],  nrow=n_p)
  risk_median_p = matrix(risk_median[1, ],  nrow=n_p)

  colnames(mspe_mean_p) =  colnames(risk_mean_p) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_p) =  colnames(risk_median_p) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_p, mspe_median=mspe_median_p,
       risk_mean=risk_mean_p, risk_median=risk_median_p)
}

## Example 2 ###
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


repfun2_vary_N = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_N = length(N)
  res = ris = matrix(NA, num_rep, 6*n_N)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=max(N), n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    Y = dat1$y.train
    X = dat1$x.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(N)){
       mu_train = dat1$mu_train[1:N[j]]
       Xj = as.matrix(X[1:N[j], ]) 
       Yj = as.vector(Y[1:N[j]])
       re1 = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T)
       res[s, j] = re1$mspe
       res[s, n_N+j ] = re2$mspe
       res[s, 2*n_N+j] = re3$mspe
       res[s, 3*n_N+j] = re2$mspe_ew  
       res[s, 4*n_N+j] = re3$mspe_ew
       res[s, 5*n_N+j] = re4$mspe

       ris[s, j] = mean((mu_train-re1$mu_w_est)^2)
       ris[s, n_N+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 2*n_N+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 3*n_N+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 4*n_N+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 5*n_N+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
   
  mspe_mean_N = matrix(mspe_mean[1, ],  nrow=n_N)
  risk_mean_N =  matrix(risk_mean[1, ],  nrow=n_N)
  mspe_median_N = matrix(mspe_median[1, ],  nrow=n_N)
  risk_median_N = matrix(risk_median[1, ],  nrow=n_N)

  colnames(mspe_mean_N) =  colnames(risk_mean_N) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_N) =  colnames(risk_median_N) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_N, mspe_median=mspe_median_N,
       risk_mean=risk_mean_N, risk_median=risk_median_N)
}




repfun2_vary_J = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_J = length(J)
  res = ris = matrix(NA, num_rep, 1+5*n_J)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train
    re1 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)

    for(j in 1:length(J)){
       re2 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T)
       res[s, 1+j ] = re2$mspe
       res[s, 1+n_J+j] = re3$mspe
       res[s, 1+2*n_J+j] = re2$mspe_ew  
       res[s, 1+3*n_J+j] = re3$mspe_ew
       res[s, 1+4*n_J+j] = re4$mspe

       ris[s, 1+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 1+n_J+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 1+2*n_J+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 1+3*n_J+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 1+4*n_J+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd)) 
  
  mspe_mean_J = cbind(mspe_mean[1, 1],  matrix(mspe_mean[1, 2:(1+5*n_J)],  nrow=n_J))
  risk_mean_J = cbind(risk_mean[1, 1],  matrix(risk_mean[1, 2:(1+5*n_J)],  nrow=n_J))
  mspe_median_J =  cbind(mspe_median[1, 1],  matrix(mspe_median[1, 2:(1+5*n_J)],  nrow=n_J))
  risk_median_J = cbind(risk_median[1, 1],  matrix(risk_median[1, 2:(1+5*n_J)],  nrow=n_J))

  colnames(mspe_mean_J) =  colnames(risk_mean_J) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_J) =  colnames(risk_median_J) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_J, mspe_median=mspe_median_J,
       risk_mean=risk_mean_J, risk_median=risk_median_J)
}


repfun2_vary_p = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_p = length(p)
  res = ris = matrix(NA, num_rep, 6*n_p)
  #colnames(res) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  # colnames(ris) = c("DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "MA-global", "DP-single")
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    Y = dat1$y.train
    X = dat1$x.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    mu_train = dat1$mu_train
    for(j in 1:length(p)){
       Xj = as.matrix(X[, 1:p[j]]) 
       X_newj = as.matrix(X_new[, 1:p[j]]) 
       re1 = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re4 = dp_single(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T)
       res[s, j] = re1$mspe
       res[s, n_p+j ] = re2$mspe
       res[s, 2*n_p+j] = re3$mspe
       res[s, 3*n_p+j] = re2$mspe_ew  
       res[s, 4*n_p+j] = re3$mspe_ew
       res[s, 5*n_p+j] = re4$mspe

       ris[s, j] = mean((mu_train-re1$mu_w_est)^2)
       ris[s, n_p+j ] = mean((mu_train-re2$mu_w_est)^2)
       ris[s, 2*n_p+j] = mean((mu_train-re3$mu_w_est)^2)
       ris[s, 3*n_p+j] = mean((mu_train-re2$mu_w_est_ew)^2)
       ris[s, 4*n_p+j] = mean((mu_train-re3$mu_w_est_ew)^2)  
       ris[s, 5*n_p+j] = mean((mu_train-re4$mu_est)^2)    
    }

    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
   
  mspe_mean_p = matrix(mspe_mean[1, ],  nrow=n_p)
  risk_mean_p =  matrix(risk_mean[1, ],  nrow=n_p)
  mspe_median_p = matrix(mspe_median[1, ],  nrow=n_p)
  risk_median_p = matrix(risk_median[1, ],  nrow=n_p)

  colnames(mspe_mean_p) =  colnames(risk_mean_p) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  colnames(mspe_median_p) =  colnames(risk_median_p) = c("MA-global",  "DMAP-SA", "DMAP-SL", "DMAP-SA-ew", "DMAP-SL-ew", "DP-single") 
  list(mspe_mean=mspe_mean_p, mspe_median=mspe_median_p,
       risk_mean=risk_mean_p, risk_median=risk_median_p)
}








#####################################################################################
#### compare  phi=2 and phi=log(n) for DMAP-SL, DMAP-SA, gMAP methods #############

### Example 1 ###
repfun1_compare_acrossT = function(num_rep, N, n_test, J, p, p0, rho, sigma0, 
             T_vec, cm_type="nested", n_group=NULL){
  n_T = length(T_vec)
  res = ris = matrix(NA, num_rep, 4+2*n_T)
  name_method = c("gMAP(phi=2)", "DMAP-SA(phi=2)", paste("DMAP-SL(phi=2)T=", T_vec, sep=""),  
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", paste("DMAP-SL(phi=log(n))T=", T_vec, sep=""))
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

    re1 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type=NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    re2 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type=NULL, cm_type=cm_type, n_group=n_group)
    res[s, 2] = re2$mspe
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    for(i in 1:n_T){
      re3 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T_vec[i], phi_type=NULL, cm_type=cm_type, n_group=n_group)
      res[s, 2+i] = re3$mspe
      ris[s, 2+i] = mean((mu_train-re3$mu_w_est)^2)
    }

    re1b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type=1, cm_type=cm_type, n_group=n_group)
    res[s, 2+n_T+1] = re1b$mspe
    ris[s, 2+n_T+1] = mean((mu_train-re1b$mu_w_est)^2)
    re2b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type=1, cm_type=cm_type, n_group=n_group)
    res[s, 2+n_T+2] = re2b$mspe
    ris[s, 2+n_T+2] = mean((mu_train-re2b$mu_w_est)^2)
    for(j in 1:n_T){
      re3b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T_vec[j], phi_type=1, cm_type=cm_type, n_group=n_group)
      res[s, 4+n_T+j] = re3b$mspe
      ris[s, 4+n_T+j] = mean((mu_train-re3b$mu_w_est)^2)
    }

    cat("Repeat s = ", s, "\n") 
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}

### Example 2 ###
repfun2_compare_acrossT = function(num_rep, N, n_test, J, p, p0, rho, sigma0, 
             T_vec, cm_type="nested", n_group=NULL){
  n_T = length(T_vec)
  res = ris = matrix(NA, num_rep, 4+2*n_T)
  name_method = c("gMAP(phi=2)", "DMAP-SA(phi=2)", paste("DMAP-SL(phi=2)T=", T_vec, sep=""),  
               "gMAP(phi=log(n))", "DMAP-SA(phi=log(n))", paste("DMAP-SL(phi=log(n))T=", T_vec, sep=""))
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

    re1 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type=NULL, cm_type=cm_type, n_group=n_group)
    res[s, 1] = re1$mspe
    ris[s, 1] = mean((mu_train-re1$mu_w_est)^2)
    re2 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type=NULL, cm_type=cm_type, n_group=n_group)
    res[s, 2] = re2$mspe
    ris[s, 2] = mean((mu_train-re2$mu_w_est)^2)
    for(i in 1:n_T){
      re3 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T_vec[i], phi_type=NULL, cm_type=cm_type, n_group=n_group)
      res[s, 2+i] = re3$mspe
      ris[s, 2+i] = mean((mu_train-re3$mu_w_est)^2)
    }

    re1b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type=1, cm_type=cm_type, n_group=n_group)
    res[s, 2+n_T+1] = re1b$mspe
    ris[s, 2+n_T+1] = mean((mu_train-re1b$mu_w_est)^2)
    re2b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, phi_type=1, cm_type=cm_type, n_group=n_group)
    res[s, 2+n_T+2] = re2b$mspe
    ris[s, 2+n_T+2] = mean((mu_train-re2b$mu_w_est)^2)
    for(j in 1:n_T){
      re3b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J, T=T_vec[j], phi_type=1, cm_type=cm_type, n_group=n_group)
      res[s, 4+n_T+j] = re3b$mspe
      ris[s, 4+n_T+j] = mean((mu_train-re3b$mu_w_est)^2)
    }

    cat("Repeat s = ", s, "\n") 
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  mspe_mean = rbind(colMeans(res), apply(res, 2, sd))
  mspe_median = rbind(apply(res, 2, median), apply(res, 2, sd))
  risk_mean = rbind(colMeans(ris), apply(ris, 2, sd))
  risk_median = rbind(apply(ris, 2, median), apply(ris, 2, sd))
  list(mspe_mean=mspe_mean, mspe_median=mspe_median, risk_mean=risk_mean, risk_median=risk_median)
}




#####################################################################################
##### compare computing time  ##########################################################

### Example 1 ###
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


repfun1_compare_time_N = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_N = length(N) 
  res = matrix(NA, num_rep, 6*n_N)
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=max(N), n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(N)){
       # mu_train = dat1$mu_train[1:N[j]]
       Xj = as.matrix(X[1:N[j], ]) 
       Yj = as.vector(Y[1:N[j]])

       re1 = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_N+j] = re2$t_cost
       res[s, 2*n_N+j] = re3$t_cost

       re1b = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re3b = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 3*n_N+j] = re1b$t_cost
       res[s, 4*n_N+j] = re2b$t_cost
       res[s, 5*n_N+j] = re3b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_N = matrix(time_mean[1, ], nrow=n_N)
  time_median_N = matrix(time_median[1, ], nrow=n_N)
 
  colnames(time_mean_N) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_N) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_N, time_median=time_median_N)
}



repfun1_compare_time_J = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_J = length(J) 
  res = matrix(NA, num_rep, 2+4*n_J)
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test

    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 2*n_J+1] = re3$t_cost
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4*n_J+2] = re3b$t_cost

    for(j in 1:length(J)){
       re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_J+j] = re2$t_cost
       re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 2*n_J+1+j] = re1b$t_cost
       res[s, 3*n_J+1+j] = re2b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_J = cbind(cbind(matrix(time_mean[1,1:(2*n_J)], nrow=n_J), time_mean[1,2*n_J+1]),
                      cbind(matrix(time_mean[1,(2*n_J+2):(4*n_J+1)], nrow=n_J), time_mean[1,4*n_J+2]))

  time_median_J = cbind(cbind(matrix(time_median[1,1:(2*n_J)], nrow=n_J), time_median[1,2*n_J+1]),
                      cbind(matrix(time_median[1,(2*n_J+2):(4*n_J+1)], nrow=n_J), time_median[1,4*n_J+2]))

  colnames(time_mean_J) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_J) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_J, time_median=time_median_J)
}



repfun1_compare_time_p = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_p = length(p) 
  res = matrix(NA, num_rep, 6*n_p)
  s = 1
  repeat{
    dat1 = DGPfun1(n.train=N, n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(p)){
       Xj = as.matrix(X[, 1:p[j]]) 
       X_newj = as.matrix(X_new[, 1:p[j]]) 
       re1 = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_p+j] = re2$t_cost
       res[s, 2*n_p+j] = re3$t_cost

       re1b = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re3b = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 3*n_p+j] = re1b$t_cost
       res[s, 4*n_p+j] = re2b$t_cost
       res[s, 5*n_p+j] = re3b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_p = matrix(time_mean[1, ], nrow=n_p)
  time_median_p = matrix(time_median[1, ], nrow=n_p)
 
  colnames(time_mean_p) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_p) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_p, time_median=time_median_p)
}



### Example 2 ####
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


repfun2_compare_time_N = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_N = length(N) 
  res = matrix(NA, num_rep, 6*n_N)
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=max(N), n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(N)){
       # mu_train = dat1$mu_train[1:N[j]]
       Xj = as.matrix(X[1:N[j], ]) 
       Yj = as.vector(Y[1:N[j]])

       re1 = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_N+j] = re2$t_cost
       res[s, 2*n_N+j] = re3$t_cost

       re1b = dmap_sa(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re3b = ma_global(X=Xj, Y=Yj, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 3*n_N+j] = re1b$t_cost
       res[s, 4*n_N+j] = re2b$t_cost
       res[s, 5*n_N+j] = re3b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_N = matrix(time_mean[1, ], nrow=n_N)
  time_median_N = matrix(time_median[1, ], nrow=n_N)
 
  colnames(time_mean_N) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_N) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_N, time_median=time_median_N)
}



repfun2_compare_time_J = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_J = length(J) 
  res = matrix(NA, num_rep, 2+4*n_J)
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test

    re3 = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
    res[s, 2*n_J+1] = re3$t_cost
    re3b = ma_global(X=X, Y=Y, X_new=X_new, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
    res[s, 4*n_J+2] = re3b$t_cost

    for(j in 1:length(J)){
       re1 = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_J+j] = re2$t_cost
       re1b = dmap_sa(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=X, Y=Y, X_new=X_new, Y_new=Y_new, J=J[j], T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 2*n_J+1+j] = re1b$t_cost
       res[s, 3*n_J+1+j] = re2b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_J = cbind(cbind(matrix(time_mean[1,1:(2*n_J)], nrow=n_J), time_mean[1,2*n_J+1]),
                      cbind(matrix(time_mean[1,(2*n_J+2):(4*n_J+1)], nrow=n_J), time_mean[1,4*n_J+2]))

  time_median_J = cbind(cbind(matrix(time_median[1,1:(2*n_J)], nrow=n_J), time_median[1,2*n_J+1]),
                      cbind(matrix(time_median[1,(2*n_J+2):(4*n_J+1)], nrow=n_J), time_median[1,4*n_J+2]))

  colnames(time_mean_J) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_J) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_J, time_median=time_median_J)
}



repfun2_compare_time_p = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  n_p = length(p) 
  res = matrix(NA, num_rep, 6*n_p)
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=max(p), p0=p0, rho=rho, sigma=sigma0)
    X = dat1$x.train
    Y = dat1$y.train
    X_new = dat1$x.test
    Y_new = dat1$y.test
    for(j in 1:length(p)){
       Xj = as.matrix(X[, 1:p[j]]) 
       X_newj = as.matrix(X_new[, 1:p[j]]) 
       re1 = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re2 = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       re3 = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = NULL, cm_type=cm_type, n_group=n_group)
       res[s, j] = re1$t_cost
       res[s, n_p+j] = re2$t_cost
       res[s, 2*n_p+j] = re3$t_cost

       re1b = dmap_sa(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re2b = dmap_sl(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, J=J, T=T, phi_type = 1, cm_type=cm_type, n_group=n_group)
       re3b = ma_global(X=Xj, Y=Y, X_new=X_newj, Y_new=Y_new, phi_type = 1, cm_type=cm_type, n_group=n_group)
       res[s, 3*n_p+j] = re1b$t_cost
       res[s, 4*n_p+j] = re2b$t_cost
       res[s, 5*n_p+j] = re3b$t_cost
    }
    cat("Repeat s = ", s, "\n")
    if(s>=num_rep){ break }else{ s=s+1 }
  } 
  time_mean = rbind(colMeans(res), apply(res, 2, sd))
  time_median = rbind(apply(res, 2, median), apply(res, 2, sd))

  time_mean_p = matrix(time_mean[1, ], nrow=n_p)
  time_median_p = matrix(time_median[1, ], nrow=n_p)
 
  colnames(time_mean_p) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  colnames(time_median_p) = c("DMAP-SA(phi=2)", "DMAP-SL(phi=2)", "MA-global(phi=2)", 
                  "DMAP-SA(phi=log(n))", "DMAP-SL(phi=log(n))", "MA-global(phi=log(n))")
  list(time_mean=time_mean_p, time_median=time_median_p)
}









##################################################################################################
##### compare CPUtime and MSPE DMAP-SL and DMAP-SL2 (two initial estiamtors) #################################

### Example 1 ###
rep_compare_MSPEandTime_fun1 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res_time = matrix(NA, num_rep, 4)
  res_mspe = matrix(NA, num_rep, 4)
  name_method = c("DMAP-SL(phi=2)(init1)",  "DMAP-SL(phi=2)(init2)",  
                 "DMAP-SL(phi=log(n))(init1)", "DMAP-SL(phi=log(n))(init2)")
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


### Example 2 ###
rep_compare_MSPEandTime_fun2 = function(num_rep, N, n_test, J, p, p0, rho, sigma0, T, cm_type="nested", n_group=NULL){
  res_time = matrix(NA, num_rep, 4)
  res_mspe = matrix(NA, num_rep, 4)
  name_method = c("DMAP-SL(phi=2)(init1)",  "DMAP-SL(phi=2)(init2)",  
                 "DMAP-SL(phi=log(n))(init1)", "DMAP-SL(phi=log(n))(init2)")
  colnames(res_time) = colnames(res_mspe) = name_method
  s = 1
  repeat{
    dat1 = DGPfun2(n.train=N, n.test=n_test, p_max=p, p0=p0, rho=rho, sigma=sigma0)
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










