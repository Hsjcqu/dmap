
# library(MASS)
# library(splines)
# library(quadprog)
# install.packages("quadprog")


CLS = function(H, b){
  # min w'*H'*H*w - 2b'*w, s.t: w>=0, sum(w)=1
  p = length(b)
  Dmat = t(H)%*%H 
  if(qr(Dmat)$rank<p) Dmat = Dmat + diag(rep(1e-9, p))
  dvec = b     
  Amat = cbind(rep(1, p), diag(rep(1, p)))
  bvec = c(1, rep(0, p))
  sc = norm(Dmat, "2")
  fit.cls = quadprog::solve.QP(Dmat/sc, dvec/sc, Amat, bvec, meq=1, factorized=FALSE)
  list(solution=fit.cls$solution, solution.unc = fit.cls$unconstrained.solution)
}





index_partition_fun =function(N, J){
  n0 = as.integer(N/J)
  n = ifelse(N!=n0*J,  n0+1, n0)
  id_set_list = list()
  for(j in 1:J){
    if(j<J){ 
      id_set_list[[j]] = (n*(j-1)+1):(n*j) 
    }else{
      id_set_list[[J]] = (n*(j-1)+1):N
    }
  }
  list(id_list=id_set_list, n_list=n)
}



candset_fun = function(L, type="nested", n_group=NULL){
  set_list = list()
  num = 0
  if(type=="nested"){
    for(i in 1:L){
      num = num + 1
      set_list[[num]] = 1:i
    }
  }
  if(type=="full"){
    for(j in 1:L){
      Mj = combn(1:L,j) 
      mp = ncol(Mj)
      for(k in 1:mp){
        num = num + 1
        set_list[[num]] = Mj[, k]
      }
    }
  }
  if(type=="marginal"){
    for(i in 1:L){
      num = num + 1
      set_list[[num]] = i
    }
  } 

  if(type=="group"){
    if(is.null(n_group)){ M=5 }else{ M=n_group }
    p0 = as.integer(L/M)
    for(j in 1:M){
      num = j
      ## set_list[[num]] = ifelse(j<M, (p0*(j-1)+1):(p0*j) , (p0*(j-1)+1):L)
      if(j<M){ 
         set_list[[num]] = (p0*(j-1)+1):(p0*j) 
      }else{
         set_list[[num]] = (p0*(j-1)+1):L
      }
    }
  }
  list(id_list = set_list, n=num, n_list = sapply(set_list, length))
}




#### Proposed DMAP-SA method #####################
dmap_sa = function(X, Y, X_new, Y_new, J, phi_type = NULL, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  n_new = nrow(X_new)

  id_list = index_partition_fun(N=N, J=J)
  n = id_list[[2]]
  xid_list = id_list[[1]]
  if(is.null(phi_type)){ phi = 2 }else{ phi=log(n) }

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]
  k_vec = cset_list[[3]]+1

  ## compute on each machine  ##
  w_est = matrix(nrow=M, ncol=J)
  # beta_mat = matrix(nrow=N, ncol=M)
  mu_est = matrix(NA, N, M)
  mu_pred = matrix(NA, n_new*J, M)
   
  time_vec= rep(NA, J)
  for(j in 1:J){  
    t1 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j] 
    for(m in 1:M){
      cmm = cset[[m]]
      Xjm = as.matrix(X[xid_j, cmm])
      dfj = data.frame(y=Yj, x=Xjm)
      lmfit_j =  lm(y~., data=dfj) 
      # betajm = coef(lmfit_j)
      mu_est[xid_j, m] = as.vector(fitted(lmfit_j))
      # mu_est[xid_j, m] = as.vector(Xjm%*%as.vector(betajm))
      mu_pred[(n_new*(j-1)+1):(n_new*j), m] = as.vector(predict(lmfit_j, newdata=data.frame(x=X_new[, cmm])))
      # mu_pred[(n_new*(j-1)+1):(n_new*j), m] = as.matrix(X_new[, cmm])%*%betajm  
    }
    ## find weights ##
    muj_est = as.matrix(mu_est[xid_j, ])
    df_full = data.frame(y=Yj, x=X[xid_j,])
    sigmaj2_max = sigma(lm(y~., data=df_full))^2
    H0 = muj_est
    b0 = as.vector(t(Yj)%*%muj_est) - phi*sigmaj2_max*k_vec/2
    w_est[, j] = CLS(H=H0, b=b0)$solution
    t2 = proc.time()
    time_vec[j] = (t2-t1)[3]
  } 
  
  t3 = proc.time()
  w_vec = rowMeans(w_est)
  mu_w_est = mu_est%*%w_vec
  mu_w_pred = rowMeans(matrix(as.vector(mu_pred%*%w_vec), ncol=J))
  mspe = mean( (Y_new - mu_w_pred)^2 )
  t4 = proc.time()
  t_cost = (t4-t3)[3] + time_vec[1]

  w_ew = rep(1, M)/M      ## equal weights 
  mu_w_est_ew = mu_est%*%w_ew
  mu_w_pred_ew = rowMeans(matrix(as.vector(mu_pred%*%w_ew), ncol=J))
  mspe_ew = mean( (Y_new - mu_w_pred_ew)^2 )

  list(mspe = mspe, mu_w_est=mu_w_est, mspe_ew=mspe_ew, mu_w_est_ew=mu_w_est_ew, t_cost=t_cost)
} 



#### Proposed DMAP-SL method #####################

## (1) use the first local esimator as the initial estimator ####
beta_csl_ma = function(X, Y, J, T, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  
  id_list = index_partition_fun(N=N, J=J)
  xid_list = id_list[[1]]

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]

  ## compute initial estimator of beta on first machine  ##
  t1 = proc.time()
  beta_initial_list = list()
  xid_1 = xid_list[[1]]
  Y1 = Y[xid_1] 
  for(m in 1:M){
    X1m = cbind(1, as.matrix(X[xid_1, cset[[m]]]))
    df1 = data.frame(y=Y1, x=X1m)
    beta_initial_list[[m]] = coef( lm(y~.+0, data=df1))
  }
  
  Sigma_mj = z_mj = list()
  t2 = proc.time()
  t_cost1 = (t2-t1)[3]
   
  t_record = rep(NA, J)
  for(j in 1:J){
    t3 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j]
    for(m in 1:M){
      Xjm = cbind(1, as.matrix(X[xid_j, cset[[m]]]))
      Sigma_mj[[M*(j-1)+m]] = as.matrix( t(Xjm)%*%Xjm/nrow(Xjm) )
      z_mj[[M*(j-1)+m]] = t(Xjm)%*%Yj/nrow(Xjm)
    }
    t4 = proc.time()
    t_record[j] = (t4-t3)[3] 
  }
  t_cost2 = max(t_record)

  t5 = proc.time()
  Sigma_m = z_m = list()
  t6 = proc.time()
  t_cost3 = (t6-t5)[3]
 
  tim = 0
  for(m in 1:M){
     t_in1 = proc.time()
     Sigm = 0;  zm = 0
     t_in2 = proc.time()
     t_record2 = rep(NA, J)
     for(j in 1:J){
       t_in3 = proc.time()
       Sigm = Sigm + Sigma_mj[[m+M*(j-1)]]
       zm = zm + z_mj[[m+M*(j-1)]]
       t_in4 = proc.time()
       t_record2[j] = (t_in4 - t_in3)[3] 
     }
     t_in5 = proc.time()
     Sigma_m[[m]] = Sigm/J
     z_m[[m]] = zm/J 
     t_in6 = proc.time()
     tim = tim+ (t_in2-t_in1)[3] + max(t_record2) + (t_in6-t_in5)[3]
  }
  
  t_in7 = proc.time()
  beta_list =  list()
  for(m in 1:M){ 
    beta_list[[m]] = beta_initial_list[[m]]
    Sigma_1m_inverse = MASS::ginv(Sigma_mj[[m]])
    A = diag(rep(1, length(cset[[m]])+1)) - Sigma_1m_inverse%*%Sigma_m[[m]]
    B =  Sigma_1m_inverse%*%z_m[[m]]
    for(t in 1:T){         
      beta_list[[m]] = A%*%beta_list[[m]] +  B 
    }
  } 
  t_in8 = proc.time()
  t_cost = t_cost1 + t_cost2 + t_cost3 + tim + (t_in8 - t_in7)[3]

  list(beta_list=beta_list, t_cost=t_cost) 
} 




dmap_sl = function(X, Y, X_new, Y_new, J, T=1, phi_type = NULL, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  n_new = nrow(X_new)
  res_beta = beta_csl_ma(X=X, Y=Y, J=J, T=T, cm_type=cm_type, n_group=n_group)
  beta_cls_list =res_beta$beta_list
  t_cost1 = res_beta$t_cost   

  id_list = index_partition_fun(N=N, J=J)
  n = id_list[[2]]
  xid_list = id_list[[1]]
  if(is.null(phi_type)){ phi = 2 }else{ phi=log(n) }

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]
  k_vec = cset_list[[3]]+1

  ## compute on each machine  ##
  t1 = proc.time()
  w_est = matrix(nrow=M, ncol=J)
  mu_est = matrix(NA, N, M)
  mu_pred = matrix(NA, n_new*J, M)
  t2 = proc.time()
  t_record = rep(NA, J)
  for(j in 1:J){  
    t3 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j] 
    for(m in 1:M){
      cmm = cset[[m]]
      Xjm = cbind(1, as.matrix(X[xid_j, cmm]))
      mu_est[xid_j, m] = as.vector(Xjm%*%beta_cls_list[[m]])
      mu_pred[(n_new*(j-1)+1):(n_new*j), m] = cbind(1, as.matrix(X_new[, cmm]))%*%beta_cls_list[[m]]
    }
    ## find weights ##
    muj_est = as.matrix(mu_est[xid_j, ])
    # Xfull = X[xid_j, ]
    df_full = data.frame(y=Yj, x=X[xid_j, ])
    sigmaj2_max = sigma(lm(y~., data=df_full))^2
    H0 = muj_est
    b0 = as.vector(t(Yj)%*%muj_est) - phi*sigmaj2_max*k_vec/2
    w_est[, j] = CLS(H=H0, b=b0)$solution
    t4 = proc.time()
    t_record[j] = (t4-t3)[3]
  } 

  t5 = proc.time()
  w_vec = rowMeans(w_est)
  mu_w_est = mu_est%*%w_vec
  mu_w_pred = rowMeans(matrix(as.vector(mu_pred%*%w_vec), ncol=J))
  mspe = mean( (Y_new - mu_w_pred)^2 )
  t6 = proc.time()
  t_cost = t_cost1 + (t2-t1)[3] + max(t_record) + (t6-t5)[3]

  w_ew = rep(1, M)/M      ## equal weights 
  mu_w_est_ew = mu_est%*%w_ew
  mu_w_pred_ew = rowMeans(matrix(as.vector(mu_pred%*%w_ew), ncol=J))
  mspe_ew = mean( (Y_new - mu_w_pred_ew)^2 )

  list(mspe = mspe, mu_w_est=mu_w_est, mspe_ew=mspe_ew, mu_w_est_ew=mu_w_est_ew, t_cost=t_cost)
} 



## (2) use the average of local esimators as the initial estimator ####

beta_csl_ma2 = function(X, Y, J, T, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  
  id_list = index_partition_fun(N=N, J=J)
  xid_list = id_list[[1]]

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]

  ## compute initial estimator: the average of local estimators  ##
  t0a = proc.time()
  beta_initial_list = list()
  for(m in 1:M){ beta_initial_list[[m]] = 0  }
  t_1_vec = rep(NA, J)
  t0b = proc.time()

  for(j in 1:J){
    t_init1 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j] 
    Xj = X[xid_j, ]
    for(m in 1:M){
      Xjm = cbind(1, as.matrix(Xj[, cset[[m]]]))
      df1 = data.frame(y=Yj, x=Xjm)
      beta_m_init = as.vector(coef( lm(y~.+0, data=df1)))
      beta_initial_list[[m]] = beta_initial_list[[m]] + beta_m_init/J 
    }
    t_init2 = proc.time()
    t_1_vec[j] = (t_init2- t_init1)[3]  
  }
  t_cost1= max(t_1_vec) + (t0b - t0a)[3]

  t1a = proc.time()
  Sigma_mj = z_mj = list()
  t_record = rep(NA, J)
  t1b = proc.time()
  t_cost1b = (t1b - t1a)[3]

  for(j in 1:J){
    t3 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j]
    for(m in 1:M){
      Xjm = cbind(1, as.matrix(X[xid_j, cset[[m]]]))
      Sigma_mj[[M*(j-1)+m]] = as.matrix( t(Xjm)%*%Xjm/nrow(Xjm) )
      z_mj[[M*(j-1)+m]] = t(Xjm)%*%Yj/nrow(Xjm)
    }
    t4 = proc.time()
    t_record[j] = (t4-t3)[3] 
  }
  t_cost2 = max(t_record)

  t5 = proc.time()
  Sigma_m = z_m = list()
  t6 = proc.time()
  t_cost3 = (t6-t5)[3]
 
  tim = 0
  for(m in 1:M){
     t_in1 = proc.time()
     Sigm = 0;  zm = 0
     t_in2 = proc.time()
     t_record2 = rep(NA, J)
     for(j in 1:J){
       t_in3 = proc.time()
       Sigm = Sigm + Sigma_mj[[m+M*(j-1)]]
       zm = zm + z_mj[[m+M*(j-1)]]
       t_in4 = proc.time()
       t_record2[j] = (t_in4 - t_in3)[3] 
     }
     t_in5 = proc.time()
     Sigma_m[[m]] = Sigm/J
     z_m[[m]] = zm/J 
     t_in6 = proc.time()
     tim = tim+ (t_in2-t_in1)[3] + max(t_record2) + (t_in6-t_in5)[3]
  }
  
  t_in7 = proc.time()
  beta_list =  list()
  for(m in 1:M){ 
    beta_list[[m]] = beta_initial_list[[m]]
    Sigma_1m_inverse = MASS::ginv(Sigma_mj[[m]])
    A = diag(rep(1, length(cset[[m]])+1)) - Sigma_1m_inverse%*%Sigma_m[[m]]
    B =  Sigma_1m_inverse%*%z_m[[m]]
    for(t in 1:T){         
      beta_list[[m]] = A%*%beta_list[[m]] +  B 
    }
  } 
  t_in8 = proc.time()
  t_cost = t_cost1 + t_cost1b + t_cost2 + t_cost3 + tim + (t_in8 - t_in7)[3]

  list(beta_list=beta_list, t_cost=t_cost) 
} 

dmap_sl2 = function(X, Y, X_new, Y_new, J, T=1, phi_type = NULL, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  n_new = nrow(X_new)
  res_beta = beta_csl_ma2(X=X, Y=Y, J=J, T=T, cm_type=cm_type, n_group=n_group)
  beta_cls_list =res_beta$beta_list
  t_cost1 = res_beta$t_cost   

  id_list = index_partition_fun(N=N, J=J)
  n = id_list[[2]]
  xid_list = id_list[[1]]
  if(is.null(phi_type)){ phi = 2 }else{ phi=log(n) }

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]
  k_vec = cset_list[[3]]+1

  ## compute on each machine  ##
  t1 = proc.time()
  w_est = matrix(nrow=M, ncol=J)
  mu_est = matrix(NA, N, M)
  mu_pred = matrix(NA, n_new*J, M)
  t2 = proc.time()
  t_record = rep(NA, J)
  for(j in 1:J){  
    t3 = proc.time()
    xid_j = xid_list[[j]]
    Yj = Y[xid_j] 
    for(m in 1:M){
      cmm = cset[[m]]
      Xjm = cbind(1, as.matrix(X[xid_j, cmm]))
      mu_est[xid_j, m] = as.vector(Xjm%*%beta_cls_list[[m]])
      mu_pred[(n_new*(j-1)+1):(n_new*j), m] = cbind(1, as.matrix(X_new[, cmm]))%*%beta_cls_list[[m]]
    }
    ## find weights ##
    muj_est = as.matrix(mu_est[xid_j, ])
    # Xfull = X[xid_j, ]
    df_full = data.frame(y=Yj, x=X[xid_j, ])
    sigmaj2_max = sigma(lm(y~., data=df_full))^2
    H0 = muj_est
    b0 = as.vector(t(Yj)%*%muj_est) - phi*sigmaj2_max*k_vec/2
    w_est[, j] = CLS(H=H0, b=b0)$solution
    t4 = proc.time()
    t_record[j] = (t4-t3)[3]
  } 

  t5 = proc.time()
  w_vec = rowMeans(w_est)
  mu_w_est = mu_est%*%w_vec
  mu_w_pred = rowMeans(matrix(as.vector(mu_pred%*%w_vec), ncol=J))
  mspe = mean( (Y_new - mu_w_pred)^2 )
  t6 = proc.time()
  t_cost = t_cost1 + (t2-t1)[3] + max(t_record) + (t6-t5)[3]

  w_ew = rep(1, M)/M      ## equal weights 
  mu_w_est_ew = mu_est%*%w_ew
  mu_w_pred_ew = rowMeans(matrix(as.vector(mu_pred%*%w_ew), ncol=J))
  mspe_ew = mean( (Y_new - mu_w_pred_ew)^2 )

  list(mspe = mspe, mu_w_est=mu_w_est, mspe_ew=mspe_ew, mu_w_est_ew=mu_w_est_ew, t_cost=t_cost)
} 



## Global Model Averaging Method ############
ma_global = function(X, Y, X_new, Y_new, phi_type = NULL, cm_type="nested", n_group=NULL){
  N = nrow(X) 
  p_max = ncol(X)
  n_new = nrow(X_new)
  if(is.null(phi_type)){ phi = 2 }else{ phi=log(N) }

  ## candidate model sets  ##
  cset_list = candset_fun(L=p_max, type=cm_type, n_group=n_group)
  cset = cset_list[[1]]
  M = cset_list[[2]]
  k_vec = cset_list[[3]]+1

  t1 = proc.time()
  mu_est = matrix(NA, N, M)
  mu_pred = matrix(NA, n_new, M)
  for(m in 1:M){
    cmm = cset[[m]]
    Xm = as.matrix(X[, cmm])
    dfm = data.frame(y=Y, x=Xm)
    lmfit_full = lm(y~., data=dfm)  
    # betam = coef( lm(y~, data=dfm) )
    mu_est[, m] = as.vector(fitted(lmfit_full))
    # mu_est[, m] = as.vector(Xm%*%betam)
    mu_pred[, m]  = as.vector(predict(lmfit_full, newdata=data.frame(x=X_new[, cmm])))
    # mu_pred[, m] = as.matrix(X_new[, cmm])%*%betam
  }
  ## find weights ##
  mu_est = as.matrix(mu_est)
  sigma2_max = sigma(lm(Y~X))^2
  H0 = mu_est
  b0 = as.vector(t(Y)%*%mu_est) - phi*sigma2_max*k_vec/2
  w_est = CLS(H=H0, b=b0)$solution

  mu_w_est = as.vector(mu_est%*%w_est)
  mu_w_pred = as.vector(mu_pred%*%w_est)
  mspe = mean( (Y_new - mu_w_pred)^2 )
  t2 = proc.time()
  t_cost = (t2-t1)[3]
  list(mspe = mspe, mu_w_est=mu_w_est, t_cost=t_cost)
}



##### Distributed prediction method based on a single model ##############

beta_csl_single = function(X, Y, J, T){
  N = nrow(X) 
  p_max = ncol(X)
  
  id_list = index_partition_fun(N=N, J=J)
  xid_list = id_list[[1]]

  ## compute initial estimator of beta on first machine  ##
  xid_1 = xid_list[[1]]
  Y1 = Y[xid_1] 
  X1 = as.matrix(X[xid_1, ])
  df1 = data.frame(y=Y1, x=X1)
  beta_initial = coef( lm(y~., data=df1))
  
  Sigma_j = z_j = list()
  for(j in 1:J){
    xid_j = xid_list[[j]]
    Yj = Y[xid_j]
    Xj = cbind(1, as.matrix(X[xid_j, ]))
    Sigma_j[[j]] = as.matrix( t(Xj)%*%Xj/nrow(Xj) )
    z_j[[j]] = t(Xj)%*%Yj/nrow(Xj)
  }

  Sig = 0;  z = 0
  for(j in 1:J){
    Sig = Sig + Sigma_j[[j]]
    z = z + z_j[[j]]
  }
  Sigma = Sig/J
  z = z/J 
  
  beta_est = beta_initial
  Sigma_1_inverse = MASS::ginv(Sigma_j[[1]])
  A = diag(rep(1, p_max+1)) - Sigma_1_inverse%*%Sigma
  B = Sigma_1_inverse%*%z
  for(t in 1:T){    
    beta_est = A%*%beta_est + B  
  } 
  return(beta_est) 
} 



dp_single = function(X, Y, X_new, Y_new, J, T=1){
  N = nrow(X) 
  p_max = ncol(X)
  n_new = nrow(X_new)

  beta_est = beta_csl_single(X=X, Y=Y, J=J, T=T)

  id_list = index_partition_fun(N=N, J=J)
  xid_list = id_list[[1]]

  mu_est = rep(NA, N)
  mu_pred = rep(NA, n_new*J)
  for(j in 1:J){  
    xid_j = xid_list[[j]]
    Xj = cbind(1, X[xid_j, ])
    mu_est[xid_j] = as.vector(Xj%*%beta_est)
    mu_pred[(n_new*(j-1)+1):(n_new*j)] = as.vector(cbind(1,X_new)%*%beta_est)
  }
  mspe = mean( (Y_new - mu_pred)^2 )
  list(mspe = mspe, mu_est=mu_est)
} 












