### Generate DGP #################


## Example 1: y is generated from an infinite-order model (p0=1000) p = [2n^{1/3}]###
DGPfun1 = function(n.train=2^12, n.test=50, p_max=NULL, p0=1000, rho=0.8, sigma=1){
  n = n.train + n.test
  if(is.null(p_max)){ p_max = 10 }
  # M = as.integer(2*n.train^{1/3})
  # R square: R^2 = 1/(1+sigma^2) >> sigma = sqrt((1-R^2)/R^2)
  xx = matrix(rnorm(n*p0), n, p0)
  corrmat = toeplitz(rho^(0:(p0-1)))
  cholmat = chol(corrmat)            #chol(A) = A'A
  xx0 = xx%*%cholmat 

  # system.time({
  #   covmat = outer(1:p0, 1:p0, function(u,v){ rho^(abs(u-v)) } )
  #   corrmat = toeplitz(rho^(0:(p0-2)))
  #   xx0 = MASS::mvrnorm(n, mu=rep(0, p0-1), Sigma=toeplitz(rho^(0:(p0-2))))
  #  })

  X = xx0
  linp = X%*%(1/(1:p0))
  mu = linp/sd(linp)
  Y = mu + rnorm(n, 0, sigma)
  id_train = 1:n.train
  id_test = (n.train+1):n
  list(x.train = X[id_train, 1:p_max], y.train=Y[id_train], 
       x.test = X[id_test, 1:p_max], y.test=Y[id_test], 
       mu_train = mu[id_train], mu_test = mu[id_test])
}

# dat1 = DGPfun1(n.train=2^12, n.test=50, p_max=12, p0=1000, rho=0.5, sigma=1)



## Example 2:  y is generated from a sparse model p=20 #####
DGPfun2 = function(n.train=2^12, n.test=50, p_max=15, p0=20, rho=0.8, sigma=1){
  n = n.train + n.test
  if(is.null(p_max)){ p_max = 10 }
  # R square: R^2 = 1/(1+sigma^2) >> sigma = sqrt((1-R^2)/R^2)
  xx = matrix(rnorm(n*(p0-1)), n, p0-1)
  corrmat = toeplitz(rho^(0:(p0-2)))
  cholmat = chol(corrmat)            #chol(A) = A'A
  xx0 = xx%*%cholmat

  X = xx0         #cbind(1, xx0)
  linp = 0.8*X[,1] +0.8*X[,6] + 0.2*X[,11]^2 + 0.2*sin(0.2*pi*X[, 16]) + 0.5*X[, 29]
  mu = linp/sd(linp)
  Y = mu + rnorm(n, 0, sigma)
  id_train = 1:n.train
  id_test = (n.train+1):n
  list(x.train = X[id_train, 1:p_max], y.train=Y[id_train], 
       x.test = X[id_test, 1:p_max], y.test=Y[id_test], 
       mu_train = mu[id_train], mu_test = mu[id_test])
}








































