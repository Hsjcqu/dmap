library(glmnet)
#install.packages("doParallel")
library(doParallel)
#install.packages("foreach")
library(foreach)


#######
####### data import and data processing
source("Hfunctions.R")


df1 = read.csv(file="housing.csv", header=T)
df2 = df1[complete.cases(df1), -10] ##data cleaning

dat = as.data.frame(scale(df2))
names(dat) = c("x1","x2","x3","x4","x5","x6","x7","x8","y")


set.seed(1)
id0 = 1:200
dat0 = dat[id0, ]
id = setdiff(1:dim(dat)[1], id0)
dat1 = dat[id, ]


p = dim(dat0)[2]-1

R2 = rep()
for(j in 1:p){
   lm_model = lm(y~., data=dat0[, c(j, 9)])
   sumlm = summary(lm_model)
   R2 = c(R2, sumlm$r.squared)
}


R_order = order(R2, decreasing = T)  #  8 1 2 4 7 5 6 3

dat1a = dat1[, c(R_order, 9)]


####### Lasso
library(glmnet)

x0 = dat0[, -9]
y0 = dat0[, 9]
f0 = glmnet(x0, y0, data=dat0, family="gaussian", nlambda=60, alpha=1)  
coef(f0)[-1, ]
bm = f0$beta
a = apply(as.matrix(bm[, ]), 2, function(u){which(u!=0)})
b = unique(a)
length(b)

R_order2 = c(8, 1, 4, 2, 7, 6, 3, 5) #first 200 obs

dat1b = dat1[, c(R_order2, 9)]


###### the set of the number of machines
J <- c(2, 4, 8, 16, 32, 64)
n_J <- length(J)

####### the results of the first candidate model
####### the MSPEs of the first candidate model

n_train = 19968

S <- 1000 ## repetition number

# Set up parallel backend
cores <- 48
cl <- makeCluster(cores) # Use one less than the total number of cores
registerDoParallel(cl)

result <- foreach(i = 1:S) %dopar% {
  set.seed(i)
  sub <- sample(1:nrow(dat1a), n_train)
  train1 <- as.matrix(dat1a[sub,])
  y.train1 <- train1[, 9]
  x.train1 <- train1[,-9]
  test1 <- as.matrix(dat1a[-sub,])
  y.test1 <- test1[, 9]
  x.test1 <- test1[,-9]

  train2 <- as.matrix(dat1b[sub,])
  y.train2 <- train2[, 9]
  x.train2 <- train2[,-9]
  test2 <- as.matrix(dat1b[-sub,])
  y.test2 <- test2[, 9]
  x.test2 <- test2[,-9]

  mspe_dp <- rep(0,n_J)
  mspe_sa1 <- mspe_sa2 <- mspe_sa3 <- rep(0,n_J)
  mspe_sl1 <- mspe_sl2 <- mspe_sl3 <- rep(0,n_J)
  time_sa1 <- time_sa2 <- time_sa3 <- rep(0, n_J)
  time_sl1 <- time_sl2 <- time_sl3 <- rep(0, n_J)

  j <- 1
  repeat{
    sa1 <- dmap_sa(x.train1, y.train1, x.test1, y.test1, J=J[j])
    mspe_sa1[j] <- sa1$mspe
    time_sa1[j] <- sa1$t_cost[[1]]
    
    sa2 <- dmap_sa(x.train2, y.train2, x.test2, y.test2, J=J[j])
    mspe_sa2[j]<- sa2$mspe
    time_sa2[j]<- sa2$t_cost[[1]]
    
    sa3 <- dmap_sa(x.train1, y.train1, x.test1, y.test1, J=J[j], phi_type = NULL, cm_type="group", n_group=4)
    mspe_sa3[j]<- sa3$mspe
    time_sa3[j]<- sa3$t_cost[[1]]
    
    sl1 <- dmap_sl(x.train1, y.train1, x.test1, y.test1, J=J[j], T=3)
    mspe_sl1[j]<- sl1$mspe
    time_sl1[j]<- sl1$t_cost[[1]]
    
    sl2 <- dmap_sl(x.train2, y.train2, x.test2, y.test2, J=J[j], T=3)
    mspe_sl2[j]<- sl2$mspe
    time_sl2[j] <- sl2$t_cost[[1]]
    
    sl3 <- dmap_sl(x.train1, y.train1, x.test1, y.test1, J=J[j], T=3, phi_type = NULL,cm_type="group", n_group=4)
    mspe_sl3[j]<- sl3$mspe
    time_sl3[j] <- sl3$t_cost[[1]]
    
    dp1 <- dp_single(x.train1, y.train1, x.test1, y.test1, J=J[j], T=3)
    mspe_dp[j]<- dp1$mspe
   
    cat("J = ", J[j], "\n")
    if(j >=n_J){ break}else{ j=j+1 }
  }

  no1 <- ma_global(x.train1, y.train1, x.test1, y.test1)
  mspe_no1<- no1$mspe
  time_no1<- no1$t_cost[[1]]
    
  no2 <- ma_global(x.train2, y.train2, x.test2, y.test2)
  mspe_no2<- no2$mspe
  time_no2 <- no2$t_cost[[1]]
    
  no3 <- ma_global(x.train1, y.train1, x.test1, y.test1, phi_type = NULL, cm_type="group", n_group=4)
  mspe_no3<- no3$mspe
  time_no3<- no3$t_cost[[1]]    
    
  c(mspe_sa1, mspe_sa2, mspe_sa3, 
    mspe_sl1, mspe_sl2, mspe_sl3, 
    mspe_no1, mspe_no2, mspe_no3, mspe_dp, 
    time_sa1, time_sa2, time_sa3, 
    time_sl1, time_sl2, time_sl3, 
    time_no1, time_no2, time_no3)
}

stopCluster(cl)

# resultM <- colMeans(matrix(unlist(result), nrow = S, byrow = TRUE))

re = matrix(NA, S, 13*n_J+6)
for(s in 1:S){
  re[s, ] = result[[s]]
}

re_mean = apply(re, 2, mean)
rem = cbind(matrix(re_mean[1:(6*n_J)], nrow=n_J), 
       matrix(rep(re_mean[(6*n_J+1):(6*n_J+3)], each=n_J), nrow=n_J),
       matrix(re_mean[(6*n_J+3+1):(7*n_J+3)], nrow=n_J),
       matrix(re_mean[(7*n_J+4):(13*n_J+3)], nrow=n_J),
       matrix(rep(re_mean[(13*n_J+4):(13*n_J+6)], each=n_J), nrow=n_J))  
mspe_J = rem[, 1:10]
time_J = rem[, 11:19]
colnames(mspe_J) = c("DMAP-SA(CM1)", "DMAP-SA(CM2)", "DMAP-SA(CM3)", 
                     "DMAP-SL(CM1)", "DMAP-SL(CM2)", "DMAP-SL(CM3)",
                     "gMAP(CM1)", "gMAP(CM2)", "gMAP(CM3)", "DP")
rownames(mspe_J) = paste("J=", J, sep="")
colnames(time_J) = c("DMAP-SA(CM1)", "DMAP-SA(CM2)", "DMAP-SA(CM3)", 
                     "DMAP-SL(CM1)", "DMAP-SL(CM2)", "DMAP-SL(CM3)",
                     "gMAP(CM1)", "gMAP(CM2)", "gMAP(CM3)")
rownames(time_J) = paste("J=", J, sep="")

mspe1 = t(mspe_J[, c(1, 4, 7, 2, 5, 8, 3, 6, 9, 10)])
time1 = time_J[, c(1, 4, 7, 2, 5, 8, 3, 6, 9)]

write.csv(mspe1, file = "mspe.csv", row.names = TRUE)

write.csv(time1, file = "time.csv", row.names = TRUE)
