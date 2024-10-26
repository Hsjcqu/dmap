library(glmnet)
#install.packages("doParallel")
library(doParallel)
#install.packages("foreach")
library(foreach)


#######
####### data import and data processing
source("main functions.R")

df1 = read.csv(file="housing.csv", header=T)
df2 = df1[complete.cases(df1), -10] ##data cleaning
str(df2)

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
J = c(8, 32)
n_J = length(J)


####### the results of the first candidate model
####### the MSPEs of the first candidate model

n_train = 19968

S <- 1000 ## repetition number


# Set up parallel backend
cores <- 48
cl <- makeCluster(cores) # Use one less than the total number of cores
registerDoParallel(cl)

result <- foreach(i = 1:S) %dopar% {
  sub <- sample(1:nrow(dat1a), n_train)
  train1 <- as.matrix(dat1a[sub,])
  y.train1 <- train1[, 9]
  x.train1 <- train1[,-9]

  train2 <- as.matrix(dat1b[sub,])
  y.train2 <- train2[, 9]
  x.train2 <- train2[,-9]

  beta_w_sa1 = beta_w_sa2 = beta_w_sa3 = c()
  beta_w_sl1 = beta_w_sl2 = beta_w_sl3 = c()
  
  j <- 1
  repeat{
    sa1 = beta_dmap_sa(X=x.train1, Y= y.train1, J=J[j], phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
    beta_w_sa1 = c(beta_w_sa1, sa1)
    sa2 = beta_dmap_sa(X=x.train2, Y= y.train2, J=J[j], phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
    beta_w_sa2 = c(beta_w_sa2, sa2)
    sa3 = beta_dmap_sa(X=x.train1, Y= y.train1, J=J[j], phi_type = NULL, cm_type="group", n_group=4)$beta_w[,1]
    beta_w_sa3 = c(beta_w_sa3, sa3)
    
    sl1 = beta_dmap_sl(X=x.train1, Y= y.train1, J=J[j], T=3, phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
    beta_w_sl1 = c(beta_w_sl1, sl1)
    sl2 = beta_dmap_sl(X=x.train2, Y= y.train2, J=J[j], T=3, phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
    beta_w_sl2 = c(beta_w_sl2, sl2)
    sl3 = beta_dmap_sl(X=x.train1, Y= y.train1, J=J[j], T=3, phi_type = NULL, cm_type="group", n_group=4)$beta_w[,1]
    beta_w_sl3 = c(beta_w_sl3, sl3)

    #cat("J = ", J[j], "\n")
    if(j >=n_J){ break}else{ j=j+1 }
  }

  beta_w_gmap1 = beta_ma_global(X=x.train1, Y= y.train1, phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
  beta_w_gmap2 = beta_ma_global(X=x.train2, Y= y.train2, phi_type = NULL, cm_type="nested", n_group=NULL)$beta_w[,1]
  beta_w_gmap3 = beta_ma_global(X=x.train1, Y= y.train1, phi_type = NULL, cm_type="group", n_group=4)$beta_w[,1]
    
  c(beta_w_sa1, beta_w_sa2, beta_w_sa3, 
    beta_w_sl1, beta_w_sl2, beta_w_sl3, 
    beta_w_gmap1, beta_w_gmap2, beta_w_gmap3)
}

stopCluster(cl)


re = matrix(NA, S, p*6*n_J + p*3)
for(s in 1:S){
  re[s, ] = result[[s]]
}

re_mean = apply(re, 2, mean)
re1 = matrix(re_mean, nrow=p)  

re_sd = apply(re, 2, sd)
re2 = matrix(re_sd, nrow=p)  


colnames(re1) = c(paste("DMAP-SA(CM1)J=", J, sep=""), 
                  paste("DMAP-SA(CM2)J=", J, sep=""), 
                  paste("DMAP-SA(CM3)J=", J, sep=""), 
                  paste("DMAP-SL(CM1)J=", J, sep=""), 
                  paste("DMAP-SL(CM2)J=", J, sep=""), 
                  paste("DMAP-SL(CM3)J=", J, sep=""), 
                  "gMAP(CM1)", "gMAP(CM2)", "gMAP(CM3)")
colnames(re2) = c(paste("DMAP-SA(CM1)J=", J, sep=""), 
                  paste("DMAP-SA(CM2)J=", J, sep=""), 
                  paste("DMAP-SA(CM3)J=", J, sep=""), 
                  paste("DMAP-SL(CM1)J=", J, sep=""), 
                  paste("DMAP-SL(CM2)J=", J, sep=""), 
                  paste("DMAP-SL(CM3)J=", J, sep=""), 
                  "gMAP(CM1)", "gMAP(CM2)", "gMAP(CM3)")


re3 = cbind(re1, re2)
beta_cm1 = re3[, c(1, 1+15, 7, 7+15, 2, 2+15, 8, 8+15, 13, 13+15)]
beta_cm2 = re3[, c(3, 3+15, 9, 9+15, 4, 4+15, 10, 10+15, 14, 14+15)]
beta_cm3 = re3[, c(5, 5+15, 11, 11+15, 6, 6+15, 12, 12+15, 15, 15+15)]

write.csv(beta_cm1, file="beta cm1.csv")
write.csv(beta_cm2, file="beta cm2.csv")
write.csv(beta_cm3, file="beta cm3.csv")

