library(data.table)
library(ggplot2)
library(tidyverse)
library(scales)
library(glmnet)
library(doParallel)
library(foreach)

#######
####### data import and data processing
source("mainfunctions.R")
df1 <- read.csv(file="housing.csv", header=T)
df2 <- df1[complete.cases(df1),] ##data cleaning
df2 <- df2[,-10]
names(df2) <- c("x1","x2","x3","x4","x5","x6","x7","x8","y")
data <- as.data.frame(round(scale(df2, center = T, scale = T), 4))
y <- data[,9]
x <- data[,-9] 
data <- data.frame(y,x)
N <- nrow(data)

###### the set of the number of machines
H <- c(2,4,8,16,32,64)
J <- length(H)
#######
S <- 1000 ## repetition number
msesa1 <- rep(0,J)
msesa2 <- rep(0,J)
msesa3 <- rep(0,J)

msesl1 <- rep(0,J)
msesl2 <- rep(0,J)
msesl3 <- rep(0,J)

msedp1 <- rep(0,J)

timesa1 <- rep(0,J)
timesa2 <- rep(0,J)
timesa3 <- rep(0,J)

timesl1 <- rep(0,J)
timesl2 <- rep(0,J)
timesl3 <- rep(0,J)


# Set up parallel backend
cores <- 12
cl <- makeCluster(cores) # Use one less than the total number of cores
registerDoParallel(cl)

result <- foreach(i = 1:S) %dopar% {
  #####sample
  s <- sample(c(1:nrow(data)),200)
  xs <- x[s,]
  ys <- y[s]
  
  xt <- x[-s,]
  yt <- y[-s]
  ####### R-Squared
  R_order <- func_lm(xs,ys)
  x1 <- xt[,R_order]
  data1 <- as.data.frame(cbind(yt, x1))
  
  ####### Lasso
  R_order2 <- func_Lasso(xs,ys)
  x2 <- xt[,R_order2]
  data2 <- as.data.frame(cbind(yt,x2))

  sub <- sample(c(1:nrow(data1)), 20032)
  train1 <- as.matrix(data1[sub,])
  y.train1 <- train1[,1]
  x.train1 <- train1[,-1]
  
  test1 <- as.matrix(data1[-sub,])
  y.test1 <- test1[,1]
  x.test1 <- test1[,-1]
  
  train2 <- as.matrix(data2[sub,])
  y.train2 <- train2[,1]
  x.train2 <- train2[,-1]
  
  test2 <- as.matrix(data2[-sub,])
  y.test2 <- test2[,1]
  x.test2 <- test2[,-1]
  j <- 0
  for (h in H) {
    j <- j+1
    sa1 <- dmap_sa(x.train1,y.train1,x.test1,y.test1,h)
    msesa1[j] <- sa1$mspe
    timesa1[j] <- sa1$t_cost[[1]]
    
    sa2 <- dmap_sa(x.train2,y.train2,x.test2,y.test2,h)
    msesa2[j]<- sa2$mspe
    timesa2[j]<- sa2$t_cost[[1]]
    
    sa3 <- dmap_sa(x.train1,y.train1,x.test1,y.test1,h,phi_type = NULL, cm_type="group", n_group=4)
    msesa3[j]<- sa3$mspe
    timesa3[j]<- sa3$t_cost[[1]]
    
    sl1 <- dmap_sl(x.train1,y.train1,x.test1,y.test1,h,T=3)
    msesl1[j]<- sl1$mspe
    timesl1[j]<- sl1$t_cost[[1]]
    
    sl2 <- dmap_sl(x.train2,y.train2,x.test2,y.test2,h,T=3)
    msesl2[j]<- sl2$mspe
    timesl2[j] <- sl2$t_cost[[1]]
    
    sl3 <- dmap_sl(x.train1,y.train1,x.test1,y.test1,h,T=3,phi_type = NULL,cm_type="group", n_group=4)
    msesl3[j]<- sl3$mspe
    timesl3[j] <- sl3$t_cost[[1]]
    
    dp1 <- dp_single(x.train1,y.train1,x.test1,y.test1,h,T=3)
    msedp1[j]<- dp1$mspe
  }
    no1 <- ma_global(x.train1,y.train1,x.test1,y.test1)
    mseno1<- no1$mspe
    timeno1<- no1$t_cost[[1]]
    
    no2 <- ma_global(x.train2,y.train2,x.test2,y.test2)
    mseno2<- no2$mspe
    timeno2 <- no2$t_cost[[1]]
    
    no3 <- ma_global(x.train1,y.train1,x.test1,y.test1,phi_type = NULL, cm_type="group", n_group=4)
    mseno3<- no3$mspe
    timeno3<- no3$t_cost[[1]]    
    
  c(msesa1,msesl1,mseno1,msesa2,msesl2,mseno2,msesa3,msesl3,mseno3,msedp1,timesa1,timesa2,timesa3,timesl1,timesl2,timesl3,timeno1,timeno2,timeno3)
}

stopCluster(cl)


resultM <- colMeans(matrix(unlist(result), nrow = S, byrow = TRUE))

result_df <- data.frame(t(resultM))
