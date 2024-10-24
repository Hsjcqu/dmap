library(data.table)
library(ggplot2)
library(tidyverse)
library(scales)
library(glmnet)
library(doParallel)
library(foreach)

#######
####### data import and data processing
source("main functions.R")
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
h <- 32
J <- length(h)
#######
S <- 1000 ## repetition number
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
  names(data1) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8")
  
  ####### Lasso
  R_order2 <- func_Lasso(xs,ys)
  x2 <- xt[,R_order2]
  data2 <- as.data.frame(cbind(yt,x2))
  names(data2) <- c("y","x1","x2","x3","x4","x5","x6","x7","x8")
  
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
  
  sa1 <- dmap_sa(x.train1,y.train1,x.test1,y.test1,h)
  saw1 <- sa1$w_vec[8]
  betasa <- beta_csa_ma(x.train1,y.train1,h)$beta_list
  betasa1 <- betasa[[8]][-1]
  betasaw <- saw1*betasa1
  
  sl1 <- dmap_sl(x.train1,y.train1,x.test1,y.test1,h,T=3)
  slw1 <- sl1$w_vec[8]
  betasl <- beta_csl_ma(x.train1,y.train1,h,T=3)$beta_list
  betasl1 <- betasl[[8]][-1]
  betaslw <- slw1*betasl1
  
  no1 <- ma_global(x.train1,y.train1,x.test1,y.test1)
  now1 <- no1$w_vec[8]
  betano <- beta_ma(x.train1,y.train1)$beta_list
  betano1 <- betano[[8]][-1]
  betanow <- now1*betano1
  
  c(betasaw,betaslw,betanow)
}

stopCluster(cl)

result_matrix <- matrix(unlist(result), nrow = S, byrow = TRUE)

resultM <- colMeans(result_matrix)
std_devs <- apply(result_matrix, 2, sd)


result_df <- data.frame(cbind(resultM,std_devs))

