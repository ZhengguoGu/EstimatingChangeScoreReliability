#######################################################
### Testing file
#######################################################

#### This is to test a few functions in populationRel.cpp. 
Rcpp::sourceCpp("populationRel.cpp")

#### Test  GRMc_1theta
set.seed(1)
theta <- rnorm(100000, 0, 1)
num_items <- 9
itempar <- matrix(NA,num_items,5)
itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
avg_beta <- runif(num_items, -1.5, 1.5)
itempar[,2] <- avg_beta - .75
itempar[,3] <- avg_beta - .25
itempar[,4] <- avg_beta + .25
itempar[,5] <- avg_beta + .75

result <- matrix(0, 100000, 10)
for (n in 1:100000){
  
  result[n,] <- GRMc_1theta(theta[n], itempar[,1], itempar[,2:5])
  
}

library(ltm)

parameteres<-ltm::grm(data.matrix(result))  #pretest scores

est_par <- matrix(unlist(parameteres$coefficients), ncol = 5, byrow = TRUE)
est_par <- cbind(est_par[, 5], est_par[, -5])
plot(itempar[, 1], est_par[, 1])  #slope
plot(itempar[, 2], est_par[, 2])  #item-difficulty 1
plot(itempar[, 3], est_par[, 3])  #item-difficulty 2
plot(itempar[, 4], est_par[, 4])  #item-difficulty 3
plot(itempar[, 5], est_par[, 5])  #item-difficulty 4 


###########################################################

#### Test function GRMc_1thetaMD
EMP <- FALSE
num_items <- 21
itempar <- matrix(NA,num_items,5)
itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
avg_beta <- runif(num_items, -1.5, 1.5)
itempar[,2] <- avg_beta - .75
itempar[,3] <- avg_beta - .25
itempar[,4] <- avg_beta + .25
itempar[,5] <- avg_beta + .75
N_pop <- 1000000
dimension <- 3
cov_pretest <- .5
sd_change <- .14
theta_pop <- Mulchange_sim(N_pop, dimension, cov_pretest, sd_change, EMP) #from ZhengguoFunction.R
theta_pre_pop <- theta_pop[[1]]

id <- vector()
for(d in 1: dimension){
  id <- cbind(id, rep(d, num_items/dimension))
}
id <- as.vector(id)
pretest_obs <- matrix(NA, N_pop, num_items)

for(i in 1:N_pop){
  abil_pre <- theta_pre_pop[i, ][id]
  pretest_obs[i, ] <- GRMc_1thetaMD(abil_pre, itempar[,1], itempar[, 2:5])
}


library(ltm)

parameteres<-ltm::grm(pretest_obs[, 1:7]) 
est_par <- matrix(unlist(parameteres$coefficients), ncol = 5, byrow = TRUE)
est_par <- cbind(est_par[, 5], est_par[, -5])
plot(itempar[1:7, 1], est_par[, 1])  #slope
plot(itempar[1:7, 2], est_par[, 2])  #item-difficulty 1
plot(itempar[1:7, 3], est_par[, 3])  #item-difficulty 2
plot(itempar[1:7, 4], est_par[, 4])  #item-difficulty 3
plot(itempar[1:7, 5], est_par[, 5])  #item-difficulty 4 

parameteres2<-ltm::grm(pretest_obs[, 8:14])
est_par2 <- matrix(unlist(parameteres2$coefficients), ncol = 5, byrow = TRUE)
est_par2 <- cbind(est_par2[, 5], est_par2[, -5])
plot(itempar[8:14, 1], est_par2[, 1])  #slope
plot(itempar[8:14, 2], est_par2[, 2])  #item-difficulty 1
plot(itempar[8:14, 3], est_par2[, 3])  #item-difficulty 2
plot(itempar[8:14, 4], est_par2[, 4])  #item-difficulty 3
plot(itempar[8:14, 5], est_par2[, 5])  #item-difficulty 4 

parameteres3<-ltm::grm(pretest_obs[, 15:21])
est_par3 <- matrix(unlist(parameteres3$coefficients), ncol = 5, byrow = TRUE)
est_par3 <- cbind(est_par3[, 5], est_par3[, -5])
plot(itempar[15:21, 1], est_par3[, 1])  #slope
plot(itempar[15:21, 2], est_par3[, 2])  #item-difficulty 1
plot(itempar[15:21, 3], est_par3[, 3])  #item-difficulty 2
plot(itempar[15:21, 4], est_par3[, 4])  #item-difficulty 3
plot(itempar[15:21, 5], est_par3[, 5])  #item-difficulty 4 



#### Test function Carryover
pretest_obs <- 1:20
posttest_obs <- 20:1
Eff <- 'W'

Carryover(pretest_obs, posttest_obs, Eff)  #C++ 
pretest_obs
posttest_obs

carry_over(pretest_obs, posttest_obs, proc_N = 1)  #R function in ZhengguoFunction.R



