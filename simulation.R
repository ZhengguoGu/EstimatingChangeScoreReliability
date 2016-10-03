# data simulation

################################################

set.seed(110)


#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

parallel_items <- 1 # 1: parallel, otherwise nonparallel

num_items <- 10

if (parallel_items == 1) {
  itempar <- matrix(NA,num_items,3)
  itempar[,1] <- runif(1,1,3)   # discrimination
  itempar[,2] <- runif(1,-2,0)  # difficulty  
  itempar[,3] <- runif(1,0.5,2) # difficulty
} else {
  itempar <- matrix(NA,num_items,3)
  itempar[,1] <- runif(num_items,1,3)  #the parameters are set according to Wilco
  itempar[,2] <- runif(num_items0,-2,0)
  itempar[,3] <- runif(num_items,0.5,2)
}

#-------------------------------------------------
# simulate multidimensional theta's
#-------------------------------------------------

dimension <- 3
cov_pretest <- 0.8 # need to justify why 0.8  
mean_change <- 1
sd_change <- 0.3 # need to justify why 0.3
EMP <- FALSE
num_persons <- 1000
theta <- Mulchange_sim(num_persons, dimension, cov_pretest, mean_change, sd_change, EMP)

theta_pre <- theta[[1]]
theta_post <- theta[[2]]

#-------------------------------------------------
# simulate graded response data
#-------------------------------------------------

response_pre <- GRM_sim(as.vector(theta_pre), itempar)
response_post <- GRM_sim(as.vector(theta_post), itempar)
# estimate the response data by means of function grm from package itm

library(ltm)
estimates_pre <- grm(response_pre)
print(estimates_pre)
estimates_post <- grm(response_post)
print(estimates_post)
