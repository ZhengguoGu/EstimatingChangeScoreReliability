# data simulation

################################################

set.seed(110)

#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

parallel_items <- 2 # 1: parallel, otherwise nonparallel
dimension <- c(1, 2, 4) # number of dimensions in theta
num_persons <- c(4000, 2000, 1000)

num_items <- 10

if (parallel_items == 1) {
  itempar <- matrix(NA,num_items,3)
  itempar[,1] <- runif(1,1,3)   # discrimination
  itempar[,2] <- runif(1,-2,0)  # difficulty  
  itempar[,3] <- runif(1,0.5,2) # difficulty
} else {
  itempar <- matrix(NA,num_items,3)
  itempar[,1] <- runif(num_items,1,3)  #the parameters are set according to Wilco
  itempar[,2] <- runif(num_items,-2,0)
  itempar[,3] <- runif(num_items,0.5,2)
}

r_forplots <- list()

for (L in 1:length(dimension)){
  
  #-------------------------------------------------
  # simulate multidimensional theta's
  #-------------------------------------------------

  if (dimension[L] == 1){
  
    
    n_sub <- num_persons[1]
  
    sd_pre <- 1
    mean_change <- 1
    sd_change <- 0.3
  
    theta <- Unichange_sim(n_sub, sd_pre, mean_change, sd_change)
    
    theta_pre <- theta[[1]]
    theta_post <- theta[[2]]
  
  } else{
  
    cov_pretest <- 0.8 # need to justify why 0.8  
    mean_change <- 1
    sd_change <- 0.3 # need to justify why 0.3
    EMP <- FALSE
    
    n_sub <- num_persons[L]
    theta <- Mulchange_sim(n_sub, dimension[L], cov_pretest, mean_change, sd_change, EMP)

  }

  theta_pre <- theta[[1]]
  theta_post <- theta[[2]]

  #-------------------------------------------------------------------

  n_sim <- 1000 # simulate 1000 datasets
  r_simresuls <- matrix(NA, n_sim, 3)  # 2 methods. 
  for (i in 1:n_sim){
  #-------------------------------------------------
  # simulate graded response data
  #-------------------------------------------------

  response_pre <- GRM_sim(as.vector(theta_pre), itempar)
  response_post <- GRM_sim(as.vector(theta_post), itempar)
  # estimate the response data by means of function grm from package itm

  #library(ltm)
  #estimates_pre <- grm(response_pre)
  #print(estimates_pre)
  #estimates_post <- grm(response_post)
  #print(estimates_post)

  #--------------------------------------------------
  # Calculate inter-item covariance matrix
  # for observed scores
  #--------------------------------------------------

  #cov_obsmatrix <- cov(response_pre, response_post)

  #-------------------------------------------------
  # calculate change-score reliability
  #-------------------------------------------------

  # sum scores
  sum_pre <- rowSums(response_pre)
  sum_post <- rowSums(response_post)

  ######## method 1  #########################

  r_pre <- psychometric::alpha(response_pre)  # ! cronback alpha is used here. 
  r_post <- psychometric::alpha(response_post)

  r_simresuls[i, 1] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))

  ######## method 2  #########################

  response_change <- response_post - response_pre
  r_simresuls[i, 2] <- psychometric::alpha(response_change)

  ######## method 3: lambda 2 ################

  r_simresuls[i, 3] <- Lambda4::lambda2(response_change)

  }

  r_forplots[[L]] <- r_simresuls

}

#plot1v2[[L]] <- plot(r_simresuls[,1], r_simresuls[, 2], xlab = "Estimated change-score reliability based on Method 1", ylab = "Estimated change-score reliability based on Method 2")
#lines(c(0,1),c(0,1),col="red")

#plot1v3[[L]] <- plot(r_simresuls[,1], r_simresuls[, 3], xlab = "Estimated change-score reliability based on Method 1", ylab = "Estimated change-score reliability based on Method 3")
#lines(c(0,1),c(0,1),col="red")