library(Lambda4)

# data simulation

################################################

set.seed(110)

#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

parallel_items <- 2 # 1: parallel, otherwise nonparallel
dimension <- c(1,2) # number of dimensions in theta
num_persons <- 1000

num_items <- 50

if (parallel_items == 1) {
  itempar <- matrix(NA,num_items,3)
  itempar[,1] <- runif(1,1,3)   # discrimination
  itempar[,2] <- runif(1,-2,0)  # difficulty  
  itempar[,3] <- runif(1,0.5,2) # difficulty
} else {
  if (dimension == 1){
    
    itempar <- matrix(NA,num_items,3)
    itempar[,1] <- runif(num_items,1,3)  #the parameters are set according to Wilco
    itempar[,2] <- runif(num_items,-2,0)
    itempar[,3] <- runif(num_items,0.5,2)
  } else if (dimension == 2){
    itempar <- matrix(NA,num_items,4)
    itempar[,1] <- runif(num_items,1,3)  #the parameters are set according to Wilco
    itempar[,2] <- itempar[,1] # for this moment, let discimination for each item to be the same across dimenions. 
    itempar[,3] <- runif(num_items,-2,0)
    itempar[,4] <- runif(num_items,0.5,2)
  }
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
    
    n_sub <- num_persons
    theta <- Mulchange_sim(n_sub, dimension[L], cov_pretest, mean_change, sd_change, EMP)

  }

  theta_pre <- theta[[1]]
  theta_post <- theta[[2]]

  #-------------------------------------------------------------------

  n_sim <- 500 # simulate 1000 datasets
  r_simresults <- matrix(NA, n_sim, 4)  # 3 methods. 
  for (i in 1:n_sim){
  #-------------------------------------------------
  # simulate graded response data
  #-------------------------------------------------

    responses <- GRM_sim(theta_pre, itempar)
    response_pre <- responses[[1]]
    true_pre <- responses[[2]]
  
    responses <- GRM_sim(theta_post, itempar)
    response_post <- responses[[1]]
    true_post <- responses[[2]]
    #Ides = rep(1,num_items)
    #X1 = FsimMDGRM(as.matrix(theta_pre),itempar,Ides)  # pretest
    #X2 = FsimMDGRM(as.matrix(theta_post),itempar,Ides)  # posttest
    #true_pre = X1[[1]]
    #true_post = X2[[1]]
    #response_pre = X1[[2]]
    #response_post= X2[[2]]
    
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
  sum_true_pre <- rowSums(true_pre)
  sum_post <- rowSums(response_post)
  sum_true_post <- rowSums(true_post)
  
  ######## method 0: the true reliability ###########
  
  truechange_sumscores <- sum_true_post - sum_true_pre
  change_sumscores <- sum_post - sum_pre
  
  r_simresults[i, 1] <- var(truechange_sumscores)/var(change_sumscores)
  
  ######## method 1  #########################

  r_pre <- psychometric::alpha(response_pre)  # ! cronback alpha is used here. 
  r_post <- psychometric::alpha(response_post)

  r_simresults[i, 2] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))

  ######## method 2  #########################

  response_change <- response_post - response_pre
  r_simresults[i, 3] <- psychometric::alpha(response_change)

  ######## method 3: lambda 2 ################

  r_simresults[i, 4] <- Lambda4::lambda2(response_change)


  }

  r_forplots[[L]] <- r_simresults

}

#------------------------
results <- r_forplots[[1]]
colnames(results) <- c('true_reliability', 'method1', 'method2', 'method3')

plot(results[,2], results[, 3], xlab = "Estimated change-score reliability based on Method 1",
                     ylab = "Estimated change-score reliability based on Method 2", 
                     main = 'Dimension = 1')
lines(c(0,1),c(0,1),col="red")

plot(results[,2], results[, 4], xlab = "Estimated change-score reliability based on Method 1", 
     ylab = "Estimated change-score reliability based on Method 3",
     main = 'Dimension = 1')
lines(c(0,1),c(0,1),col="red")

plot(results[,2], results[, 1], xlab = "Estimated change-score reliability based on Method 1",
     ylab = "True reliability", 
     main = 'Dimension = 1')
lines(c(0,1),c(0,1),col="red")

#------------------------
results <- r_forplots[[2]]

plot(results[,2], results[, 3], xlab = "Estimated change-score reliability based on Method 1",
     ylab = "Estimated change-score reliability based on Method 2", 
     main = 'Dimension = 2')
lines(c(0,1),c(0,1),col="red")

plot(results[,2], results[, 4], xlab = "Estimated change-score reliability based on Method 1", 
     ylab = "Estimated change-score reliability based on Method 3",
     main = 'Dimension = 2')
lines(c(0,1),c(0,1),col="red")

plot(results[,2], results[, 1], xlab = "Estimated change-score reliability based on Method 1",
     ylab = "True reliability", 
     main = 'Dimension = 2')
lines(c(0,1),c(0,1),col="red")

#-----------------------------

#results <- r_forplots[[3]]

#plot(results[,2], results[, 3], xlab = "Estimated change-score reliability based on Method 1",
#     ylab = "Estimated change-score reliability based on Method 2", 
#     main = 'Dimension = 4', ylim = c(0.2, 0.4), xlim = c(0.2, 0.4))
#lines(c(0,1),c(0,1),col="red")

#plot(results[,2], results[, 4], xlab = "Estimated change-score reliability based on Method 1", 
#     ylab = "Estimated change-score reliability based on Method 3",
#     main = 'Dimension = 4', ylim = c(0.2, 0.4), xlim = c(0.2, 0.4))
#lines(c(0,1),c(0,1),col="red")

#plot(results[,2], results[, 1], xlab = "Estimated change-score reliability based on Method 1",
#     ylab = "True reliability", 
#     main = 'Dimension = 4', ylim = c(0, 0.4), xlim = c(0, 0.4))
#lines(c(0,1),c(0,1),col="red")
#-----------------------------

#results <- r_forplots[[4]]

#plot(results[,2], results[, 3], xlab = "Estimated change-score reliability based on Method 1",
#     ylab = "Estimated change-score reliability based on Method 2", 
#     main = 'Dimension = 8', ylim = c(0.2, 0.4), xlim = c(0.2, 0.4))
#lines(c(0,1),c(0,1),col="red")

#plot(results[,2], results[, 4], xlab = "Estimated change-score reliability based on Method 1", 
#     ylab = "Estimated change-score reliability based on Method 3",
#     main = 'Dimension = 8', ylim = c(0.2, 0.4), xlim = c(0.2, 0.4))
#lines(c(0,1),c(0,1),col="red")

#plot(results[,2], results[, 1], xlab = "Estimated change-score reliability based on Method 1",
#     ylab = "True reliability", 
#     main = 'Dimension = 8', ylim = c(0, 0.4), xlim = c(0, 0.4))
#lines(c(0,1),c(0,1),col="red")
