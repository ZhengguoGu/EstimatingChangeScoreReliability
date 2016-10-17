library(Lambda4)
library(psychometric)

# data simulation

################################################

set.seed(110)

#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

parallel_items <- 2 # 1: parallel, otherwise nonparallel
dimension <- c(1, 2, 4) # number of dimensions in theta
num_persons <- 1000

num_items <- 10

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
  } else if (dimension == 4){
    itempar <- matrix(NA,num_items,6)
    itempar[,1] <- runif(num_items,1,3)  #the parameters are set according to Wilco
    itempar[,2] <- itempar[,1] # for this moment, let discimination for each item to be the same across dimenions. 
    itempar[,3] <- itempar[,1]
    itempar[,4] <- itempar[,1]
    itempar[,5] <- runif(num_items,-2,0)
    itempar[,6] <- runif(num_items,0.5,2)
  }
}


#####################################################
#
#  generate random samples from population (theta)
#
#####################################################

r_pupulation <- list()

p <- 1
while(p<=100) {
  
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
    r_simresults <- matrix(NA, n_sim, 8)  # 8 methods.
  
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
    
      #####################################################################################################################################################
      # Methods for calculating reliability
      #
      # method 0.1: ture-change reliability - var(true change)/var(observed change)
      # method 0.2: ture-change reliability - cor(true change, observed change)
      # 
      # method 1.1: estimated reliability - alpha (i.e. pre and post reliability estimated by alpha )
      # method 1.2: estimated reliability - lambda2 (i.e. pre and post reliability estimated by lambda2)
      # method 1.3: estimated reliability - lambda4 (i.e. pre and post reliability estimated by lambda4)
      #
      # method 2.1: estimated reliability (item-level) - alpha (i.e. change scores at item level are used to estimate reliability by means of alpha) 
      # method 2.2: estimated reliability (item-level) - lambda2 (i.e. change scores at item level are used to estimate reliability by means of lambda2) 
      # method 2.3: estimated reliability (item-level) - lambda4 (i.e. change scores at item level are used to estimate reliability by means of lambda4) 
      ###################################################################################################################################################
    
      ######## method 0.1: ture-change reliability - var(true change)/var(observed change) ###########
    
      truechange_sumscores <- sum_true_post - sum_true_pre
      change_sumscores <- sum_post - sum_pre
    
      r_simresults[i, 1] <- var(truechange_sumscores)/var(change_sumscores)
    
      ######## method 0.2: ture-change reliability - cor(true change, observed change) ###########
    
      r_simresults[i, 2] <- (cor(truechange_sumscores, change_sumscores))^2
    
      ######## method 1.1: estimated reliability - alpha (i.e. pre and post reliability estimated by alpha )  #########################
    
      r_pre <- psychometric::alpha(response_pre)  # ! cronback alpha is used here. 
      r_post <- psychometric::alpha(response_post)
    
      r_simresults[i, 3] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
    
      ######## method 1.2: estimated reliability - lambda2 (i.e. pre and post reliability estimated by lambda2) #######
    
      r_pre <- Lambda4::lambda2(response_pre)  # ! lambda2 is used here. 
      r_post <- Lambda4::lambda2(response_post)
    
      r_simresults[i, 4] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
    
      ######## method 1.3: estimated reliability - lambda4 (i.e. pre and post reliability estimated by lambda4)
    
      r_simresults[i, 5] <- 0
    
      ######## method 2.1: estimated reliability (item-level) - alpha (i.e. change scores at item level are used to estimate reliability by means of alpha)   #########################
    
      response_change <- response_post - response_pre
      r_simresults[i, 6] <- psychometric::alpha(response_change)
    
      ######## method 2.2: estimated reliability (item-level) - lambda2 (i.e. change scores at item level are used to estimate reliability by means of lambda2) ################
    
      r_simresults[i, 7] <- Lambda4::lambda2(response_change)
    
      ########method 2.3: estimated reliability (item-level) - lambda4 (i.e. change scores at item level are used to estimate reliability by means of lambda4) 
    
      r_simresults[i, 8] <- 0
    
    }
  
    r_forplots[[L]] <- r_simresults
  
  }
  r_pupulation[[p]] <- r_forplots
  P <- P+1
}

###################################################################################
#
#  A few plots
#
###################################################################################

#results <- r_forplots[[1]]
#colnames(results) <- c('true_reliability_var', 'true_reliability_cor', 'method1', 'method2', 'method3')

#plot(results[,1], results[, 2], xlab = "true reliability - proportion variance",
#                     ylab = "true reliability - correlation", 
#                     main = 'Dimension = 2', xlim=c(0, 1), ylim=c(0, 1))
#lines(c(0,1),c(0,1),col="red")


