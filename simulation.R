library(Lambda4)
library(psychometric)

# data simulation

################################################

set.seed(110)

#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

num_items <- 9 # test length: 9, 18, 27, 36 items

parallel_items <- 1 # 1: parallel, otherwise nonparallel

dimension <- 3 # number of dimensions in theta

cov_pretest <- 0.8  # 0.3, 0.5, 0.8
sd_change <- 0.14 # 0.14, 0.5

existence_carryover <- 1 # 1 = yes, 0 = no
strong_weak <- 1 # 1 = strong, -1 = weak carry-over effects

num_persons <- 1000 # number of subjects
n_sim <- 1000 # simulate 1000 datasets

if (parallel_items == 1) {
  
  itempar <- matrix(NA,num_items,5)
  itempar[,1] <- runif(1,1.5,2.5)   # discrimination
  avg_beta <- runif(1, 0, 1.25)
  itempar[,2] <- avg_beta - 1
  itempar[,3] <- avg_beta - .5
  itempar[,4] <- avg_beta + .5
  itempar[,5] <- avg_beta + 1
  
} else {
  
  itempar <- matrix(NA,num_items,5)
  itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
  avg_beta <- runif(num_items, 0, 1.25)
  itempar[,2] <- avg_beta - 1
  itempar[,3] <- avg_beta - .5
  itempar[,4] <- avg_beta + .5
  itempar[,5] <- avg_beta + 1
  
}

id <- vector()
for(d in 1: dimension){
  id <- cbind(id, rep(d, num_items/dimension))
}
id <- as.vector(id)
#####################################################
#
#  generate random samples from population (theta)
#
#####################################################

sample_results <- list()
maxp <- 10
num_methods <- 8 #see below, method 0.1 to 2.3
r_avg <- matrix(NA, maxp, num_methods)
p <- 1

while(p<=maxp) {

  #-------------------------------------------------
  # simulate multidimensional theta's
  #-------------------------------------------------
  
    if (dimension == 1){
      
      theta <- Unichange_sim(num_persons, sd_change)
      theta_pre <- theta[[1]]
      theta_post <- theta[[2]]
    
    } else{
    
      EMP <- FALSE
      theta <- Mulchange_sim(num_persons, dimension, cov_pretest, sd_change, EMP)
    
    }
  
    theta_pre <- theta[[1]]
    theta_post <- theta[[2]]
  
  #-------------------------------------------------------------------
  
  
    r_simresults <- matrix(NA, n_sim, num_methods)  # 8 methods.
  
    for (i in 1:n_sim){
    #-------------------------------------------------
    # simulate graded response data
    #-------------------------------------------------
    
      responses <- GRM_sim(theta_pre, itempar, id)
      response_pre <- responses[[1]]
      true_pre <- responses[[2]]
    
      responses <- GRM_sim(theta_post, itempar, id)
      response_post <- responses[[1]]
      true_post <- responses[[2]]
      
      carryover_results <- carry_over(response_pre, response_post)
      response_post_strong <- carryover_results[[1]]
      response_post_weak <- carryover_results[[2]]
   
      
      if (existence_carryover == 1){
        if (strong_weak == 1){
          response_post <- response_post_strong #replace with scores with strong carryover effects
        } else if (strong_weak == -1){
          response_post <- response_post_weak #replace with scores with weak carryover effects
        }
      }
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
    
      r_pre <- Lambda4::quant.lambda4(response_pre)$lambda4.quantile
      r_post <- Lambda4::quant.lambda4(response_post)$lambda4.quantile
      
      r_simresults[i, 5] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
    
      ######## method 2.1: estimated reliability (item-level) - alpha (i.e. change scores at item level are used to estimate reliability by means of alpha)   #########################
    
      response_change <- response_post - response_pre
      
      r_simresults[i, 6] <- psychometric::alpha(response_change)
    
      ######## method 2.2: estimated reliability (item-level) - lambda2 (i.e. change scores at item level are used to estimate reliability by means of lambda2) ################
    
      r_simresults[i, 7] <- Lambda4::lambda2(response_change)
    
      ########method 2.3: estimated reliability (item-level) - lambda4 (i.e. change scores at item level are used to estimate reliability by means of lambda4) 
     
      r_simresults[i, 8] <- Lambda4::quant.lambda4(response_change)$lambda4.quantile
    
    }
    
    sample_results[[p]] <- r_simresults
    
    r_avg[p, ] <- colSums(r_simresults)/n_sim
    p <- p+1
  
}
  
  
  
  


###################################################################################
#
#  A few plots
#
###################################################################################

##########################################################
# cross all samples of persons
results <- r_avg

plot(results[,1], results[, 2], xlab = "true reliability - proportion variance",
                     ylab = "true reliability - correlation", 
                     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 3], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by alpha )", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 4], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by lambda2)", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 6], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - alpha", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")


#######################
# one sample 

p <- sample(1:maxp, 1)

results_oneS <- sample_results[[p]]

plot(results_oneS[,1], results_oneS[, 2], xlab = "true reliability - proportion variance",
     ylab = "true reliability - correlation", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 3], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by alpha )", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 4], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by lambda2)", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 6], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - alpha", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")


######## trial #################

plot(results[, 2], results[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0.3, 0.5), ylim=c(0.3, 0.5),asp=1)
lines(c(0,1),c(0,1),col="red")

for (i in 1:maxp){
  
  r_02 <- sample_results[[i]][, 2]
  r_22 <- sample_results[[i]][, 7]
  
  points(r_02, r_22, pch=8, col = 'grey')
}

## some pretty colors
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

for (i in 1:maxp){
  
  r_02 <- sample_results[[i]][, 2]
  r_22 <- sample_results[[i]][, 7]
  
  ## compute 2D kernel density, see MASS book, pp. 130-131
  z <- kde2d(r_02, r_22, n=100)
  
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
}





