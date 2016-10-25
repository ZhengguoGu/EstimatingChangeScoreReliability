library(Lambda4)
library(psychometric)

# data simulation

################################################

set.seed(110)

#------------------------------------------------- 
# simulate item parameters
#-------------------------------------------------

parallel_items <- 2 # 1: parallel, otherwise nonparallel
dimension <- 4 # number of dimensions in theta
num_persons <- 1000

num_items <- 20

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
maxp <- 20
num_methods <- 8 #see below, method 0.1 to 2.3
r_avg <- matrix(NA, maxp, num_methods)
r_sd <- matrix(NA, maxp, num_methods)
p <- 1

while(p<=maxp) {

  #-------------------------------------------------
  # simulate multidimensional theta's
  #-------------------------------------------------
  
    if (dimension == 1){
    
    
      n_sub <- num_persons
    
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
    
      theta <- Mulchange_sim(num_persons, dimension, cov_pretest, mean_change, sd_change, EMP)
    
    }
  
    theta_pre <- theta[[1]]
    theta_post <- theta[[2]]
  
  #-------------------------------------------------------------------
  
    n_sim <- 1000 # simulate 1000 datasets
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
    
    sample_results[[p]] <- r_simresults
    
    r_avg[p, ] <- colSums(r_simresults)/n_sim
    r_sd[p, ] <- apply(r_simresults, 2, sd)
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





