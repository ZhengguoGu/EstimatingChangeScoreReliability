library(Lambda4)
library(psychometric)
#library(foreach)
#library(doSNOW)
#library(doRNG)

#Update: 2017-11-13

#Note: On Jan 21, 2018, Zhengguo checked this file to add a few annotated descriptions. 

######################### START ###############################################################################################################

set.seed(110)

#------------------------------------------------------------------------------ 
# parameteres 
#------------------------------------------------------------------------------

test_length <- c(9, 21, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .1, .5) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
magnitude_sd <- c(sqrt(.14), sqrt(.5))  # .14 == small variance, .5 == large variance.
strongweak_carry <- c(0, 10, 1) #0 == no, 10 == strong, 1 == weak

num_condition <- length(test_length)*length(parallel_item)*length(correlated_facets)*length(magnitude_sd)*length(strongweak_carry)
conditions <- list()

p <- 1
for(i in 1:length(test_length)){
  for (j in 1:length(parallel_item)){
    for (k in 1:length(correlated_facets)){
      for(l in 1:length(magnitude_sd)){
        for(m in 1:length(strongweak_carry)){
          
          conditions[[p]] <- c(test_length[i], parallel_item[j], correlated_facets[k], magnitude_sd[l], strongweak_carry[m])
          p <- p + 1
        }
      }
    }
  }
}


df <- data.frame(matrix(unlist(conditions), nrow=num_condition, byrow = T))
#colnames(df) <- c('test length', 'parallel item', 'correlated facets',
                 # 'maginute of sd', 'carry-over effects')
#library(xlsx)
#library(rJava)
#write.csv(df, 'conditions.csv', sep=',')



# --------------------------------------------------------------------------------
#
# Run the simulations following the order of the test conditions --> df dataframe. 
#
#---------------------------------------------------------------------------------

restuls_conditions <- list()
simulatedRawdata <- list()  # this include the simulated theta's and responses. 
r_pop <- array()            # this is to save the population-level reliability. 
num_test <- 1

while (num_test <= nrow(df)){

  num_items <- df[num_test, 1]
  
  parallel_items <- df[num_test, 2]
  
  if(df[num_test, 3] == 1){
      dimension <- 1
  } else if (df[num_test, 3] == .1){
      dimension <- 3
      cov_pretest <- .1
  } else if (df[num_test, 3] == .5){
      dimension <- 3
      cov_pretest <- .5
  }

  sd_change <- df[num_test, 4]
  
  if(df[num_test, 5] == 0){
      existence_carryover <- 0
  } else if (df[num_test, 5] == 10){
      existence_carryover <- 1
      strong_weak <- 1  # 1: strong carry-over effect; Please dont be confused: Later a function carry_over() will be called, and this function recognize 1 as strong effect. 
                        # I know that this setup (i.e., first coded as 10 and then swtich to 1) is a bit tedious. 
  } else if (df[num_test, 5] == 1){
      existence_carryover <- 1
      strong_weak <- -1  # -1: weak carry-over effect (please see the comment above.)
  }



    #------------------------------------------------------------------------------

  if (parallel_items == 1) {
  
    itempar <- matrix(NA,num_items,5)
    itempar[,1] <- runif(1,1.5,2.5)   # discrimination
    avg_beta <- runif(1, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
  
  } else {
  
    itempar <- matrix(NA,num_items,5)
    itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
    avg_beta <- runif(num_items, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
  
  }

  #ITEM_PAR[[num_test]] <- itempar
  
  id <- vector()
  for(d in 1: dimension){
    id <- cbind(id, rep(d, num_items/dimension))
  }
  id <- as.vector(id)  #this is for identifying the facets/dimensions in theta. 
  
  #####################################################
  #
  #  generate random samples from population (theta)
  #
  #####################################################

  
  maxp <- 20  # thus 20 samples of persons are simulated from the population
  num_persons <- 1000 # number of persons
  n_sim <- 50 # each sample generates 50 item-score datasets
  num_methods <- 7 #see below, method 0.1 to 2.3
  r_avg <- matrix(NA, maxp, num_methods)
  r_sd <- matrix(NA, maxp, num_methods)
  sample_results <- list()
  sample_theta <- list()
  simResponses <- list()  #this contains a list of p lists. 
  
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
      theta_pre <- theta[[1]]
      theta_post <- theta[[2]]
      
    }
  
    
    sample_theta[[p]] <- list(theta_pre, theta_post)
  #-------------------------------------------------------------------
  
  
    r_simresults <- matrix(NA, n_sim, num_methods)  # 7 methods.
  
    
    sumpre_response <- list()
    sum_Truepre_response <- list()
    sumpost_response <- list()
    sum_Truepost_response <- list()
    
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
      
      #post_response_carry[[i]] <- response_post 
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
      sumpre_response[[i]] <- sum_pre
      sum_true_pre <- rowSums(true_pre)
      sum_Truepre_response[[i]] <- sum_true_pre
      sum_post <- rowSums(response_post)
      sumpost_response[[i]] <- sum_post
      sum_true_post <- rowSums(true_post)
      sum_Truepost_response[[i]] <- sum_true_post
      #####################################################################################################################################################
      # Methods for calculating reliability
      #
      # method 0.1: ture-change reliability - cor(true change, observed change)^2  (Note, because later we approximate the true change-score reliability 
      #                                                                             in the population, method 0.1 is not used in the paper. But it is nice
      #                                                                             to check it.)
      # 
      # method 1.1: estimated reliability - alpha (i.e. pre and post reliability estimated by alpha )
      # method 1.2: estimated reliability - lambda2 (i.e. pre and post reliability estimated by lambda2)
      # method 1.3: estimated reliability - lambda4 (i.e. pre and post reliability estimated by lambda4)
      #
      # method 2.1: estimated reliability (item-level) - alpha (i.e. change scores at item level are used to estimate reliability by means of alpha) 
      # method 2.2: estimated reliability (item-level) - lambda2 (i.e. change scores at item level are used to estimate reliability by means of lambda2) 
      # method 2.3: estimated reliability (item-level) - lambda4 (i.e. change scores at item level are used to estimate reliability by means of lambda4) 
      ###################################################################################################################################################
    

    
      ######## method 0.1: ture-change reliability - cor(true change, observed change)^2 ###########
      
      truechange_sumscores <- sum_true_post - sum_true_pre
      change_sumscores <- sum_post - sum_pre
    
      r_simresults[i, 1] <- (cor(truechange_sumscores, change_sumscores))^2
    
      ######## method 1.1: estimated reliability - alpha (i.e. pre and post reliability estimated by alpha )  #########################
    
      r_pre <- psychometric::alpha(response_pre)  # ! cronback alpha is used here. 
      r_post <- psychometric::alpha(response_post)
    
      r_simresults[i, 2] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
      
      ######## method 1.2: estimated reliability - lambda2 (i.e. pre and post reliability estimated by lambda2) #######
    
      r_pre <- Lambda4::lambda2(response_pre)  # ! lambda2 is used here. 
      r_post <- Lambda4::lambda2(response_post)
    
      r_simresults[i, 3] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
    
      ######## method 1.3: estimated reliability - lambda4 (i.e. pre and post reliability estimated by lambda4)
    
      r_pre <- Lambda4::quant.lambda4(response_pre)$lambda4.quantile
      r_post <- Lambda4::quant.lambda4(response_post)$lambda4.quantile
      
      r_simresults[i, 4] <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))
    
      ######## method 2.1: estimated reliability (item-level) - alpha (i.e. change scores at item level are used to estimate reliability by means of alpha)   #########################
    
      response_change <- response_post - response_pre
      
      r_simresults[i, 5] <- psychometric::alpha(response_change)
    
      ######## method 2.2: estimated reliability (item-level) - lambda2 (i.e. change scores at item level are used to estimate reliability by means of lambda2) ################
    
      r_simresults[i, 6] <- Lambda4::lambda2(response_change)
    
      ########method 2.3: estimated reliability (item-level) - lambda4 (i.e. change scores at item level are used to estimate reliability by means of lambda4) 
     
      r_simresults[i, 7] <- Lambda4::quant.lambda4(response_change)$lambda4.quantile
    
    }
    
    sample_results[[p]] <- r_simresults
    
    r_avg[p, ] <- colSums(r_simresults)/n_sim
    r_sd[p, ] <- apply(r_simresults, 2, sd)
    
    simResponses[[p]] <- list(sumpre_response, sum_Truepre_response, sumpost_response, sum_Truepost_response) 
    
    p <- p+1
    
    
  }
  
  restuls_conditions[[num_test]] <- list(sample_results, r_avg, r_sd)
  simulatedRawdata[[num_test]] <- list(sample_theta, simResponses)
  
  
  ############ Calculate population reliability ##################################
  ### Here we use 1 million people to approximate the population reliability
  
  #note that I planned to use Monte Carlo, but Wilco suggest that simulating a super 
  #large sample suffices. But I keep the parallel computing commands here (in case)
  
  #c2 <- makeCluster(12)
  #registerDoSNOW(c2)
  
  #set.seed(112)
  #sim_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
   
  if (dimension == 1){
    
    theta_pop <- Unichange_sim(1000000, sd_change)
    theta_pre_pop <- theta_pop[[1]]
    theta_post_pop <- theta_pop[[2]]
    
  } else{
      
    EMP <- FALSE
    theta_pop <- Mulchange_sim(1000000, dimension, cov_pretest, sd_change, EMP)
    theta_pre_pop <- theta_pop[[1]]
    theta_post_pop <- theta_pop[[2]]
  }
  
  
  #-------------------------------------------------
  # simulate graded response data
  #-------------------------------------------------
  
  responses_pop <- GRM_sim(theta_pre_pop, itempar, id)
  response_pre_pop <- responses_pop[[1]]
 
  true_pre_pop <- responses_pop[[2]]
   
  responses_pop <- GRM_sim(theta_post_pop, itempar, id)
  response_post_pop <- responses_pop[[1]]
 
  true_post_pop <- responses_pop[[2]]
  
  carryover_results_pop <- carry_over(response_pre_pop, response_post_pop)
  response_post_strong_pop <- carryover_results_pop[[1]]
  response_post_weak_pop <- carryover_results_pop[[2]]
  
   
  if (existence_carryover == 1){
    if (strong_weak == 1){
      response_post_pop <- response_post_strong_pop #replace with scores with strong carryover effects
    } else if (strong_weak == -1){
      response_post_pop <- response_post_weak_pop #replace with scores with weak carryover effects
    }
  }
  
  
    # sum scores
  
  sum_pre_pop <- rowSums(response_pre_pop)
  sum_true_pre_pop <- rowSums(true_pre_pop)
  sum_post_pop <- rowSums(response_post_pop)
  sum_true_post_pop <- rowSums(true_post_pop)
   
  truechange_sumscores_pop <- sum_true_post_pop - sum_true_pre_pop
  change_sumscores_pop <- sum_post_pop - sum_pre_pop
 
  #r_population <- (cor(truechange_sumscores_pop, change_sumscores_pop))^2
  
    #return(r_population)
  #}
  #stopCluster(c2)
  
  #r_pop[num_test] <- mean(sim_result)
  r_pop[num_test] <- (cor(truechange_sumscores_pop, change_sumscores_pop))^2
  
  num_test <- num_test + 1
  
  print(num_test)
} # END OF WHILE



save(r_pop, restuls_conditions, simulatedRawdata, file = "D:/ZG/LargeSample20171122.RData")

