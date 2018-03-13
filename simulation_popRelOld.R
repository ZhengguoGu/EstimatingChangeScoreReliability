library(Lambda4)
library(psychometric)

#### This is to simulate population reliability


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
  
 
  ############ Calculate population reliability ##################################
  ### Here we use 1,000,000,000 people to approximate the population reliability
  
  #note that I planned to use Monte Carlo, but Wilco suggest that simulating a super 
  #large sample suffices. But I keep the parallel computing commands here (in case)
  
  #c2 <- makeCluster(12)
  #registerDoSNOW(c2)
  
  #set.seed(112)
  #sim_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
  
  if (dimension == 1){
    
    theta_pop <- Unichange_sim(1000000000, sd_change)
    theta_pre_pop <- theta_pop[[1]]
    theta_post_pop <- theta_pop[[2]]
    
  } else{
    
    EMP <- FALSE
    theta_pop <- Mulchange_sim(1000000000, dimension, cov_pretest, sd_change, EMP)
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



save(r_pop, file = "D:/ZG/PopulationRel.RData")