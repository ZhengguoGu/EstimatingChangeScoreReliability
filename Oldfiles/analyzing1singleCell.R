## simulate one cell in the simulation ###
# Note: This is for diagnosing the simulations by taking a closer look at, for example, what kinda
# responses have been generated, what the estimated reliability is like. 




#------------------------------------------------------------------------------ 
# parameteres 
#------------------------------------------------------------------------------

test_length <- c(9, 21, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .6, .8) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
magnitude_sd <- c(sqrt(.14), sqrt(.5))
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
colnames(df) <- c('test length', 'parallel item', 'correlated facets',
                  'maginute of sd', 'carry-over effects')

#-----------------------------------------------


num_test <- 17  #this is the 17th low of the df file (which contains 108 rows, thus 108 cells)
  
num_items <- df[num_test, 1]
  
parallel_items <- df[num_test, 2]
  
if(df[num_test, 3] == 1){
  dimension <- 1
} else if (df[num_test, 3] == .6){
  dimension <- 3
  cov_pretest <- .6
} else if (df[num_test, 3] == .8){
  dimension <- 3
  cov_pretest <- .8
}
  
sd_change <- df[num_test, 4]
  
if(df[num_test, 5] == 0){
  existence_carryover <- 0
} else if (df[num_test, 5] == 10){
  existence_carryover <- 1
  strong_weak <- 1  # 1: strong !!!note that this part might be a bit confusing
} else if (df[num_test, 5] == 1){
  existence_carryover <- 1
  strong_weak <- -1  # -1: weak !!!note that this part might be a bit confusing
}
  
  
  
  
#------------------------------------------------------------------------------
  
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
id <- as.vector(id)  #this is for identifying the facets/dimensions in theta. 
  
  #####################################################
  #
  #  generate random samples from population (theta)
  #
  #####################################################
  
  
maxp <- 1  # thus 100 samples (each with 1000 persons) were simulated
num_persons <- 1000 # number of subjects
n_sim <- 5 # simulate 1000 datasets
num_methods <- 8 #see below, method 0.1 to 2.3
r_avg <- matrix(NA, maxp, num_methods)
r_sd <- matrix(NA, maxp, num_methods)
sample_results <- list()
simulated_data <- list() #this contains simulated responses etc.
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
  pre_response <- list()
  pre_response_true <- list()
  post_response <- list()
  post_response_true <- list()
  post_response_carry <- list()
  pre_sum <- list()
  pre_sum_true <- list()
  post_sum <- list()
  post_sum_true <- list()
  
  for (i in 1:n_sim){
    #-------------------------------------------------
    # simulate graded response data
    #-------------------------------------------------
      
    responses <- GRM_sim(theta_pre, itempar, id)
    response_pre <- responses[[1]]
    true_pre <- responses[[2]]
    pre_response[[i]] <- response_pre
    pre_response_true[[i]] <- true_pre
      
    responses <- GRM_sim(theta_post, itempar, id)
    response_post <- responses[[1]]
    true_post <- responses[[2]]
    post_response[[i]] <- response_post
    post_response_true[[i]] <- true_post
    
      
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
    
    post_response_carry[[i]] <- response_post
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
    pre_sum[[i]] <- sum_pre
    sum_true_pre <- rowSums(true_pre)
    pre_sum_true[[i]] <- sum_true_pre
    sum_post <- rowSums(response_post)
    post_sum[[i]] <- sum_post
    sum_true_post <- rowSums(true_post)
    post_sum_true[[i]] <- sum_true_post
      
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
  simulated_data[[p]] <- list(pre_response, pre_response_true, post_response, post_response_true, post_response_carry, 
                              pre_sum, pre_sum_true, post_sum, post_sum_true)
  
  r_avg[p, ] <- colSums(r_simresults)/n_sim
  r_sd[p, ] <- apply(r_simresults, 2, sd)
  p <- p+1
    
}

restuls_conditions <- list(sample_results, r_avg, r_sd, simulated_data)

#### Let's see what happened to each simulated sample of responses (note here i simulated one sample of
# persons and they are fixed).

sample_n <- 1  #note sample_n = 1,2,..., 5. Because n_sim=5

res_pre <- restuls_conditions[[4]][[1]][[1]][[sample_n]]
res_pre_true <- restuls_conditions[[4]][[1]][[2]][[sample_n]]
res_post <- restuls_conditions[[4]][[1]][[3]][[sample_n]]
res_post_true <- restuls_conditions[[4]][[1]][[4]][[sample_n]]
res_post_carry <- restuls_conditions[[4]][[1]][[5]][[sample_n]]
sum_pre <- restuls_conditions[[4]][[1]][[6]][[sample_n]]
sum_pre_true <- restuls_conditions[[4]][[1]][[7]][[sample_n]]
sum_post <- restuls_conditions[[4]][[1]][[8]][[sample_n]]
sum_post_true <- restuls_conditions[[4]][[1]][[9]][[sample_n]]

######## method 1.1: estimated reliability - alpha (i.e. pre and post reliability estimated by alpha )  #########################

r_pre <- psychometric::alpha(res_pre)  # ! cronback alpha is used here. 
r_post <- psychometric::alpha(res_post_carry)

r_est1<- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))

######## method 0.1: ture-change reliability - var(true change)/var(observed change) ###########

truechange_sumscores <- sum_post_true - sum_pre_true 
change_sumscores <- sum_post - sum_pre

r_est2<- var(truechange_sumscores)/var(change_sumscores)


#############--------------------------------------
# given the 108 cells, which cells generate negative reliabilites?

weird_index <- rep(0, 108)
for(i in 1:108){
  allestimates <- matrix(unlist(restuls_conditions[[i]][2]), nrow = 20, ncol = 8) 
  if(sum(allestimates < 0)>=1){
    weird_index[i] <- 1
  }
}

weird_index2 <- rep(0, 108)
for(i in 1:108){
  allestimates <- matrix(unlist(restuls_conditions[[i]][2]), nrow = 20, ncol = 8) 
  if(sum(allestimates[, 2:8] >1)>=1){
    weird_index2[i] <- 1
  }
}
