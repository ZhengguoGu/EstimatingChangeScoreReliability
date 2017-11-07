#########################################################
#
# reliability of change scores at the population level
# via numerical integration
#
# Zhengguo Gu
#
##########################################################

library(mvtnorm)
library(MASS)
library(Matrix)
set.seed(110)

#------------------------------------------------------------------------------ 
# parameteres (they are directly copied from "simularion.R)
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

load("D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 3 - item difference scores\\20170929 itemPar/itemparLarge.RData")  #desktop at home
load("/Users/zhengguogu/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170929 itemPar/itemparLarge.RData")  #macbook pro
load("D:\\Dropbox\\Dropbox\\tilburg office\\Research Individual change\\Project 3 - item difference scores\\20170929 itemPar/itemparLarge.RData") #laptop in office


rel <- array()  #reliability at the population level

for(sim in 1: num_condition){
  cond <- df[sim, ]

  if(as.numeric(cond[3]) == 1){
    # in this case, unidimensional
    if(as.numeric(cond[5])==0){
      # in this case, no carry-over effects
      qpoint <- Qbipoints(n_person = 100, mu = c(0, 0), sigma = matrix(c(1, 0, 0, as.numeric(cond[4])^2), 2, 2), bd=3)  #simulate quadrature points for pretest and change
      variance_d <- 0  # variance of observed change scores over the entire distribution
      expectation_trueD <-0  #expectation of true change score
      expectation_trueD2 <- 0 ##expectation of true change score^2
        for(q in 1:length(qpoint$weights)){
          pre_theta <- qpoint$points[q,][1]
          post_theta <- sum(qpoint$points[q,])  #i.e., pretest + change
          
          pre_obs_prob <- GRM_sim_1theta(pre_theta, itempar = ITEM_PAR[[sim]])  #the probabilities of answering each category
          post_obs_prob <- GRM_sim_1theta(post_theta, itempar = ITEM_PAR[[sim]])
          
          pre_distribution <- Phi_X(pre_obs_prob)
          post_distribution <- Phi_X(post_obs_prob)
          
          d_distribution <- Phi_D(pre_distribution, post_distribution)
          
          expectation_d <- sum(d_distribution[1,] * d_distribution[2,])
          expectation_d2 <- sum(d_distribution[1,]^2 * d_distribution[2,])
          variance_d <- variance_d + (expectation_d2 - expectation_d^2) * qpoint$weights[q]
          
          expectation_trueD <- expectation_trueD + expectation_d * qpoint$weights[q]
          expectation_trueD2 <- expectation_trueD2 + expectation_d^2 * qpoint$weights[q]
        }
      variance_trueD <- expectation_trueD2 - expectation_trueD^2
      rel[sim] <- variance_trueD/variance_d 
    }else{
      #in this case, carry-over effect
      #we simulate 1m persons 
      theta_matrix <- rmvnorm(n = 1000000, mean = c(0, 0), sigma = matrix(c(1, 0, 0, as.numeric(cond[4]))^2, 2, 2))
      
      responses <- GRM_sim(theta_matrix[, 1], ITEM_PAR[[sim]], id = 1)
      response_pre <- responses[[1]]
      true_pre <- responses[[2]]
      
      responses <- GRM_sim(theta_matrix[, 2], ITEM_PAR[[sim]], id = 1)
      response_post <- responses[[1]]
      true_post <- responses[[2]]
      
      carryover_results <- carry_over(response_pre, response_post)
      response_post_strong <- carryover_results[[1]]
      response_post_weak <- carryover_results[[2]]
      
      if (as.matrix(cond[5]) == 10){
          response_post <- response_post_strong #replace with scores with strong carryover effects
      } else if (as.matrix(cond[5]) == 1){
          response_post <- response_post_weak #replace with scores with weak carryover effects
      }
      
      sum_pre <- rowSums(response_pre)
      sum_true_pre <- rowSums(true_pre)
      sum_post <- rowSums(response_post)
      sum_true_post <- rowSums(true_post)
      truechange_sumscores <- sum_true_post - sum_true_pre
      change_sumscores <- sum_post - sum_pre
      rel[sim] <- var(truechange_sumscores)/var(change_sumscores)
    }
  } else{
    # in this case, 3 dimensional
    dimension <- 3
    id <- vector()
    for(d in 1: dimension){
      id <- cbind(id, rep(d, as.numeric(cond[1])/dimension))
    }
    id <- as.vector(id)
    
     if(as.numeric(cond[5])==0){
      # in this case, no carry-over effects
       sigma_pre <- diag(1, 3)
       sigma_pre[which(sigma_pre==0)] <- as.numeric(cond[3])
       sigma_d <- diag(as.numeric(cond[4])^2, 3)
       sigma_d[which(sigma_d==0)] <- as.numeric(cond[4])^2 * as.numeric(cond[3])
       
       sigma_pre_d <- as.matrix(Matrix::bdiag(sigma_pre, sigma_d))
       qpoint <- mgauss.hermite(n = 10, mu = c(0, 0, 0, 0, 0, 0), sigma = sigma_pre_d, prune=.2)
       
       variance_d <- 0  # variance of observed change scores over the entire distribution
       expectation_trueD <-0  #expectation of true change score
       expectation_trueD2 <- 0 ##expectation of true change score^2
       
       for(q in 1:length(qpoint$weights)){
         pre_theta <- qpoint$points[q,][1:3]
         post_theta <- qpoint$points[q,][1:3] + qpoint$points[q,][4:6] #pretest + change
         
         pre_obs_prob <- GRM_sim_1theta(pre_theta, itempar = ITEM_PAR[[sim]], id = id)
         post_obs_prob <- GRM_sim_1theta(post_theta, itempar = ITEM_PAR[[sim]], id = id)
         
         pre_distribution <- Phi_X(pre_obs_prob)
         post_distribution <- Phi_X(post_obs_prob)
         
         d_distribution <- Phi_D(pre_distribution, post_distribution)
         
         expectation_d <- sum(d_distribution[1,] * d_distribution[2,])
         expectation_d2 <- sum(d_distribution[1,]^2 * d_distribution[2,])
         variance_d <- variance_d + (expectation_d2 - expectation_d^2) * qpoint$weights[q]
         
         expectation_trueD <- expectation_trueD + expectation_d * qpoint$weights[q]
         expectation_trueD2 <- expectation_trueD2 + expectation_d^2 * qpoint$weights[q]
       }
       
       variance_trueD <- expectation_trueD2 - expectation_trueD^2
       rel[sim] <- variance_trueD/variance_d 
       
     }else{
       #in this case, carry-over effect
       #we simulate 1m persons 
       
       if (as.numeric(cond[3]) == .6){
         cov_pretest <- .6
       } else if (as.numeric(cond[3]) == .8){
         cov_pretest <- .8
       }

       theta <- Mulchange_sim(num_persons=10000000, dimension = 3, cov_pretest, as.numeric(cond[4]), EMP = FALSE)
       
       theta_pre <- theta[[1]]
       theta_post <- theta[[2]]
       
       responses <- GRM_sim(theta_pre, ITEM_PAR[[sim]], id)
       response_pre <- responses[[1]]
       true_pre <- responses[[2]]
       
       responses <- GRM_sim(theta_post, ITEM_PAR[[sim]], id)
       response_post <- responses[[1]]
       true_post <- responses[[2]]
       
       carryover_results <- carry_over(response_pre, response_post)
       response_post_strong <- carryover_results[[1]]
       response_post_weak <- carryover_results[[2]]
       
       
       if (as.numeric(cond[5]) == 10){
           response_post <- response_post_strong #replace with scores with strong carryover effects
        } else if (as.numeric(cond[5]) == 1){
           response_post <- response_post_weak #replace with scores with weak carryover effects
        }
       
       
       
       sum_pre <- rowSums(response_pre)
       sum_true_pre <- rowSums(true_pre)
       sum_post <- rowSums(response_post)
       sum_true_post <- rowSums(true_post)
       truechange_sumscores <- sum_true_post - sum_true_pre
       change_sumscores <- sum_post - sum_pre
       rel[sim] <- var(truechange_sumscores)/var(change_sumscores)
     }
    
  }
  
  
}

