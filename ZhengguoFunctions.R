
#-----------------------------------------------------------------
#
#  A collection of funtions for psychometrics
#
#  Zhengguo Gu, Tilburg University
#
#-----------------------------------------------------------------

# Last check ???????????????????????????


library(MASS)

#-----------------------------------------------------------------
# GRM_sim: simulate response data based on the graded response
# model (Samejima, 1969)
# 
# 2016.10.14. Updated. An "Id" variable is added, so as to 
# generate multidimensional responses. (Based on Wilco's code)
# However, it should be noted that this is not strickly a typical
# multidimensional irt model. Here I assume that items have 
# different facets - that is, a person might have, say 3 dimenions 
# in theta, and each dimension persons to a few items. But there
# is no crossloadings between items and thetas
#
#-----------------------------------------------------------------

GRM_sim <- function(ability, itempar, id){
  
  # descrption:
  #
  # ability = ability parameter
  # itempar = item parameter
  # id = which theta goes with which discriminating parameter
  #
  # This function is designed for items with more than 2 answer alternatives!
  
  
  if(sum(id == 1)==length(id)) {
    
    n_sub <- length(ability)
    response <- matrix(NA, n_sub, nrow(itempar))
    true_response <- matrix(NA, n_sub, nrow(itempar))
    
    # if true, then unidimensional 
    for(i in 1:n_sub){
    
      numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
      P_star <- numeritor/(1+numeritor) # this is the "true response"
    
      response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
      true_response[i, ] <- rowSums(P_star)
    
    }
  } else{
    
    n_sub <- nrow(ability)
    response <- matrix(NA, n_sub, nrow(itempar))
    true_response <- matrix(NA, n_sub, nrow(itempar))
    
    for(i in 1:n_sub){
      
      abil <- ability[i, ][id]
      numeritor <- exp(sweep((abil-itempar[, -1]), 1, itempar[, 1], "*"))
      P_star <- numeritor/(1+numeritor) # this is the "true response"
    
      response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
      true_response[i, ] <- rowSums(P_star)
    }
  }  
  
  return(list(response, true_response))
}

#-----------------------------------------------------------------
# GRM_sim_1theta: simulate response data based on the graded response
# model (Samejima, 1969), given 1 theta
#-----------------------------------------------------------------

GRM_sim_1theta <- function(theta, itempar){
  
  # description:
  #
  # theta: person parameter; a scaler 
  # itempar = item parameter
  #
  # This function is designed for items with more than 2 answer alternatives!
  
  true_response <- matrix(NA, 1, nrow(itempar))
  
  numeritor <- exp(itempar[, 1]*(theta - itempar[, -1]))
  P_star <- numeritor/(1+numeritor) # this is the "true response"
  P <- cbind(1, P_star)
  P[, 1:4] <- P[, 1:4] - P_star    
  return(P)  #returns the probability of answering each category 
}


#----------------------------------------------------------------
# MulChange_sim: simulate multidimensional change in theta
# Note that this is based on the code by Dr. Wilco Emons at TiU
#
# Last update: 2016/12/08
#----------------------------------------------------------------

Mulchange_sim <- function(n_sub, dimension, covar, sd_change, EMP){
  
  # generate theta values for pretest
  Sigma_theta <- matrix(covar, dimension, dimension)
  diag(Sigma_theta)  <- 1
  theta_pre <- mvrnorm(n_sub,mu=rep(0,dimension),Sigma=Sigma_theta,empirical=EMP)
  
  cov_change <- covar*sd_change^2 # We set it in this way, but there are of course other ways of setting up the covariance
  Sigma_change <- matrix(cov_change,dimension,dimension)
  diag(Sigma_change) <- sd_change^2
  theta_change <- mvrnorm(n_sub,mu=rep(0,dimension),Sigma=Sigma_change,empirical=EMP)
  theta_post <- theta_pre + theta_change
  
  return(list(theta_pre, theta_post, theta_change))
  
}

#-------------------------------------------------------------
# UniChange_sim: simulate unidimensional change in theta
#
# Last update: 2016/12/08
#-------------------------------------------------------------

Unichange_sim <- function(n_sub, sd_change){
  # generate theta values for pretest
  
  theta_pre <- rnorm(n_sub, mean = 0, sd=1)
  theta_change <- rnorm(n_sub, mean = 0, sd = sd_change) #note that because mean does not influnce reliability
                                                    #we set mean = 0
  theta_post <- theta_pre + theta_change
  
  return(list(theta_pre, theta_post, theta_change))
  
}


#--------------------------------------------------------------
# Mimicking carry-over effects
#
# Last update: 2016/12/08
#--------------------------------------------------------------

carry_over <- function(pre, post){
  
  strong_post <- post
  weak_post <- post  
  
  ind1 <- (pre - post < -1) 
  strong_post[ind1] <- pre[ind1] + 1
  weak_post[ind1] <- post[ind1] - 1
  
  ind2 <- (pre - post > 1)
  strong_post[ind2] <- pre[ind2] - 1
  weak_post[ind2] <- post[ind2] + 1
  
  ind3 <- (abs(pre-post)==1 )
  strong_post[ind3] <- pre[ind3]
  weak_post[ind3] <- post[ind3]
  

  
  return(list(strong_post, weak_post))
}
  

#----------------------------------------------------------------------
# Recursive algorithm for expanding Phi(X+|theta) for polytomous items,
# which include GRM_sim_1theta, 
#
# Last update 2017-9-29
#----------------------------------------------------------------------

Phi_X <- function(P){
  
  # description:
  # 
  # P: the probabilities of answering each category per item correctly. 
  #    It is a matrix, with rows representing items, and columns representing
  #    answer categories.

  # initial step:
  p1 <- P[1, ]
  p2 <- P[2, ]
  cate_P <- 0:(dim(P)[2] - 1) #answer categories: 0, 1, 2, ....
  
  p_expand <- expand.grid(p1, p2)
  p_expand <- p_expand[, 1] * p_expand[, 2]
  cate_expand <- rowSums(expand.grid(cate_P, cate_P))
  w <- matrix(NA, 1, max(cate_expand)+1)
  for(i in 1:(max(cate_expand)+1)){
    w[1, i] <- sum(p_expand[cate_expand == (i-1)])
  }
    
  for(i in 3:dim(P)[1]){
    p_i <- P[i, ]
    p_expand <- expand.grid(p_i, w)
    p_expand <- p_expand[, 1] * p_expand[, 2]
    cate_expand <- rowSums(expand.grid(cate_P, 0:(length(w) - 1)))
    w <- matrix(NA, 1, max(cate_expand)+1)
    for(j in 1:(max(cate_expand)+1)){
      w[1, j] <- sum(p_expand[cate_expand == (j-1)])
    }
  }
  
  return(w) #this is the probability of answering X+ = 0, 1, 2, 3, 4, ..... 
}

#----------------------------------------------------------------------
# Calculat the full distribution of change scores. 
#
# Last update 2017-9-29
#----------------------------------------------------------------------

Phi_D <- function(Prob_pre, Prob_post){
  
  # Description
  #
  # Prob_pre: the probabilities of answering X+ = 0, 1, 2, ..., of a person at pretest
  # Prob_post: the probabilities of answering X+ = 0, 1, 2, ..., of a person at posttest
  
  
  cate_pre <- 0:(length(Prob_pre)-1)
  cate_post <- 0:(length(Prob_post)-1)
  
  cate_expand <- expand.grid(cate_pre, cate_post)
  D_scores <- cate_expand[, 2] - cate_expand[, 1]
  prob_expand <- expand.grid(Prob_pre, Prob_post)
  prob_D <- prob_expand[, 1] * prob_expand[, 2]
  
  prob_Dfinal <- matrix(NA, 1, length(unique(D_scores)))
  D_index <- sort(unique(D_scores)) #change score sorted from min to max
  for(j in 1:dim(prob_Dfinal)[2]){
    prob_Dfinal[1, j] <- sum(prob_D[D_scores == D_index[j]]) # probabilities corresponding to D_index
  }
  return(D_index, prob_Dfinal)
}


