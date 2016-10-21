
#-----------------------------------------------------------------
#
#  A collection of funtions for psychometrics
#
#  Zhengguo Gu, Tilburg University
#
#-----------------------------------------------------------------

# Last update Sep. 2016

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
#-----------------------------------------------------------------

GRM_sim <- function(ability, itempar, id){
  
  # descrption:
  #
  # ability = ability parameter
  # itempar = item parameter
  # id = which theta goes with whitch discriminating parameter
  #
  # This function is designed for items with more than 2 answer alternatives!
  # (Thus, it needs to be extended to incorporate dichotomous items)
  
  n_sub <- nrow(ability)
  response <- matrix(NA, n_sub, nrow(itempar))
  true_response <- matrix(NA, n_sub, nrow(itempar))
  
  if(is.na(id[1])){

    for(i in 1:n_sub){
    
      numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
      P_star <- numeritor/(1+numeritor) # this is the "true response"
    
      response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
      true_response[i, ] <- rowSums(P_star)
    
    }
  } else{
    
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


#----------------------------------------------------------------
# MulChange_sim: simulate multidimensional change in theta
# Note that this is based on the code by Dr. Wilco Emons at TiU
#----------------------------------------------------------------

Mulchange_sim <- function(n_sub, dimension, covar, mean_change, sd_change, EMP){
  
  # generate theta values for pretest
  Sigma_theta <- matrix(covar, dimension, dimension)
  diag(Sigma_theta)  <- 1
  theta_pre <- mvrnorm(n_sub,mu=rep(0,dimension),Sigma=Sigma_theta,empirical=EMP)
  
  cov_change <- covar*sd_change^2 # We set it in this way, but there are of course other ways of setting up the covariance
  Sigma_change <- matrix(cov_change,dimension,dimension)
  diag(Sigma_change) <- sd_change^2
  theta_change <- mvrnorm(n_sub,mu=rep(mean_change,dimension),Sigma=Sigma_change,empirical=EMP)
  theta_post <- theta_pre + theta_change
  
  return(list(theta_pre, theta_post, theta_change))
  
}

#-------------------------------------------------------------
# UniChange_sim: simulate unidimensional change in theta
#-------------------------------------------------------------

Unichange_sim <- function(n_sub, sd_pre, mean_change, sd_change){
  # generate theta values for pretest
  theta_pre <- rnorm(n_sub, mean = 0, sd=1)

  sd_change <- sd_pre * sd_change^2  # need to justify why it is set like this
  theta_change <- rnorm(n_sub, mean_change, sd_change)
  
  theta_post <- theta_pre + theta_change
  
  return(list(theta_pre, theta_post, theta_change))
  
}
  