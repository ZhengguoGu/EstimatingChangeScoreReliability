
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
  