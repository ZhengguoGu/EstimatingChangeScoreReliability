
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
#-----------------------------------------------------------------

GRM_sim <- function(ability, itempar){
  
  # descrption:
  # ability = ability parameter, aka, theta
  # itempar = item parameter. Note that the first column is discrimination (alpha)
  #
  # This function is designed for items with more than 2 answer alternatives!
  # (Thus, it needs to be extended to incorporate dichotomous items)

  n_sub <- length(ability)
  response <- matrix(NA, n_sub, nrow(itempar))


  for(i in 1:n_sub){

    numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
    P_star<- numeritor/(1+numeritor)

    response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
  }

  return(response)
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
  