
#-----------------------------------------------------------------
#
#  A collection of funtions for psychometrics
#
#  (Note: Dear reviewers, not all functions in this documents
# are relevant. For example, the numerical integration for estimating
# the reliability: later on we realize that for the carry-over effects
# conditions, numerical integration doesn't work; we mention this also
# in the paper.)
#-----------------------------------------------------------------


library(MASS)
library(mvtnorm)
#-----------------------------------------------------------------
# GRM_sim: simulate response data based on the graded response
# model (Samejima, 1969)
# 
# 2016.10.14. Updated. An "Id" variable is added, so as to 
# generate multidimensional responses. (Based on the second author's code)
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

GRM_sim_1theta <- function(theta, itempar, id){
  
  # description:
  #
  # theta: person parameter; a scaler 
  # itempar = item parameter
  #
  # This function is designed for items with more than 2 answer alternatives!
  
  theta <- as.numeric(theta) 
  if(missing(id)){
    
    true_response <- matrix(NA, 1, nrow(itempar))
    
    numeritor <- exp(itempar[, 1]*(theta - itempar[, -1]))
    P_star <- numeritor/(1+numeritor) # this is the "true response"
    P <- cbind(1, P_star)
    P[, 1:4] <- P[, 1:4] - P_star  
    
  
    }else{
    
    theta <- theta[id]
    numeritor <- exp(sweep((theta-itempar[, -1]), 1, itempar[, 1], "*"))
    P_star <- numeritor/(1+numeritor) # this is the "true response"
    
    P <- cbind(1, P_star)
    P[, 1:4] <- P[, 1:4] - P_star
  
  }
    
  return(P)  #returns the probability of answering each category 
}


#----------------------------------------------------------------
# MulChange_sim: simulate multidimensional change in theta
# Note that this is based on the code by the second author
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
# Recursive algorithm for expanding Phi(X+|theta) for polytomous items
# 
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
  # Prob_pre: the probabilities of answering X+ = 0, 1, 2, ..., of a person at pretest, 
  #           which can be obtained by Phi_X()
  # Prob_post: the probabilities of answering X+ = 0, 1, 2, ..., of a person at posttest,
  #           which can be obtained by Phi_X()
  
  cate_pre <- 0:(length(Prob_pre)-1)
  cate_post <- 0:(length(Prob_post)-1)
  cate_expand <- expand.grid(cate_pre, cate_post)
  Dtrue_scores <- cate_expand[, 2] - cate_expand[, 1] 
  
  prob_expand <- expand.grid(Prob_pre, Prob_post)
  prob_D <- prob_expand[, 1] * prob_expand[, 2]
  
  
  prob_DtrueFinal <- matrix(NA, 1, length(unique(Dtrue_scores)))
  D_index <- sort(unique(Dtrue_scores))
  for(j in 1:dim(prob_DtrueFinal)[2]){
    prob_DtrueFinal[1, j] <- sum(prob_D[Dtrue_scores == D_index[j]]) # probabilities corresponding to D_index
  }
  
  D_true <- rbind(D_index, prob_DtrueFinal)  
  row.names(D_true) <- c("score", "probability")

  return(D_true)
}

#-----------------------------------------------------------------------
# quadrature points, gaussian 
#-----------------------------------------------------------------------
Qpoints <- function(n_person, mu, sd, bd){
  
  # description:
  #
  # n_person: number of thetas to be simulated
  # mu: the average of theta distribution
  # sd: the sd of theta distribution
  # bd: for defining the upper and lower bound of theta: 
  #     upper bound: mu + bd*sd
  #     lower bound: mu - bd*sd
  
  q_point <- seq(mu - bd*sd, mu + bd*sd, length.out = n_person)
  q_weight <- dnorm(q_point, mean = mu, sd = sd)*abs(q_point[1]-q_point[2])/sum(dnorm(q_point, mean = mu, sd = sd)*abs(q_point[1]-q_point[2]))

  q_result <- list(q_point, q_weight)
  names(q_result) <- c("points", "weights")
  return(q_result)
}



#---------------------------------------------------------------------
# quadrature points, bivariate gaussian (for pretest and change) 
# based on Rieman Rule
#---------------------------------------------------------------------
Qbipoints <- function(n_person, mu, sigma, bd){
  # description:
  # requires package 'mvtnorm'
  #
  # n_person: number of theta's
  # mu: the mean vector
  # sigma: the covariance matrix (between pretest and change); Note: covariance between pretest and change == 0!
  # bd: bound. The the upper/lower bound is defind along the first dimension 
  #     is mu[1] +/- sqrt(sigma[1,1])*bd; and for the second dimension
  #     mu[2] +/- sqrt(sigma[2,2])*bd
  
  if(length(mu) != 2){
    stop("This fuction is for bivariate mvtnormal!")
  }
  q_point <- matrix(NA, n_person, 2)
  
  for(i in 1:2){
    q_point[, i] <- seq(mu[i] - sqrt(sigma[i,i])*bd, mu[i] + sqrt(sigma[i,i])*bd, length.out = n_person)
  }
  
  q_point_final <- expand.grid(q_point[, 1], q_point[, 2])
  
  densities <- dmvnorm(q_point_final, mean = mu, sigma = sigma)
  q_weight <- densities*abs(q_point[1, 1]-q_point[2, 1])*abs(q_point[1, 2]-q_point[2, 2])/sum(densities*abs(q_point[1, 1]-q_point[2, 1])*abs(q_point[1, 2]-q_point[2, 2]))
  
  q_result <- list(q_point_final, q_weight)
  names(q_result) <- c("points", "weights")
  
  return(q_result)
  
}

#-----------------------------------------------------------------------------
# Gauss-Hermite quadrature
# The following code is obtained from https://biostatmatt.com/archives/2754
#
#------------------------------------------------------------------------------

hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}

# compute multivariate Gaussian quadrature points
## n     - number of points each dimension before pruning
## mu    - mean vector
## sigma - covariance matrix
## prune - NULL - no pruning; [0-1] - fraction to prune
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}
