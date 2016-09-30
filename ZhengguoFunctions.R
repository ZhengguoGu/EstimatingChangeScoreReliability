
#-----------------------------------------------------------------
#
#  A collection of funtions for psychometrics
#
#  Zhengguo Gu, Tilburg University
#
#-----------------------------------------------------------------

# Last update Sep. 2016


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
