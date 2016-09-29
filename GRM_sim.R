# this is a program for simulating item reponse data from the graded response model (Samejima, 1969)

GRM_sim <- function(ability, itempar){
  
  n_sub <- length(ability)
  response <- matrix(NA, n_sub, nrow(itempar))
  
  compare_sub <- runif(n_sub, min=0, max=1)
  
  for(i in 1:n_sub){
    
    
    numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
    P_star<- numeritor/(1+numeritor)
    
    response[i, ] <- rowSums(P_star >= compare_sub[i])
  }
 
  
 return(response) 
  
}