# this is a program for simulating item reponse data from the graded response model (Samejima, 1969)

GRM_sim <- function(ability, itempar){
  
  n_sub <- length(ability)
  #P <- matrix(NA, nrow(itempar), ncol(itempar))
  response <- matrix(NA, n_sub, nrow(itempar))
  
  for(i in 1:n_sub){
    
    P <- matrix(NA, nrow(itempar), ncol(itempar))
    
    
    numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
    P_star<- numeritor/(1+numeritor)
    
    P[, ncol(P_star)+1] <- P_star[,ncol(P_star)]
    
    for(k in ncol(P_star):2){
      P[,k] <- P_star[,k-1] - P_star[,k]
    }
    
    P[, 1] <- 1 - P_star[, 1]
    
    for(j in 1:nrow(itempar)) {
      response[i, j] <- which(P[j,] == max(P[j,]))
    }
  }
 
 return(response) 
  
}