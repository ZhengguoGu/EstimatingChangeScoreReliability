##########################################################################################
#
# This file contains the code for analyzing the results
#
##########################################################################################

# 1. load simulation results
load('results20170122.RData')

# 2. Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells 
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:8] 

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> reliability estimates of 20 samples from persons under condition cell 1 
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8]

# 3. analysis at the population level: that is, we average reliability estimates across the 20 samples of persons.
# The idea is that the estimates reflect the reliability at the population level. 

pop_re <- matrix(NA, nrow = 108, ncol = 8)
for (i in 1:108){
  pop_re[i, ] <- colMeans(restuls_conditions[[i]][[2]])
}

# plot pop_re
plot(pop_re[, 2], xlim = c(0, 108), ylim = c(0,1), xlab = "108 cells", ylab = "True Reliability", col='red', pch=19)

plot(pop_re[, 2], xlim = c(0, 108), ylim = c(0,1), xlab = "108 cells", ylab = "(Estimated) Reliability", col='red', pch=19)
points(pop_re[, 3], pch=0)
points(pop_re[, 4], pch=1)
points(pop_re[, 5], pch=2)
points(pop_re[, 6], pch=15, col='blue')
points(pop_re[, 7], pch=19, col='blue')
points(pop_re[, 8], pch=17, col='blue')

for(i in 1:108){
  abline(v=i, lty="dotted")
}
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

df[, 6] <- weird_index

########################Below: not useful##########################

# explore between-person and within-person variance for each cell

# explore the simulated raw data
#(simulatedRawdata)  # 108 cells
#length(simulatedRawdata[[1]]) # 2 lists sample_theta, simResponses
#length(simulatedRawdata[[1]][[1]]) # 20 samples of persons
#length(simulatedRawdata[[1]][[2]]) # 20 samples of persons

# Focus on simulated responses now
#length(simulatedRawdata[[1]][[2]][[1]]) #4 lists, given the 1st cell, 2nd list (i.e. simResponses), and 1st sample (out of 20 samples of persons)
#length(simulatedRawdata[[1]][[2]][[1]][[1]]) # 50 response datasets of sumpre_response, given the 1st cell, 2nd list (i.e. simResponses), and 1st sample (out of 20 samples of persons)
#length(simulatedRawdata[[1]][[2]][[1]][[1]][[50]]) #sumpre_response data (the 50th dataset), given the 1st cell, 2nd list (i.e. simResponses), and 1st sample (out of 20 samples of persons)

# For each sample of persons, 50 simulated responses are generated. --> these persons are measured 50 times with brain washing. 
# Now we generate a matrix, the rows are the persons, the columns are the 50 change scores 


#BTinform <- list() #this list contains the between/within variances of all the 108 cases. 
#for (c in 1:108){
#  Matrixlist <- list()
#  varDecomp <- matrix(NA, 20, 4)
#  for (p in 1:20){
#    wbMatrix <- matrix(0, 1000, 50)
#    for (s in 1:50){
#      wbMatrix[, s] <- simulatedRawdata[[c]][[2]][[p]][[3]][[s]] - simulatedRawdata[[c]][[2]][[p]][[1]][[s]]
#    }
#    Matrixlist[[p]]<- wbMatrix
#    grandmean <- mean(wbMatrix)
#    personmean <- rowMeans(wbMatrix)
#    totalVar <- sum((wbMatrix-grandmean)^2)
#    betweenPVar <- sum(50*(personmean-grandmean)^2)
#    withinPVar <- sum((wbMatrix-personmean)^2)
#    btRatio <- betweenPVar/totalVar 
#    varDecomp[p, 1] <- totalVar
#    varDecomp[p, 2] <- betweenPVar
#    varDecomp[p, 3] <- withinPVar
#    varDecomp[p, 4] <- btRatio
#  }
#  BTinform[[c]] <- varDecomp
#}

#since each cell contains 20 samples of persons, we average them per cell
#ratioCell <- array()
#for (c in 1:108){
#  ratioCell[c] <- mean(BTinform[[c]][, 4])
#}

#lines(ratioCell, type="l")
#########################above: not useful####################


######## New analysis: 2017-2-2 between-within variances 

# Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells 
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:8] 

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> reliability estimates of 20 samples from persons under condition cell 1 
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8]

# a function for calculating within/between var
BTvar <- function(mat){
  #each column is a group
  num_row <- dim(mat)[1]
  gmean <- mean(mat)
  cmean <- colMeans(mat)
  tvar <- sum((mat-gmean)^2)
  betvar <- sum(num_row*(cmean-gmean)^2)
  
  sum_a <- array() 
  for(i in 1:num_row){
    sum_a[i] <- sum((mat[i, ]-cmean)^2)
  }
  
  #witvar <- sum((mat-cmean)^2)
  witvar <- sum(sum_a)
  btRatio <- betvar/tvar
  
  return_list <- list()
  return_list$TotalVar <- tvar
  return_list$BetweenVar <- betvar
  return_list$WithinVar <- witvar
  return_list$BTratio <- btRatio
  return(return_list)
}


BT_tru_relia <- list()
BT_tradi_alpha <- list()
BT_tradi_l2 <- list()
BT_tradi_l4 <- list()
BT_new_alpha <- list()
BT_new_l2 <- list()
BT_new_l4 <- list()

for (c in 1:108){
  
  tru_relia <- matrix(0, 50, 20)
  tradi_alpha <- matrix(0, 50, 20)
  tradi_l2 <- matrix(0, 50, 20)
  tradi_l4 <- matrix(0, 50, 20)
  new_alpha <- matrix(0, 50, 20)
  new_l2 <- matrix(0, 50, 20)
  new_l4 <- matrix(0, 50, 20)

  for (i in 1:20){
   tru_relia[,i] <- restuls_conditions[[c]][[1]][[i]][, 2]
   tradi_alpha[,i] <- restuls_conditions[[c]][[1]][[i]][, 3]
   tradi_l2[,i] <- restuls_conditions[[c]][[1]][[i]][, 4]
   tradi_l4[,i] <- restuls_conditions[[c]][[1]][[i]][, 5]
   new_alpha[,i] <- restuls_conditions[[c]][[1]][[i]][, 6]
   new_l2[,i] <- restuls_conditions[[c]][[1]][[i]][, 7]
   new_l4[,i] <- restuls_conditions[[c]][[1]][[i]][, 8]
  }

  BT_tru_relia[[c]] <- BTvar(tru_relia)
  BT_tradi_alpha[[c]] <- BTvar(tradi_alpha)
  BT_tradi_l2[[c]] <- BTvar(tradi_l2)
  BT_tradi_l4[[c]] <- BTvar(tradi_l4)
  BT_new_alpha[[c]] <- BTvar(new_alpha)
  BT_new_l2[[c]] <- BTvar(new_l2)
  BT_new_l4[[c]] <- BTvar(new_l4)
}

BTRmatrix <- matrix(0, 108, 7)
for(i in 1:108){
  BTRmatrix[i, 1] <- unlist(BT_tru_relia[[i]][4])
  BTRmatrix[i, 2] <- unlist(BT_tradi_alpha[[i]][4])
  BTRmatrix[i, 3] <- unlist(BT_tradi_l2[[i]][4])
  BTRmatrix[i, 4] <- unlist(BT_tradi_l4[[i]][4])
  BTRmatrix[i, 5] <- unlist(BT_new_alpha[[i]][4])
  BTRmatrix[i, 6] <- unlist(BT_new_l2[[i]][4])
  BTRmatrix[i, 7] <- unlist(BT_new_l4[[i]][4])
}


plot(pop_re[, 2], xlim = c(0, 108), ylim = c(0,1), xlab = "108 cells", ylab = "True Reliability", col='red', pch=19)
lines(BTRmatrix[, 1], type = "b", col='green', pch=19)
points(BTRmatrix[, 2], pch=0, col='grey')
points(BTRmatrix[, 3], pch=1, col='grey')
points(BTRmatrix[, 4], pch=2, col='grey')
points(BTRmatrix[, 5], pch=15, col='purple')
points(BTRmatrix[, 6], pch=19, col='purple')
points(BTRmatrix[, 7], pch=17, col='purple')

#let's focus on carry over effects
carryindex <- df[, 5]!=0
plot(pop_re[, 2][carryindex], ylim = c(-5,1), ylab = "True Reliability", col='red', pch=19)
lines(pop_re[, 3][carryindex], pch=0, type = 'b')
abline(h=0)
abline(h=.5)
points(BTRmatrix[, 2][carryindex], pch=0, col='grey')

###### New Analysis: distribution of esimated reliability of a sample (out of 20 samples)
# The idea is that, since brain-washing is not possible, we have to estimate reliability based
# on the reponses (which might not be the same if we measure again). The estimated reliability 
# is expected to be close to the true reliability of that sample, but due to measurement errors, 
# there will be deviation. 

# Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells 
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:8] 

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> reliability estimates of 20 samples from persons under condition cell 1 
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8] --> standard deviation

allplots <- list()
for(cel in 1:108){
  plotarray <- list()
  y_min <- min(restuls_conditions[[cel]][[2]][, 2] - restuls_conditions[[cel]][[3]][, 2]/sqrt(50))
  y_max <- max(restuls_conditions[[cel]][[2]][, 2] + restuls_conditions[[cel]][[3]][, 2]/sqrt(50))
  plot(restuls_conditions[[cel]][[2]][, 2], xlab = "20 samples from the population", ylab = "True reliability +/- 1SE",
     ylim = c(y_min, y_max), 
     type = "p")
  arrows(c(1:20), restuls_conditions[[cel]][[2]][, 2] - restuls_conditions[[cel]][[3]][, 2]/sqrt(50), c(1:20), restuls_conditions[[cel]][[2]][, 2] + restuls_conditions[[cel]][[3]][, 2]/sqrt(50), 
       length = 0.05, angle = 90, code = 3)
  abline(h=mean(restuls_conditions[[cel]][[2]][, 2]), lty=2)
  plotarray[[1]] <- recordPlot()
  for(i in 3:8){
    if(i == 3){
      differentN <- "Estimated reliability traditional - alpha"
    }else if(i==4){
      differentN <- "Estimated reliability traditional - lambda2"
    }else if(i==5){
      differentN <- "Estimated reliability traditional - lambda4"
    }else if(i==6){
      differentN <- "Estimated reliability new - alpha"
    }else if(i==7){
      differentN <- "Estimated reliability new - lambda2"
    }else if(i==8){
      differentN <- "Estimated reliability new - lambda4"
    }
  
    y_min <- min(restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]/sqrt(50))
    y_max <- max(restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]/sqrt(50))
    plot(restuls_conditions[[cel]][[2]][, i], xlab = "20 samples from the population", ylab = differentN,
     ylim = c(y_min, y_max), 
     type = "p")
    arrows(c(1:20), restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]/sqrt(50), c(1:20), restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]/sqrt(50), 
       length = 0.05, angle = 90, code = 3)
    abline(h=mean(restuls_conditions[[cel]][[2]][, i]), lty=2)
    plotarray[[i-1]] <- recordPlot()
  }
  allplots[[cel]] <- plotarray
}

#display plots
cel_index = 108
for(i in 1:7){
  print(allplots[[cel_index]][[i]])
}


############## all 6 methods are ploted in one pic
   
allplots <- list()
for(cel in 1:108){  

  y_min <- min(restuls_conditions[[cel]][[2]][, 2] - restuls_conditions[[cel]][[3]][, 2]/sqrt(50))
  y_max <- max(restuls_conditions[[cel]][[2]][, 2] + restuls_conditions[[cel]][[3]][, 2]/sqrt(50))
  
  for(i in 6:8){
    y_min <- min(y_min, min(restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]/sqrt(50)))
    y_max <- max(y_max, max(restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]/sqrt(50)))
  }
  plot(restuls_conditions[[cel]][[2]][, 2], xlab = "20 samples from the population", ylab = "True reliability +/- 1SE",
     ylim = c(y_min,y_max), 
     type = "p")
  arrows(c(1:20), restuls_conditions[[cel]][[2]][, 2] - restuls_conditions[[cel]][[3]][, 2]/sqrt(50), c(1:20), restuls_conditions[[cel]][[2]][, 2] + restuls_conditions[[cel]][[3]][, 2]/sqrt(50), 
       length = 0.05, angle = 90, code = 3)
  abline(h=mean(restuls_conditions[[cel]][[2]][, 2]), lty=2)
  
  for(i in 6:8){
    if(i == 3){
      col <- "blue"
      pch <- 15
    }else if(i==4){
      col <- "blue"
      pch <- 18
    }else if(i==5){
      col <- "blue"
      pch <- 17
    }else if(i==6){
      col <- "Red"
      pch <- 15
    }else if(i==7){
      col <- "Red"
      pch <- 18
    }else if(i==8){
      col <- "Red"
      pch <- 17
    }
    points(restuls_conditions[[cel]][[2]][, i], col=col, pch=pch)
    arrows(c(1:20), restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]/sqrt(50), c(1:20), restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]/sqrt(50), 
           length = 0.05, angle = 90, code = 3)
    abline(h=mean(restuls_conditions[[cel]][[2]][, i]), lty=2, col=col)
  }
  allplots[[cel]] <- recordPlot()
}  
