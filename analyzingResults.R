##########################################################################################
#
# This file contains the code for analyzing the results
#
##########################################################################################

# 1. load simulation results
load('results20161223.RData')

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

#-- generate a data matrix for MANOVA
datamanova <- matrix(NA, nrow = 1, ncol=(8+5))
results <- data.matrix(df)
for (i in 1:108){
  conditions <- matrix(NA, nrow = 20, ncol = 5)
  for (j in 1:20){
    conditions[j, ] <- t(results[i, ])
  }
  datamanova <- rbind(datamanova, cbind(restuls_conditions[[i]][[2]], conditions))
}
datamanova <- datamanova[-1, ] #the first line is NA

# remove the first column, as it is not in our analysis. 
datamanova <- datamanova[, -1]

# calculate the distance between each estimate and the true reliability
str(datamanova)
reliaManova <- matrix(NA, nrow=2160, ncol=11)
for (k in 1:6){
  reliaManova[, k] <- datamanova[, k+1] - datamanova[, 1]
}
reliaManova[, 7:11] <- datamanova[, 8:12]
reliaManova <- data.frame(reliaManova)
# columns 7:11 --> factors
for (i in 7:11){
  reliaManova[, i] <- factor(reliaManova[, i])
}
reliaManova[, 7]
reliaManova[, 8]
reliaManova[, 9]
reliaManova[, 10]
reliaManova[, 11]

# MANOVA, no interactions
fit <- manova(as.matrix(reliaManova[, 1:6]) ~ reliaManova[, 7] + reliaManova[, 8] + reliaManova[, 9] + reliaManova[, 10] + reliaManova[, 11])
summary(fit)
