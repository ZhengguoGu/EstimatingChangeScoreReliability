###################################################################################
#
#  A few plots
#
###################################################################################

##########################################################
# cross all samples of persons
results <- r_avg

plot(results[,1], results[, 2], xlab = "true reliability - proportion variance",
     ylab = "true reliability - correlation", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 3], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by alpha )", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 4], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by lambda2)", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 6], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - alpha", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results[, 2], results[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")


#######################
# one sample 

p <- sample(1:maxp, 1)

results_oneS <- sample_results[[p]]

plot(results_oneS[,1], results_oneS[, 2], xlab = "true reliability - proportion variance",
     ylab = "true reliability - correlation", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 3], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by alpha )", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 4], xlab = "true reliability - correlation",
     ylab = "estimated reliability (pre and post reliability estimated by lambda2)", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 6], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - alpha", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")

###

plot(results_oneS[, 2], results_oneS[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0, 1), ylim=c(0, 1),asp=1)
lines(c(0,1),c(0,1),col="red")


######## trial #################

plot(results[, 2], results[, 7], xlab = "true reliability - correlation",
     ylab = "estimated reliability (item-level) - lambda2", 
     xlim=c(0.3, 0.5), ylim=c(0.3, 0.5),asp=1)
lines(c(0,1),c(0,1),col="red")

for (i in 1:maxp){
  
  r_02 <- sample_results[[i]][, 2]
  r_22 <- sample_results[[i]][, 7]
  
  points(r_02, r_22, pch=8, col = 'grey')
}

## some pretty colors
library(RColorBrewer)
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

for (i in 1:maxp){
  
  r_02 <- sample_results[[i]][, 2]
  r_22 <- sample_results[[i]][, 7]
  
  ## compute 2D kernel density, see MASS book, pp. 130-131
  z <- kde2d(r_02, r_22, n=100)
  
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
}

