##########################################################################################
#
# This file contains the code for analyzing the results
#
# Note that because I re-organize the paper and discuss bia and precision first at 
# Level 1 and then at Level 2, I rewite the code for reproducing the results in the paper.
#
# Zhengguo Gu, Tilburg University
# Last update: 30/03/2017
##########################################################################################

# 1. load simulation results
#load('results20170122.RData')
load('results20170222smallsample.RData')

# 2. Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells (Note that there are two datasets, one for small sample one for large sample, and thus 108x2=216 cells)
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:8] 

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> reliability estimates of 20 samples from persons under condition cell 1 
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8]

#####################################################################
# 3. Bias and precision at Level 1 (sample level)
#####################################################################

############## all 6 methods are ploted in one pic


allplots <- list()

for(cel in 1:108){
  
  y_min <- 1
  y_max <- 1
  
  for(i in 3:8){
    
    y_min <- min(y_min, min(restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]))  #note that restuls_conditions[[cel]][[3]][, i]) contains the sd (corrected with n-1)
    y_max <- max(y_max, max(restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]))
    
  }
  layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
  xtitle <- sprintf("20 samples from the population: Cell %d",cel)
  plot(restuls_conditions[[cel]][[2]][, 2], xlab = xtitle, ylab = "True and estimated change-score reliability +/- 1SD",
       ylim = c(y_min,y_max), 
       type = "p",
       pch = 8,
       col = "red",
       cex = 1.5)
  arrows(c(1:20), restuls_conditions[[cel]][[2]][, 2] - restuls_conditions[[cel]][[3]][, 2], c(1:20), restuls_conditions[[cel]][[2]][, 2] + restuls_conditions[[cel]][[3]][, 2], 
         length = 0.05, angle = 90, code = 3, col = 'red')
  abline(h=mean(restuls_conditions[[cel]][[2]][, 2]), lty=2, col="red")
  
  for(i in 3:8){
    if(i == 3){
      col <- "black"
      pch <- 0
    }else if(i==4){
      col <- "black"
      pch <- 1
    }else if(i==5){
      col <- "black"
      pch <- 2
    }else if(i==6){
      col <- "blue"
      pch <- 15
    }else if(i==7){
      col <- "blue"
      pch <- 19
    }else if(i==8){
      col <- "blue"
      pch <- 17
    }
    points(restuls_conditions[[cel]][[2]][, i], col=col, pch=pch)
    arrows(c(1:20), restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i], c(1:20), restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i], 
           length = 0.05, angle = 90, code = 3, col=col)
    #abline(h=mean(restuls_conditions[[cel]][[2]][, i]), lty=2, col=col)
    #segments(20.2, mean(restuls_conditions[[cel]][[2]][, i]), 21, mean(restuls_conditions[[cel]][[2]][, i]), lty=2, col=col, cex=1.5, pch=pch)
    #segments(0, mean(restuls_conditions[[cel]][[2]][, i]), .8 , mean(restuls_conditions[[cel]][[2]][, i]), lty=2, col=col, cex=1.5, pch=pch)
  }
  
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center", "groups",
         c("true reliability", "traditional method + alpha", "traditional method + lambda2","traditiona method + lambda4", 
           "item-score method + alpha", "item-score method + lambda2", "item-score method + lambda4", "averaged reliability"),
         pch=c( 8, 0, 1, 2, 15, 19, 17, NA),
         lty = c(NA, NA, NA, NA, NA, NA, NA, 3),
         col=c( "red", "black", "black", "black", "blue", "blue", "blue", "red"),
         ncol=4, bty = "n")
  
  allplots[[cel]] <- recordPlot() #note that at this moment, I couldnt figure out how to print out all the plots automatically, with 6 plots on the same page. (save pics manually)
  
  dev.off()
}

allplots[[108]] #manually save all plots from 1 to 108.


#################################################################
#### 4. Bias and precision at Level 2
#################################################################

# analysis at the population level: that is, we average reliability estimates across the 20 samples of persons.
# The idea is that the estimates reflect the reliability at the population level. 

pop_re <- matrix(NA, nrow = 108, ncol = 8)
pop_sd <- matrix(NA, nrow = 108, ncol = 8)
for (i in 1:108){
  pop_re[i, ] <- colMeans(restuls_conditions[[i]][[2]])
  pop_sd[i, ] <- apply(restuls_conditions[[i]][[2]], 2, sd)
}

# plot pop_re
layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
plot(pop_re[, 2], xlim = c(0, 108), ylim = c(0,1), xlab = "108 cells", ylab = "(Estimated) Change-Score Reliability +/- 1SD", col='red', pch=8)
points(pop_re[, 3], pch=0)
points(pop_re[, 4], pch=1)
points(pop_re[, 5], pch=2)
points(pop_re[, 6], pch=15, col='blue')
points(pop_re[, 7], pch=19, col='blue')
points(pop_re[, 8], pch=17, col='blue')
arrows(c(1:108), pop_re[, 3] - pop_sd[, 3], c(1:108), pop_re[, 3] + pop_sd[, 3], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_re[, 4] - pop_sd[, 4], c(1:108), pop_re[, 4] + pop_sd[, 4], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_re[, 5] - pop_sd[, 5], c(1:108), pop_re[, 5] + pop_sd[, 5], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_re[, 6] - pop_sd[, 6], c(1:108), pop_re[, 6] + pop_sd[, 6], length = 0.05, angle = 90, code = 3, col = 'blue')
arrows(c(1:108), pop_re[, 7] - pop_sd[, 7], c(1:108), pop_re[, 7] + pop_sd[, 7], length = 0.05, angle = 90, code = 3, col = 'blue')
arrows(c(1:108), pop_re[, 8] - pop_sd[, 8], c(1:108), pop_re[, 8] + pop_sd[, 8], length = 0.05, angle = 90, code = 3, col = 'blue')

for(i in 1:108){
  abline(v=i, lty="dotted")
}
abline(v=c(37,73), col="red")
mtext("Cell no.1 ~ 36: Short test (9 items)", side = 3, line=0, at=20)
mtext("Cell no.37 ~ 72: Medium-length test (21 items)", side = 3, line=0, at=55)
mtext("Cell no.73 ~ 108: Long test (36 items)", side = 3, line=0, at=90)
axis(1, at=c(37,73), labels=c("37","73"))

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
       c("true reliability", "traditional method + alpha", "traditional method + lambda2","traditiona method + lambda4", 
         "item-score method + alpha", "item-score method + lambda 2", "item-score method + lambda4"),
       pch=c(8, 0, 1,2, 15, 19, 17),
       col=c("red", "black", "black", "black", "blue", "blue", "blue"),
       ncol=4, bty = "n")


#################################################################
#####  ANOVA: compare the six methods. 
#####  Let the item-level method with alpha be the reference (Does it make sense?)
#################################################################
g1 <- pop_re[, 6]  # item-level, alpha
g2 <- pop_re[, 7]  # item-level, lambda2
g3 <- pop_re[, 8]  # item-level, lambda4
g4 <- pop_re[, 3]  # traditional, alpha
g5 <- pop_re[, 4]  # traditional, lambda2
g6 <- pop_re[, 5]  # traditional, lambda4 

g123456 <- c(g1, g2, g3, g4, g5, g6)
groups <- factor(rep(c('g1', 'g2', 'g3', 'g4', 'g5', 'g6'), each = 108))

fit <- lm(formula = g123456 ~ groups)
anova(fit)
