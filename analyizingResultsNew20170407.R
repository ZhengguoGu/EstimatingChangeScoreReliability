#########################################################################################################
#
# This file contains the code for analyzing the results
#
# I analyize the results in the following steps:
#
############# PART I: Bias and precision  #################
# 1. Results for unidimensional data 
#    1) No carry-over effects
#    2) Carry-over effects
# 2. Results for multidimensional data
#    1) No carry-over effects
#    2) Carry-over effects
##########################################################
############# PART II: regression analysis ###############
# 3. Design factors that influence the difference between 
#    change-score reliability and its estimate.
##########################################################
#
# Zhengguo Gu, Tilburg University
# Last update: 07/04/2017
#########################################################################################################




############### PART I: Bias and precision #######################################################################

# 0.1 Load data - large sample size
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")

# 0.2 Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells (Note that there are two datasets, one for small sample one for large sample, and thus 108x2=216 cells)
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:8] 

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> reliability estimates of 20 samples from persons under condition cell 1 
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8]

# 0.3 check design factors
df$`correlated facets`==1     # unidimensional data
df$`correlated facets`!=1     # multidimensional data
df$`carry-over effects`==0    # No Carry-over effects
df$`carry-over effects`!=0    # Carry-over effects

###### draw bias and precision plots - 4 situations
situation <- 1 # "Unidimensional + No Carryover"
situation <- 2 # "Unidimensional + Carryover"
situation <- 3 # "Multidimensional + No Carryover"
situation <- 4 #"Multidimensional + Carryover"

if(situation == 1){
  index <- df$`correlated facets`==1 & df$`carry-over effects`==0
} else if (situation == 2){
  index <- df$`correlated facets`==1 & df$`carry-over effects`!=0
} else if (situation == 3){
  index <- df$`correlated facets`!=1 & df$`carry-over effects`==0
} else if (situation == 4){
  index <- df$`correlated facets`!=1 & df$`carry-over effects`!=0
}


for(j in 1:sum(index)){
  
  
  cel <- which(index == 1)[j]
  y_min <- 0
  y_max <- 1
  
  for(i in 3:8){
    
    y_min <- min(y_min, min(restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]))  #note that restuls_conditions[[cel]][[3]][, i]) contains the sd (corrected with n-1)
    y_max <- max(y_max, max(restuls_conditions[[cel]][[2]][, i] + restuls_conditions[[cel]][[3]][, i]))
    
  }
  filename <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2/Cell_", cel, ".png", sep = "")
  png(file=filename, width = 1200, height = 800, units = "px")
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
         c("true reliability", "traditional method + alpha", "traditional method + lambda2","traditional method + lambda4", 
           "item-score method + alpha", "item-score method + lambda2", "item-score method + lambda4", "averaged true reliability"),
         pch=c( 8, 0, 1, 2, 15, 19, 17, NA),
         lty = c(NA, NA, NA, NA, NA, NA, NA, 3),
         col=c( "red", "black", "black", "black", "blue", "blue", "blue", "red"),
         ncol=4, bty = "n")
  
  
  dev.off()
}


