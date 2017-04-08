#########################################################################################################
#
# This file contains the code for analyzing the results
#
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



######### 1. draw bias and precision plots - 4 situations: LEVEL 1 ##############

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

######### 2. Bias and precision - LEVEL 2 ###################

load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")


pop_re <- matrix(NA, nrow = 108, ncol = 8)
pop_sd <- matrix(NA, nrow = 108, ncol = 8)
for (i in 1:108){
  pop_re[i, ] <- colMeans(restuls_conditions[[i]][[2]])
  pop_sd[i, ] <- apply(restuls_conditions[[i]][[2]], 2, sd)
}


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


# plot pop_re
cellseq <- seq(1:108)

pop_reCOPY <- pop_re
pop_reCOPY[setdiff(cellseq, cellseq[index]), ] <- NA
pop_sdCOPY <- pop_sd
pop_sdCOPY[setdiff(cellseq, cellseq[index]), ] <- NA


layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
 
ymin <- min(0, pop_reCOPY[, 3] - pop_sdCOPY[, 3], pop_reCOPY[, 4] - pop_sdCOPY[, 4],
               pop_reCOPY[, 5] - pop_sdCOPY[, 5], pop_reCOPY[, 6] - pop_sdCOPY[, 6],
               pop_reCOPY[, 7] - pop_sdCOPY[, 7], pop_reCOPY[, 8] - pop_sdCOPY[, 8],
            na.rm = TRUE)
plot(pop_reCOPY[ , 2], xaxt = 'n', ylim = c(ymin,1), xlab = "Simulation Cells", ylab = "(Estimated) Change-Score Reliability +/- 1SD", col='red', pch=8)
axis(1, at=cellseq[index], las=2, cex.axis = .7)
points(pop_reCOPY[, 3], pch=0)
points(pop_reCOPY[, 4], pch=1)
points(pop_reCOPY[, 5], pch=2)
points(pop_reCOPY[, 6], pch=15, col='blue')
points(pop_reCOPY[, 7], pch=19, col='blue')
points(pop_reCOPY[, 8], pch=17, col='blue')
arrows(c(1:108), pop_reCOPY[, 3] - pop_sdCOPY[, 3], c(1:108), pop_reCOPY[, 3] + pop_sdCOPY[, 3], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_reCOPY[, 4] - pop_sdCOPY[, 4], c(1:108), pop_reCOPY[, 4] + pop_sdCOPY[, 4], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_reCOPY[, 5] - pop_sdCOPY[, 5], c(1:108), pop_reCOPY[, 5] + pop_sdCOPY[, 5], length = 0.05, angle = 90, code = 3, col = 'black')
arrows(c(1:108), pop_reCOPY[, 6] - pop_sdCOPY[, 6], c(1:108), pop_reCOPY[, 6] + pop_sdCOPY[, 6], length = 0.05, angle = 90, code = 3, col = 'blue')
arrows(c(1:108), pop_reCOPY[, 7] - pop_sdCOPY[, 7], c(1:108), pop_reCOPY[, 7] + pop_sdCOPY[, 7], length = 0.05, angle = 90, code = 3, col = 'blue')
arrows(c(1:108), pop_reCOPY[, 8] - pop_sdCOPY[, 8], c(1:108), pop_reCOPY[, 8] + pop_sdCOPY[, 8], length = 0.05, angle = 90, code = 3, col = 'blue')

for(i in 1:108){
  if(i %in% cellseq[index]){
     abline(v=i, lty="dotted", lwd=.1)
  }
}
abline(v=c(37,73), col="red")
mtext("Cell no.1 ~ 36: Short test (9 items)", side = 3, line=0, at=20)
mtext("Cell no.37 ~ 72: Medium-length test (21 items)", side = 3, line=0, at=55)
mtext("Cell no.73 ~ 108: Long test (36 items)", side = 3, line=0, at=90)
#axis(1, at=c(37,73), labels=c("37","73"))

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
       c("true reliability", "traditional method + alpha", "traditional method + lambda2","traditiona method + lambda4", 
         "item-score method + alpha", "item-score method + lambda 2", "item-score method + lambda4"),
       pch=c(8, 0, 1,2, 15, 19, 17),
       col=c("red", "black", "black", "black", "blue", "blue", "blue"),
       ncol=4, bty = "n")
dev.off()

######### 3. distance measure:  ######################################

load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")

pop_re <- matrix(NA, nrow = 108, ncol = 8)
pop_sd <- matrix(NA, nrow = 108, ncol = 8)
for (i in 1:108){
  pop_re[i, ] <- colMeans(restuls_conditions[[i]][[2]])
  pop_sd[i, ] <- apply(restuls_conditions[[i]][[2]], 2, sd)
}

population_rel1SD_lower <- pop_re[, 2] - pop_sd[, 2]  # -1SD
population_rel1SD_upper <- pop_re[, 2] + pop_sd[, 2]  # +1SD

population_rel3SD_lower <- pop_re[, 2] - 3*pop_sd[, 2]  # -3SD
population_rel3SD_upper <- pop_re[, 2] + 3*pop_sd[, 2]  # +3SD

SDrange1SD <- cbind(population_rel1SD_lower, population_rel1SD_upper)
SDrange3SD <- cbind(population_rel3SD_lower, population_rel3SD_upper)


for(cel in 1:108){

  SD1sum_in_regionTraAlpha[cel] <- 0 # # of estimated reliablity within the region of +/- 1SD: traditional + alpha
  SD1sum_in_regionTraL2[cel] <- 0    # # of estimated reliablity within the region of +/- 1SD: traditional + lambda 2
  SD1sum_in_regionTraL4[cel] <- 0    # # of estimated reliablity within the region of +/- 1SD: traditional + lambda 4

  SD1sum_in_regionItemAlpha[cel] <- 0 # # of estimated reliablity within the region of +/- 1SD: Item method + alpha
  SD1sum_in_regionItemL2[cel] <- 0    # # of estimated reliablity within the region of +/- 1SD: Item method + lambda 2
  SD1sum_in_regionItemL4[cel] <- 0    # # of estimated reliablity within the region of +/- 1SD: Item method + lambda 4

  SD3sum_in_regionTraAlpha[cel] <- 0 # # of estimated reliablity within the region of +/- 3SD: traditional + alpha
  SD3sum_in_regionTraL2[cel] <- 0    # # of estimated reliablity within the region of +/- 3SD: traditional + lambda 2
  SD3sum_in_regionTraL4[cel] <- 0    # # of estimated reliablity within the region of +/- 3SD: traditional + lambda 4

  SD3sum_in_regionItemAlpha[cel] <- 0 # # of estimated reliablity within the region of +/- 3SD: Item method + alpha
  SD3sum_in_regionItemL2[cel] <- 0    # # of estimated reliablity within the region of +/- 3SD: Item method + lambda 2
  SD3sum_in_regionItemL4[cel] <- 0    # # of estimated reliablity within the region of +/- 3SD: Item method + lambda 4

  for(l in 1:20){
  
    SD1sum_in_regionTraAlpha[cel] <- SD1sum_in_regionTraAlpha[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,3] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,3] >= SDrange1SD[cel, 1])
    SD1sum_in_regionTraL2[cel] <- SD1sum_in_regionTraL2[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,4] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,4] >= SDrange1SD[cel, 1])
    SD1sum_in_regionTraL4[cel] <- SD1sum_in_regionTraL4[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,5] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,5] >= SDrange1SD[cel, 1])
  
    SD3sum_in_regionTraAlpha[cel] <- SD3sum_in_regionTraAlpha[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,3] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,3] >= SDrange3SD[cel, 1])
    SD3sum_in_regionTraL2[cel] <- SD3sum_in_regionTraL2[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,4] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,4] >= SDrange3SD[cel, 1])
    SD3sum_in_regionTraL4[cel] <- SD3sum_in_regionTraL4[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,5] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,5] >= SDrange3SD[cel, 1])
  
    SD1sum_in_regionItemAlpha[cel] <- SD1sum_in_regionItemAlpha[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,6] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,6] >= SDrange1SD[cel, 1])
    SD1sum_in_regionItemL2[cel] <- SD1sum_in_regionItemL2[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,7] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,7] >= SDrange1SD[cel, 1])
    SD1sum_in_regionItemL4[cel] <- SD1sum_in_regionItemL4[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,8] <= SDrange1SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,8] >= SDrange1SD[cel, 1])
  
    SD3sum_in_regionItemAlpha[cel] <- SD3sum_in_regionItemAlpha[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,6] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,6] >= SDrange3SD[cel, 1])
    SD3sum_in_regionItemL2[cel] <- SD3sum_in_regionItemL2[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,7] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,7] >= SDrange3SD[cel, 1])
    SD3sum_in_regionItemL4[cel] <- SD3sum_in_regionItemL4[cel] + sum(restuls_conditions[[cel]][[1]][[l]][,8] <= SDrange3SD[cel, 2] & restuls_conditions[[cel]][[1]][[l]][,8] >= SDrange3SD[cel, 1])
  }
}  

dist_measureSD <- cbind(SD1sum_in_regionTraAlpha, SD3sum_in_regionTraAlpha, SD1sum_in_regionTraL2, SD3sum_in_regionTraL2, 
                        SD1sum_in_regionTraL4, SD3sum_in_regionTraL4, SD1sum_in_regionItemAlpha, SD3sum_in_regionItemAlpha,
                        SD1sum_in_regionItemL2, SD3sum_in_regionItemL2, SD1sum_in_regionItemL4, SD3sum_in_regionItemL4)/1000
dist_measureSD <- cbind(pop_sd[, 2], dist_measureSD)  # this is the matrix containing the 
colnames(dist_measureSD)[1] <- "1SD"

# 
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

dist_measureSD[index, ]