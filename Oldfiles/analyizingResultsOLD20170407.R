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

######### 3. distance measure:  part 1 ######################################

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

SD1sum_in_regionTraAlpha <- array() # # of estimated reliablity within the region of +/- 1SD: traditional + alpha
SD1sum_in_regionTraL2 <- array()    # # of estimated reliablity within the region of +/- 1SD: traditional + lambda 2
SD1sum_in_regionTraL4 <- array()    # # of estimated reliablity within the region of +/- 1SD: traditional + lambda 4

SD1sum_in_regionItemAlpha <- array() # # of estimated reliablity within the region of +/- 1SD: Item method + alpha
SD1sum_in_regionItemL2 <- array()    # # of estimated reliablity within the region of +/- 1SD: Item method + lambda 2
SD1sum_in_regionItemL4 <- array()    # # of estimated reliablity within the region of +/- 1SD: Item method + lambda 4

SD3sum_in_regionTraAlpha <- array() # # of estimated reliablity within the region of +/- 3SD: traditional + alpha
SD3sum_in_regionTraL2 <- array()    # # of estimated reliablity within the region of +/- 3SD: traditional + lambda 2
SD3sum_in_regionTraL4 <- array()    # # of estimated reliablity within the region of +/- 3SD: traditional + lambda 4

SD3sum_in_regionItemAlpha <- array() # # of estimated reliablity within the region of +/- 3SD: Item method + alpha
SD3sum_in_regionItemL2 <- array()    # # of estimated reliablity within the region of +/- 3SD: Item method + lambda 2
SD3sum_in_regionItemL4 <- array()    # # of estimated reliablity within the region of +/- 3SD: Item method + lambda 4

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
dist_measureSD <- cbind(dist_measureSD, df)
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

resultdistmeasureSD <- dist_measureSD[index, ]

############# 4.  given the 108 cells, which cells generate negative reliabilites? ####################
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")

weird_index <- matrix(NA, 108, 6)
for(i in 1:108){
  allestimates <- matrix(unlist(restuls_conditions[[i]][2]), nrow = 20, ncol = 8)[, c(-1,-2)] 
  weird_index[i, ] = (colSums(allestimates < 0)>=1) #negative reliability
}

df[, 6:11] <- weird_index 
names(df)[6:11] <- c("neg traditional alpha", "neg traditional lambda2", "neg traditional lambda4", 
                     "neg item alpha", "neg item lambda2", "neg item lambda4")
table(df$`carry-over effects`,df$`neg traditional alpha`)
table(df$`carry-over effects`,df$`neg traditional lambda2`)
table(df$`carry-over effects`,df$`neg traditional lambda4`)

weird_index2 <- rep(0, 108)
for(i in 1:108){
  allestimates <- matrix(unlist(restuls_conditions[[i]][2]), nrow = 20, ncol = 8) 
  if(sum(allestimates[, 2:8] >1)>=1){ #reliability >1
    weird_index2[i] <- 1
  }
}
sum(weird_index2) # non of them have reliability > 1.

############# 4.Extra: Why negative estimated reliabiity? A toy example  ################
# Here we simulate toy example to see why negative estimated reliability happens. 
library(psychometric)
num_items <- 9                # short
parallel_item <- 0               # Non parallel
dimension <- 1           # Unidimensional
sd_change <- sqrt(.14)        # Small variance in change
strongweak_carry <- c(0, 10, 1)  # 0: no carryover, 1: weak, 2: strong
num_persons <- 1000              # 1000 persons

set.seed(110)
# item parameter
itempar <- matrix(NA,num_items,5)
itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
avg_beta <- runif(num_items, 0, 1.25)
itempar[,2] <- avg_beta - 1
itempar[,3] <- avg_beta - .5
itempar[,4] <- avg_beta + .5
itempar[,5] <- avg_beta + 1


theta <- Unichange_sim(num_persons, sd_change) # make use of the Unichange_sim() funtion in ZhengguoFuntions.R
theta_pre <- theta[[1]]
theta_post <- theta[[2]]

id <- vector()
for(d in 1: dimension){
  id <- cbind(id, rep(d, num_items/dimension))
}
id <- as.vector(id) 


responses <- GRM_sim(theta_pre, itempar, id)  # make use of the GRM_sim() funtion in ZhengguoFuntions.R
response_pre <- responses[[1]]
true_pre <- responses[[2]]
sum_pre <- rowSums(response_pre) 
sum_true_pre <- rowSums(true_pre)
E1 <- response_pre - true_pre
sum_E1 <- rowSums(E1) 
mean(sum_E1)
responses <- GRM_sim(theta_post, itempar, id)
response_post <- responses[[1]]
true_post <- responses[[2]]
sum_true_post <- rowSums(true_post)

truechange_sumscores <- sum_true_post - sum_true_pre

CarryoverEffect <- carry_over(response_pre, response_post) # if everyone showes carry-over

true_rel_strong <- array()
true_rel_weak <- array()
r_post_strong <- array()
r_post_weak <- array()
rel_alpha_strong <- array()
rel_alpha_weak <- array()
true_relVar_strong <- array()
true_relVar_weak <- array()
var_post_strong <- array()
corr_strong <- array()
var_post_weak <- array()
corr_weak <- array()
numeritor_strong <- array()
demominator_strong <- array()
numeritor_weak <- array()
demominator_weak <- array()
cor_E1E2strong <- array()
cor_E1E2weak <- array()
j <- 1
for(i in c(0, 100,  200,  300,  400,  500,  600,  700,  800,  900, 1000)){

  subject_index <- sample(1:1000, size = i, replace = FALSE) #number of persons showed carry-over effect
  
  response_Strong <- response_post 
  response_Strong[subject_index,] <- CarryoverEffect[[1]][subject_index,]  #strong
  E2_strong <- response_Strong - true_post
  sum_E2strong <- rowSums(E2_strong)
  
  response_Weak <- response_post
  response_Weak[subject_index,] <- CarryoverEffect[[2]][subject_index,]  #weak
  E2_weak <- response_Weak - true_post
  sum_E2weak <- rowSums(E2_weak)
  
  cor_E1E2strong[j] <- cor(sum_E1, sum_E2strong)
  cor_E1E2weak[j] <- cor(sum_E1, sum_E2weak)
  
  sum_post_strong <- rowSums(response_Strong)
  sum_post_weak <- rowSums(response_Weak)
  change_sumscores_strong <- sum_post_strong - sum_pre
  change_sumscores_weak <- sum_post_weak - sum_pre
  true_rel_strong[j] <- (cor(truechange_sumscores, change_sumscores_strong))^2
  true_rel_weak[j] <- (cor(truechange_sumscores, change_sumscores_weak))^2
  r_pre <- psychometric::alpha(response_pre)  # ! cronback alpha is used here. 
  r_post_strong[j] <- psychometric::alpha(response_Strong)
  r_post_weak[j] <- psychometric::alpha(response_Weak)
  
  var_post_strong[j] <- var(sum_post_strong)
  corr_strong[j] <- cor(sum_pre, sum_post_strong)
  numeritor_strong[j] <- (var(sum_pre) * r_pre + var(sum_post_strong) * r_post_strong[j] - 2 * cor(sum_pre, sum_post_strong) * sd(sum_pre) * sd(sum_post_strong))
  demominator_strong[j] <- (var(sum_pre) + var(sum_post_strong) - 2 * cor(sum_pre, sum_post_strong) * sd(sum_pre) * sd(sum_post_strong))
  rel_alpha_strong[j] <- numeritor_strong[j]/demominator_strong[j]
  
  var_post_weak[j] <- var(sum_post_weak)
  corr_weak[j] <- cor(sum_pre, sum_post_weak)
  numeritor_weak[j] <- (var(sum_pre) * r_pre + var(sum_post_weak) * r_post_weak[j] - 2 * cor(sum_pre, sum_post_weak) * sd(sum_pre) * sd(sum_post_weak))
  demominator_weak[j] <- (var(sum_pre) + var(sum_post_weak) - 2 * cor(sum_pre, sum_post_weak) * sd(sum_pre) * sd(sum_post_weak))
  rel_alpha_weak[j] <- numeritor_weak[j]/demominator_weak[j]
  
  j <- j+1
}


layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
ymin <- min(rel_alpha_strong)
plot(true_rel_strong, ylim = c(ymin, 1), type = "b", col = "red", xlab="Proportion of persons to whom strong/weak carry-over effects happen",
     ylab = "Reliability", xaxt="n", pch=8)
axis(1, at=1:11,labels=c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))
lines(rel_alpha_strong, type = "b", col = "blue", pch=15)
lines(rel_alpha_weak, type = "b", col = "green", pch=16)
abline(h=0, col = "black", lty=3)

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
       c("true change-score reliability", "estimated reliability (strong carry-over effects)", "estimated reliability (weak carry-over effects)"),
       pch=c(8, 15, 16),
       col=c("red", "blue", "green"),
       ncol=1, bty = "n")
dev.off()

layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
par(mar=c(5,4,4,5)+.1)
plot(numeritor_strong, type = "b", ylim = c(min(numeritor_strong, demominator_strong), max(numeritor_strong, demominator_strong)), col="blue", pch = 15, 
     xlab = "Proportion of persons to whom strong carry-over effects happen", ylab = "Numerical values of the numerator and denominator of equation (2)", xaxt="n")
axis(1, at=1:11,labels=c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))
lines(demominator_strong, type = "b", col="blue", pch = 16)
abline(h=0, col='blue', lty=3)
par(new=TRUE)
plot(rel_alpha_strong,type="b",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
abline(h=0, col="red", lty=3)
axis(4)
mtext("Estimated reliability",side=4,line=3)

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
       c("numerator", "denominator", "estimated reliability"),
       pch=c(15, 16, 1),
       col=c("blue", "blue", "red"),
       ncol=3, bty = "n")
dev.off()


############### 6. Regression analysis  ##############################
### Note that throughout the code below, many data matrices are called somethingANOVA, for example, "reANOVA" or "distance ANOVA". This is because
### I wanted to do ANOVA at first and later on decided to do regression analysis. 
### Please do not be confued by the names. 



load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
set.seed(110) #important!
index_sample <- sample(1:20, 108, replace=TRUE) # replace is True here! important

reANOVA <- matrix(NA, 1, 13)
for(cel in 1:108){
  
  reSample <- restuls_conditions[[cel]][[1]][[index_sample[cel]]] # one of the 20 samples are saved to reSample
  reSample <- cbind(reSample, t(matrix(rep(df[cel,], 50), 5, 50)))
  reANOVA <- rbind(reANOVA, reSample)
  
}
reANOVA <- reANOVA[-1,]
save(reANOVA, file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/reANOVAlarge.RData")

###
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")
set.seed(110) #important!
index_sample <- sample(1:20, 108, replace=TRUE) 

reANOVAs <- matrix(NA, 1, 13)
for(cel in 1:108){
  
  reSample <- restuls_conditions[[cel]][[1]][[index_sample[cel]]] 
  reSample <- cbind(reSample, t(matrix(rep(df[cel,], 50), 5, 50)))
  reANOVAs <- rbind(reANOVAs, reSample)
  
}
reANOVAs <- reANOVAs[-1,]
save(reANOVAs, file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/reANOVAsmall.RData")

#now the analysis
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/reANOVAlarge.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/reANOVAsmall.RData")

reANOVA <- as.data.frame(reANOVA)
names(reANOVA) <- c("useless", "true", "trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                    "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects")
reANOVA$size <- "large"
reANOVAs <- as.data.frame(reANOVAs)
names(reANOVAs) <- c("useless", "true", "trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                     "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects")
reANOVAs$size <- "small"

combANOVA <- rbind(reANOVA,reANOVAs) # this is the raw data.. 10800 rows
#calculate the distance for each row, this is estimated - true
distanceANOVA <- matrix(NA, 10800, 12)
distanceANOVA[, 1] <-  as.numeric(combANOVA$`trad alpha`) - as.numeric(combANOVA$true)
distanceANOVA[, 2] <-  as.numeric(combANOVA$`trad lambda2`) - as.numeric(combANOVA$true)
distanceANOVA[, 3] <-  as.numeric(combANOVA$`trad lambda4`) - as.numeric(combANOVA$true)
distanceANOVA[, 4] <-  as.numeric(combANOVA$`item alpha`) - as.numeric(combANOVA$true)
distanceANOVA[, 5] <-  as.numeric(combANOVA$`item lambda2`) - as.numeric(combANOVA$true)
distanceANOVA[, 6] <-  as.numeric(combANOVA$`item lambda4`) - as.numeric(combANOVA$true)
# note that distanceANOVA[, 7:12] <- combANOVA[, 9:14] is not recommended! The resulting file will 31 GB. (although I dont know the reason why)
distanceANOVA <- as.data.frame(distanceANOVA)
distanceANOVA$V7 <- factor(as.numeric(combANOVA$`test length`))
distanceANOVA$V8 <- factor(as.numeric(combANOVA$`parallel item`))
distanceANOVA$V9 <- factor(as.numeric(combANOVA$`correlated facets`))
distanceANOVA$V10 <- factor(as.numeric(combANOVA$`maginute of sd`))
distanceANOVA$V11 <- factor(as.numeric(combANOVA$`carry-over effects`))
distanceANOVA$V12 <- factor(combANOVA$size)
names(distanceANOVA) <- c("trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                      "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects", "size")
save(distanceANOVA, file="D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/distanceANOVAdata.RData")

########## 6.1 Regression traditional method with alpha? ##########
fit <- lm(distanceANOVA$'trad alpha' ~ distanceANOVA$'test length' + 
            distanceANOVA$'parallel item' + 
            distanceANOVA$'correlated facets' +
            distanceANOVA$'maginute of sd' + 
            distanceANOVA$'carry-over effects' + 
            distanceANOVA$'size')
anova(fit)
fit.summary <- summary(fit)

fit_tl2 <- lm(distanceANOVA$'trad lambda2' ~ distanceANOVA$'test length' + 
                distanceANOVA$'parallel item' + 
                distanceANOVA$'correlated facets' +
                distanceANOVA$'maginute of sd' + 
                distanceANOVA$'carry-over effects' + 
                distanceANOVA$'size' )
anova(fit_tl2)
fit.summary <- summary(fit_tl2)

fit_tl4 <- lm(distanceANOVA$'trad lambda4' ~ distanceANOVA$'test length' + 
                distanceANOVA$'parallel item' + 
                distanceANOVA$'correlated facets' +
                distanceANOVA$'maginute of sd' + 
                distanceANOVA$'carry-over effects' + 
                distanceANOVA$'size' )
anova(fit_tl4)
fit.summary <- summary(fit_tl4)

# Which factors influence the item-score method?
fit2 <- lm(distanceANOVA$'item alpha' ~ distanceANOVA$'test length' + 
             distanceANOVA$'parallel item' + 
             distanceANOVA$'correlated facets' +
             distanceANOVA$'maginute of sd' + 
             distanceANOVA$'carry-over effects' + 
             distanceANOVA$'size' )
anova(fit2)
fit.summary <- summary(fit2)

fit3 <- lm(distanceANOVA$'item lambda2' ~ distanceANOVA$'test length' + 
             distanceANOVA$'parallel item' + 
             distanceANOVA$'correlated facets' +
             distanceANOVA$'maginute of sd' + 
             distanceANOVA$'carry-over effects' + 
             distanceANOVA$'size' )
anova(fit3)
fit.summary <- summary(fit3)


############# 7. the toy example (table 1)  #############
library(psychometric)
num_items <- 3                
dimension <- 1           # Unidimensional
sd_change <- sqrt(.14)        # Small variance in change
num_persons <- 15              # 1000 persons

set.seed(120)
# item parameter
itempar <- matrix(NA,num_items,6)
itempar[,1] <- runif(1,1.5,2.5)   # discrimination
avg_beta <- runif(1, 0, 1.25)
itempar[,2] <- avg_beta - 1
itempar[,3] <- avg_beta - .5
itempar[,4] <- avg_beta
itempar[,5] <- avg_beta + .5
itempar[,6] <- avg_beta + 1


theta <- Unichange_sim(num_persons, sd_change) # make use of the Unichange_sim() funtion in ZhengguoFuntions.R
theta_pre <- theta[[1]]
theta_post <- theta[[2]]

id <- vector()
for(d in 1: dimension){
  id <- cbind(id, rep(d, num_items/dimension))
}
id <- as.vector(id) 


responses <- GRM_sim(theta_pre, itempar, id)  # make use of the GRM_sim() funtion in ZhengguoFuntions.R
response_pre <- responses[[1]]
sum_pre <- rowSums(response_pre) 


responses <- GRM_sim(theta_post, itempar, id)
response_post <- responses[[1]]
sum_post <- rowSums(response_post)

# item method
item_change <- response_post - response_pre
RELitem_method <- psychometric::alpha(item_change)

# traditional method
r_pre <- psychometric::alpha(response_pre)
r_post <- psychometric::alpha(response_post)
RELtraditional <- (var(sum_pre) * r_pre + var(sum_post) * r_post - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))/(var(sum_pre) + var(sum_post) - 2 * cor(sum_pre, sum_post) * sd(sum_pre) * sd(sum_post))

toydata <- cbind(response_pre, response_post, item_change)
toydata <- data.frame(toydata)
colnames(toydata) <- c("item 1", "item 2", "item 3", "item 1", "item 2", "item 3", "item 1", "item 2", "item 3")
write.table(toydata, file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170411 newversion2 NewSupplementary/toyexample.txt", sep = '\t')
