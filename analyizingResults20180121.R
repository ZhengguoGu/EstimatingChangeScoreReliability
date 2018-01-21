#########################################################################################################
#
# This file contains the code for analyzing the results
#
#
# Zhengguo Gu, Tilburg University
# Last update: 21/January/2018
#########################################################################################################




############### 0. Check data #######################################################################
# Note: There are in total 

# 0.1 Load data 
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20171126 Newdata/LargeSample20171126.RData")
 #load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20171126 Newdata/SmallSample20171126.RData")

# 0.2 Review the structure of results

length(restuls_conditions)  # 108 --> total number of cells (Note that there are two datasets, one for small sample one for large sample, and thus 108x2=216 cells)
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:7] --> 7 methods (the first column - i.e., the first method - is not useful, because it will be replaced 
                                       # by the population reliability; see simulation.R)

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:8] --> the average reliability for each sample of persons (columns are 7 methods). NOT USED
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:8] --> the SD. NOT USED.



############## 1. Research Question 1: Bias #########################################################

D_reliability <- matrix(NA, 108, 6) # this is for calculating the bias for the 6 change-score reliability estimates for each cell. 
Est_reliability <- matrix(NA, 108, 6) # this is useful for Research Question 2
for (c in 1:108){ # cth cell
  
  temp_recorder <- matrix(NA, 20, 6)
  temp_recorder_rel <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
     temp_recorder[l, ] <- colSums(restuls_conditions[[c]][[1]][[l]][, 2:7] - r_pop[c])  #cth cell, lth sample
     temp_recorder_rel[l, ] <- colSums(restuls_conditions[[c]][[1]][[l]][, 2:7])  
    }
  
  D_reliability[c, ] <- colSums(temp_recorder) / 1000
  Est_reliability[c, ] <- colSums(temp_recorder_rel) / 1000 # the average estimated reliability across 1000 item-score datasets. Useful for Q2 and Q4. 
}

Q1_Bias <- D_reliability

############# 2. Research Question 2: Precision #######################################################

Q2_Precision <- matrix(NA, 108, 6)
for (c in 1:108){ # cth cell
  
  temp_recorder <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
    temp_recorder[l, ] <- colSums(sweep(restuls_conditions[[c]][[1]][[l]][, 2:7], 2, Est_reliability[c, ], "-")^2)  #cth cell, lth sample
  }
  Q2_Precision[c, ] <- sqrt(colSums(temp_recorder) / 999)
 
}


############# 3. Research Question 3: Proportion of variance in bias due to sampling at Level 2 ##################################
SS_L2 <- matrix(NA, 108, 6)
SS_L1 <- matrix(NA, 108, 6)
SS_Total <- matrix(NA, 108, 6)
for (c in 1:108){ # cth cell
  
  SS_L2_temp <- matrix(NA, 20, 6)
  SS_L1_temp <- matrix(NA, 20, 6)
  SS_Total_temp <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
    SS_L2_temp[l, ] <- 50*(colSums(restuls_conditions[[c]][[1]][[l]][, 2:7])/50 - r_pop[c])^2
    SS_L1_temp[l, ] <- colSums(sweep(restuls_conditions[[c]][[1]][[l]][, 2:7], 2, colSums(restuls_conditions[[c]][[1]][[l]][, 2:7])/50, "-")^2)
    SS_Total_temp[l, ] <- colSums((restuls_conditions[[c]][[1]][[l]][, 2:7] - r_pop[c])^2)
  }
  
  SS_L2[c, ] <- colSums(SS_L2_temp)
  SS_L1[c, ] <- colSums(SS_L1_temp)
  SS_Total[c, ] <- colSums(SS_Total_temp)
}
((SS_L2 + SS_L1) - SS_Total)> .0000001  # this is to verify that SS_Total = SS_L1 + SS_L2, note that rounding errors exist.

proportion_L2_Total <- matrix(NA, 108, 6)
for(c in 1:108){
  for (m in 1:6){
    proportion_L2_Total[c, m] <- SS_L2[c, m]/SS_Total[c, m]
  }
}


############# 4. Research Question 4: Proportion of variance in precision due to sampling at Level 2  ############################
SS_L2_p <- matrix(NA, 108, 6)
SS_L1_p <- matrix(NA, 108, 6)
SS_Total_p <- matrix(NA, 108, 6)

for (c in 1:108){ # cth cell
  
  SS_L2_temp <- matrix(NA, 20, 6)
  SS_L1_temp <- matrix(NA, 20, 6)
  SS_Total_temp <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
    SS_L2_temp[l, ] <- 50*(colSums(restuls_conditions[[c]][[1]][[l]][, 2:7])/50 - Est_reliability[c, ])^2
    SS_L1_temp[l, ] <- colSums(sweep(restuls_conditions[[c]][[1]][[l]][, 2:7], 2, colSums(restuls_conditions[[c]][[1]][[l]][, 2:7])/50, "-")^2)
    SS_Total_temp[l, ] <- colSums((restuls_conditions[[c]][[1]][[l]][, 2:7] - r_pop[c])^2)
  }
  
  SS_L2_p[c, ] <- colSums(SS_L2_temp)
  SS_L1_p[c, ] <- colSums(SS_L1_temp)
  SS_Total_p[c, ] <- colSums(SS_Total_temp)
}

((SS_L2_p + SS_L1_p) - SS_Total_p)> .0000001 

proportion_L2_Total_p <- matrix(NA, 108, 6)
for(c in 1:108){
  for (m in 1:6){
    proportion_L2_Total_p[c, m] <- SS_L2_p[c, m]/SS_Total_p[c, m]
  }
}



####################### Summarize results in terms of 4 situations  ############################
# situation1 --> "Unidimensional + No Carryover"
# situation2 --> "Unidimensional + Carryover"
# situation3 --> "Multidimensional + No Carryover"
# situation4 --> "Multidimensional + Carryover"


test_length <- c(9, 21, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .1, .5) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
magnitude_sd <- c(sqrt(.14), sqrt(.5))  # .14 == small variance, .5 == large variance.
strongweak_carry <- c(0, 10, 1) #0 == no, 10 == strong, 1 == weak

num_condition <- length(test_length)*length(parallel_item)*length(correlated_facets)*length(magnitude_sd)*length(strongweak_carry)
conditions <- list()

p <- 1
for(i in 1:length(test_length)){
  for (j in 1:length(parallel_item)){
    for (k in 1:length(correlated_facets)){
      for(l in 1:length(magnitude_sd)){
        for(m in 1:length(strongweak_carry)){
          
          conditions[[p]] <- c(test_length[i], parallel_item[j], correlated_facets[k], magnitude_sd[l], strongweak_carry[m])
          p <- p + 1
        }
      }
    }
  }
}


df <- data.frame(matrix(unlist(conditions), nrow=num_condition, byrow = T))
colnames(df) <- c('test length', 'parallel item', 'correlated facets','maginute of sd', 'carry-over effects')


situ1 <- df$`correlated facets`==1 & df$`carry-over effects`==0
situ2 <- df$`correlated facets`==1 & df$`carry-over effects`!=0
situ3 <- df$`correlated facets`!=1 & df$`carry-over effects`==0
situ4 <- df$`correlated facets`!=1 & df$`carry-over effects`!=0


# summurize results regarding bias


write.csv(Q1_Bias[situ1, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ1_LargeSample.csv")
write.csv(Q1_Bias[situ2, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ2_LargeSample.csv")
write.csv(Q1_Bias[situ3, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ3_LargeSample.csv")
write.csv(Q1_Bias[situ4, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ4_LargeSample.csv")

write.csv(Q2_Precision[situ1, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ1_LargeSample.csv")
write.csv(Q2_Precision[situ2, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ2_LargeSample.csv")
write.csv(Q2_Precision[situ3, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ3_LargeSample.csv")
write.csv(Q2_Precision[situ4, ], file = "D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20180121 DataAnalysis/situ4_LargeSample.csv")


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




##################### EXTRA: not included in the paper, but can be useful

#### Extra 1: draw bias and precision plots at level 1

load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170122 rerun to record sum scores/results20170122.RData")
load("D:/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170222 results small sample/results20170222smallsample.RData")

#note run simulation <- 1, and then the analysis, then simulation <- 2....
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

######### Extra 2. plots for Bias and precision - LEVEL 2 ###################

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