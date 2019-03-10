#########################################################################################################
#
# This file contains the code for analyzing the results
#
#
# Zhengguo Gu, Tilburg University
#########################################################################################################




############### 0. Check data #######################################################################


# 0.1 Load data: Note: there are in total 3 datasets, one for the population reliability, one for N=1000 (estimated rel), and one for N=100 (estimate rel).
# To replicate the result, one has to run the code for N=1000 and for N=100 seperatly. That is, after loading the N=1000 and the population rel dataset, run the code below and save the results. 
# Afterwards, one load the N=100 dataset (i.e., N=1000 data will be replaced with N=100) and run the code below again. 

load("PopulationRel20181030/PopulationRel.RData")
#load("EstimatingReliability20181024_estimatedRel/LargeSample20181025.RData")
load("EstimatingReliability20181024_estimatedRel/SmallSample20181026.RData")


# 0.2 Review the structure of results

length(restuls_conditions)  # 180 --> total number of cells (Note that there are two datasets, one for small sample one for large sample, and thus 180x2=360 cells)
length(restuls_conditions[[1]]) # 3 (lists)--> the first list contains all the reliability estimates for the 20 samples of persons
# the second list contains the average reliability for each sample of persons (thus, 20 rows)
# the third list contains the SD (might not be useful at this moment)

length(restuls_conditions[[1]][[1]]) # 20 (lists) --> each list contains the reliability estimates for 50 samples of responses
str(restuls_conditions[[1]][[1]][[1]]) # num [1:50, 1:6] --> (6 estimates; 3 for the traditional method and 3 for the item-level method)

str(restuls_conditions[[1]][[2]]) # num [1:20, 1:6] --> the average reliability for each sample of persons. NOT USED
str(restuls_conditions[[1]][[3]]) # num [1:20, 1:6] --> the SD. NOT USED.



############## 1. Data prep: Bias #########################################################

D_reliability <- matrix(NA, 180, 6) # this is for calculating the bias for the 6 change-score reliability estimates for each cell. 
Est_reliability <- matrix(NA, 180, 6) # this is useful for later
for (c in 1:180){ # cth cell
  
  temp_recorder <- matrix(NA, 20, 6)
  temp_recorder_rel <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
     temp_recorder[l, ] <- colSums(restuls_conditions[[c]][[1]][[l]][, 1:6] - r_pop[c])  #cth cell, lth sample
     temp_recorder_rel[l, ] <- colSums(restuls_conditions[[c]][[1]][[l]][, 1:6])  
    }
  
  D_reliability[c, ] <- colSums(temp_recorder) / 1000
  Est_reliability[c, ] <- colSums(temp_recorder_rel) / 1000 # the average estimated reliability across 1000 item-score datasets. 
}

Q1_Bias <- D_reliability

############# 2. Data prep: Precision #######################################################

Q2_Precision <- matrix(NA, 180, 6)
for (c in 1:180){ # cth cell
  
  temp_recorder <- matrix(NA, 20, 6)
  for(l in 1:20){ # lth sample
    temp_recorder[l, ] <- colSums(sweep(restuls_conditions[[c]][[1]][[l]][, 1:6], 2, Est_reliability[c, ], "-")^2)  #cth cell, lth sample
  }
  Q2_Precision[c, ] <- sqrt(colSums(temp_recorder) / 999)
 
}

#Q1_Bias_N1000 <- Q1_Bias
#Q2_Precision_N1000 <- Q2_Precision
#save(Q1_Bias_N1000, Q2_Precision_N1000, file = "Bias_and_precision_N1000.RData") 

Q1_Bias_N100 <- Q1_Bias
Q2_Precision_N100 <- Q2_Precision
save(Q1_Bias_N100, Q2_Precision_N100, file = "Bias_and_precision_N100.RData")



##############  3. preparing the data for ploting BIAS
load(file = "Bias_and_precision_N1000.RData")
load(file = "Bias_and_precision_N100.RData")

colnames(Q1_Bias_N1000) <- c("trad_alpha", "trad_l2", "trad_l4", "item_alpha", "item_l2", "item_l4")
colnames(Q1_Bias_N100) <- c("trad_alpha", "trad_l2", "trad_l4", "item_alpha", "item_l2", "item_l4")
colnames(Q2_Precision_N1000) <- c("trad_alpha", "trad_l2", "trad_l4", "item_alpha", "item_l2", "item_l4")
colnames(Q2_Precision_N100) <- c("trad_alpha", "trad_l2", "trad_l4", "item_alpha", "item_l2", "item_l4")


test_length <- c(9, 21, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .1, .5) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
magnitude_sd <- c(sqrt(.14), sqrt(.5))  # .14 == small variance, .5 == large variance
strongweak_carry <- c(0, 25, 50, 125, 150) #0:"no carryover effect", 25:"25% of persons showing weak", 50:"50% showing weak", 125: "25% strong", 150: "50% strong"; 
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
colnames(df) <- c('test length', 'parallel item', 'correlated facets',
                  'maginute of sd', 'carry-over effects')


Q1_Bias_N1000 <- cbind(1000, Q1_Bias_N1000)
colnames(Q1_Bias_N1000)[1] <- "sample size"
Q1_Bias_N100 <- cbind(100, Q1_Bias_N100)
colnames(Q1_Bias_N100)[1] <- "sample size"

BIAS_data <- rbind(cbind(df, Q1_Bias_N1000), cbind(df, Q1_Bias_N100))

write.table(BIAS_data, file = "BIAS_data.csv", sep=";", row.names=FALSE) # I used SPSS 24 to obtain the plots.



################ 4. preparing the data for ploting for precision


Q2_Precision_N1000 <- cbind(1000, Q2_Precision_N1000)
colnames(Q2_Precision_N1000)[1] <- "sample size"
Q2_Precision_N100 <- cbind(100, Q2_Precision_N100)
colnames(Q2_Precision_N100)[1] <- "sample size"

Precision_data <- rbind(cbind(df, Q2_Precision_N1000), cbind(df, Q2_Precision_N100))

write.table(Precision_data, file = "Precision_data.csv", sep=";", row.names=FALSE) # I used SPSS 24 to obtain the plots.

############### 5. plots for BIAS  (added on March 09, 2019 on mac)   #################

library(reshape2)
library(ggplot2)
library(latex2exp)
library(gridExtra)


bias_result <- read.table(file = "BIAS_data.csv", header = T, sep = ";")
bias_result$carry.over.effects <- as.factor(bias_result$carry.over.effects)


# plot: carry-over effects
dat_temp <- melt(bias_result,id.vars="carry.over.effects", measure.vars=c("trad_alpha", 
                                                                        "trad_l2",
                                                                        "trad_l4", 
                                                                        "item_alpha", 
                                                                        "item_l2", 
                                                                        "item_l4"))

p <- ggplot(dat_temp) +
  geom_boxplot(aes(x=carry.over.effects, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_discrete(labels=c("No effect", "25% weak", 
         "50% weak", "25% strong",
         "50% strong")) +
  labs(x = "Carry-Over Effects", y = "Bias") +
  theme(legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 
p

# plot: test length
bias_result$test.length <- as.factor(bias_result$test.length)
dat_temp_noCarry <- melt(bias_result[bias_result$carry.over.effects == 0, ],id.vars="test.length", measure.vars=c("trad_alpha", 
                                                                                                                  "trad_l2",
                                                                                                                  "trad_l4", 
                                                                                                                  "item_alpha", 
                                                                                                                  "item_l2", 
                                                                                                                  "item_l4"))
dat_temp_Carry <- melt(bias_result[bias_result$carry.over.effects != 0, ],id.vars="test.length", measure.vars=c("trad_alpha", 
                                                                                                                  "trad_l2",
                                                                                                                  "trad_l4", 
                                                                                                                  "item_alpha", 
                                                                                                                  "item_l2", 
                                                                                                                  "item_l4"))
p1 <- ggplot(dat_temp_noCarry) +
  geom_boxplot(aes(x=test.length, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(a) Without carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("9 items", "21 items", "36 items")) +
  labs(x = "Test Length", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

p2 <- ggplot(dat_temp_Carry) +
  geom_boxplot(aes(x=test.length, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(b) With carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("9 items", "21 items", "36 items")) +
  labs(x = "Test Length", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

grid.arrange(
  p1, p2,
  nrow = 2
) 

# item parameters
bias_result$parallel.item <- as.factor(bias_result$parallel.item)
dat_temp_noCarry <- melt(bias_result[bias_result$carry.over.effects == 0, ],id.vars="parallel.item", measure.vars=c("trad_alpha", 
                                                                                                                  "trad_l2",
                                                                                                                  "trad_l4", 
                                                                                                                  "item_alpha", 
                                                                                                                  "item_l2", 
                                                                                                                  "item_l4"))
dat_temp_Carry <- melt(bias_result[bias_result$carry.over.effects != 0, ],id.vars="parallel.item", measure.vars=c("trad_alpha", 
                                                                                                                "trad_l2",
                                                                                                                "trad_l4", 
                                                                                                                "item_alpha", 
                                                                                                                "item_l2", 
                                                                                                                "item_l4"))
p1 <- ggplot(dat_temp_noCarry) +
  geom_boxplot(aes(x=parallel.item, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(a) Without carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("Not identical", "Identical")) +
  labs(x = "Item Parameters", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

p2 <- ggplot(dat_temp_Carry) +
  geom_boxplot(aes(x=parallel.item, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(b) With carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("Not identical", "Identical")) +
  labs(x = "Item Parameters", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

grid.arrange(
  p1, p2,
  nrow = 2
) 


# dimensionality of theta
bias_result$correlated.facets <- as.factor(bias_result$correlated.facets)
dat_temp_noCarry <- melt(bias_result[bias_result$carry.over.effects == 0, ],id.vars="correlated.facets", measure.vars=c("trad_alpha", 
                                                                                                                    "trad_l2",
                                                                                                                    "trad_l4", 
                                                                                                                    "item_alpha", 
                                                                                                                    "item_l2", 
                                                                                                                    "item_l4"))
dat_temp_Carry <- melt(bias_result[bias_result$carry.over.effects != 0, ],id.vars="correlated.facets", measure.vars=c("trad_alpha", 
                                                                                                                  "trad_l2",
                                                                                                                  "trad_l4", 
                                                                                                                  "item_alpha", 
                                                                                                                  "item_l2", 
                                                                                                                  "item_l4"))
p1 <- ggplot(dat_temp_noCarry) +
  geom_boxplot(aes(x=correlated.facets, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(a) Without carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("3-D, corr=0.1", "3-D, corr=0.5", "Unidimensional")) +
  labs(x = expression(paste("Dimensionality of ", theta)), y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

p2 <- ggplot(dat_temp_Carry) +
  geom_boxplot(aes(x=correlated.facets, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(b) With carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("3-D, corr=0.1", "3-D, corr=0.5", "Unidimensional")) +
  labs(x = expression(paste("Dimensionality of ", theta)), y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

grid.arrange(
  p1, p2,
  nrow = 2
) 

# magnitude of variance
bias_result$maginute.of.sd <- as.factor(bias_result$maginute.of.sd)
dat_temp_noCarry <- melt(bias_result[bias_result$carry.over.effects == 0, ],id.vars="maginute.of.sd", measure.vars=c("trad_alpha", 
                                                                                                                        "trad_l2",
                                                                                                                        "trad_l4", 
                                                                                                                        "item_alpha", 
                                                                                                                        "item_l2", 
                                                                                                                        "item_l4"))
dat_temp_Carry <- melt(bias_result[bias_result$carry.over.effects != 0, ],id.vars="maginute.of.sd", measure.vars=c("trad_alpha", 
                                                                                                                      "trad_l2",
                                                                                                                      "trad_l4", 
                                                                                                                      "item_alpha", 
                                                                                                                      "item_l2", 
                                                                                                                      "item_l4"))
p1 <- ggplot(dat_temp_noCarry) +
  geom_boxplot(aes(x=maginute.of.sd, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(a) Without carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("Small variance", "Large variance")) +
  labs(x = expression(paste("Magnitude of Variance of ", theta, " Change")), y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

p2 <- ggplot(dat_temp_Carry) +
  geom_boxplot(aes(x=maginute.of.sd, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(b) With carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("Small variance", "Large variance")) +
  labs(x = expression(paste("Magnitude of Variance of ", theta, " Change")), y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

grid.arrange(
  p1, p2,
  nrow = 2
) 

# sample size
bias_result$sample.size<- as.factor(bias_result$sample.size)
dat_temp_noCarry <- melt(bias_result[bias_result$carry.over.effects == 0, ],id.vars="sample.size", measure.vars=c("trad_alpha", 
                                                                                                                     "trad_l2",
                                                                                                                     "trad_l4", 
                                                                                                                     "item_alpha", 
                                                                                                                     "item_l2", 
                                                                                                                     "item_l4"))
dat_temp_Carry <- melt(bias_result[bias_result$carry.over.effects != 0, ],id.vars="sample.size", measure.vars=c("trad_alpha", 
                                                                                                                   "trad_l2",
                                                                                                                   "trad_l4", 
                                                                                                                   "item_alpha", 
                                                                                                                   "item_l2", 
                                                                                                                   "item_l4"))
p1 <- ggplot(dat_temp_noCarry) +
  geom_boxplot(aes(x=sample.size, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(a) Without carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("N=100", "N=1,000")) +
  labs(x = "Sample Size", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

p2 <- ggplot(dat_temp_Carry) +
  geom_boxplot(aes(x=sample.size, y=value, color=variable)) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("(b) With carry-over effects") +
  ylim(-2, 0.5) +
  scale_x_discrete(labels=c("N=100", "N=1,000")) +
  labs(x = "Sample Size", y = "Bias") +
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) +
  scale_color_manual(labels = c(expression(r[alpha]), expression(r[lambda*2]), expression(r[lambda*4]), 
                                expression(a[D]), expression(l[2*D]), expression(l[4*D])), 
                     values = c("red", "green", "blue", "darkred", "darkgreen", "darkblue")) 

grid.arrange(
  p1, p2,
  nrow = 2
) 
