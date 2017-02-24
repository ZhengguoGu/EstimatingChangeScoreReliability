##########################################################################################
#
# This file contains the code for analyzing the results
#
##########################################################################################

# 1. load simulation results
#load('results20170122.RData')
load('results20170222smallsample.RData')

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
layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
plot(pop_re[, 2], xlim = c(0, 108), ylim = c(0,1), xlab = "108 cells", ylab = "(Estimated) Change-Score Reliability", col='red', pch=8)
points(pop_re[, 3], pch=0)
points(pop_re[, 4], pch=1)
points(pop_re[, 5], pch=2)
points(pop_re[, 6], pch=15, col='blue')
points(pop_re[, 7], pch=19, col='blue')
points(pop_re[, 8], pch=17, col='blue')

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

#############--------------------------------------
#4. standard deviation of true change-score reliability and the biasness
# of each estimates 


sdTrue <- array()
for (cel in 1:108){
  sdTrue[cel] <- sd(restuls_conditions[[cel]][[2]][, 2]) # true reliability's SD 
}

True1sdmax <- pop_re[, 2] + sdTrue #True change score reliability + 1SD
True1sdmin <- pop_re[, 2] - sdTrue #True change score reliability 1 1SD

True3sdmax <- pop_re[, 2] + 3*sdTrue
True3sdmin <- pop_re[, 2] - 3*sdTrue

prop1sd <- array()
prop3sd <- array()
for(i in 3:8){
  prop1sd[i-2] <- sum(pop_re[, i] <= True1sdmax & pop_re[, i] >= True1sdmin) #tells us how many cells out of 108 cells...
  prop3sd[i-2] <- sum(pop_re[, i] <= True3sdmax & pop_re[, i] >= True3sdmin)
}


# 5. sample level: for each cell, we draw plots reliabilities 

############## all 6 methods are ploted in one pic


allplots <- list()

for(cel in 1:108){

  y_min <- 1
  y_max <- 1
  
  for(i in 3:8){
    
    y_min <- min(y_min, min(restuls_conditions[[cel]][[2]][, i] - restuls_conditions[[cel]][[3]][, i]))
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

#6. pick a few representative plots to be included in the paper. the rest are in the supplementary matrial
# We need to rerun the plots because this time ylim must be the same so that we can compare. 

cel <- 24
y_min <- 0
y_max <- 1

layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
xtitle <- sprintf("20 samples from the population: Cell %d; Small sample size",cel)
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

dev.off()




#############--------------------------------------
#7.  given the 108 cells, which cells generate negative reliabilites?

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

#8. ANOVAs. The idea is to see how different conditions influence the bias of reliability estimates. 
#Note that for each cell (in total 108), there are 20 samples; and for each sample, there are 50 reponse datasets. For each cell,
#I randomly pick one sample, and combine the 108 samples together for ANOVAs.

load('results20170122.RData') #large sample case

index_sample <- sample(1:20, 108, replace=TRUE) 

reANOVA <- matrix(NA, 1, 13)
for(cel in 1:108){
  
  reSample <- restuls_conditions[[cel]][[1]][[index_sample[cel]]] 
  reSample <- cbind(reSample, t(matrix(rep(df[cel,], 50), 5, 50)))
  reANOVA <- rbind(reANOVA, reSample)

}
reANOVA <- reANOVA[-1,]
save(reANOVA, file = "reANOVAlarge.RData")

load('results20170222smallsample.RData') #small sample case

index_sample <- sample(1:20, 108, replace=TRUE) 

reANOVAs <- matrix(NA, 1, 13)
for(cel in 1:108){
  
  reSample <- restuls_conditions[[cel]][[1]][[index_sample[cel]]] 
  reSample <- cbind(reSample, t(matrix(rep(df[cel,], 50), 5, 50)))
  reANOVAs <- rbind(reANOVAs, reSample)
  
}
reANOVAs <- reANOVAs[-1,]
save(reANOVAs, file = "reANOVAsmall.RData")

load("~/EstimatingChangeScoreReliability/reANOVAlarge.RData")
load("~/EstimatingChangeScoreReliability/reANOVAsmall.RData")
names(df)

reANOVA <- as.data.frame(reANOVA)
names(reANOVA) <- c("useless", "true", "trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                    "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects")
reANOVA$size <- "large"
reANOVAs <- as.data.frame(reANOVAs)
names(reANOVAs) <- c("useless", "true", "trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                    "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects")
reANOVAs$size <- "small"

combANOVA <- rbind(reANOVA,reANOVAs) # this is the raw data.. 
#calculate the bias for each row, this is true - estimated
biasANOVA <- matrix(NA, 10800, 12)
biasANOVA[, 1] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`trad alpha`)
biasANOVA[, 2] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`trad lambda2`)
biasANOVA[, 3] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`trad lambda4`)
biasANOVA[, 4] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`item alpha`)
biasANOVA[, 5] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`item lambda2`)
biasANOVA[, 6] <- as.numeric(combANOVA$true) - as.numeric(combANOVA$`item lambda4`)
# note that biasANOVA[, 7:12] <- combANOVA[, 9:14] is not recommended! The resulting file will 31 GB. (although I dont know the reason why)
biasANOVA <- as.data.frame(biasANOVA)
biasANOVA$V7 <- factor(as.numeric(combANOVA$`test length`))
biasANOVA$V8 <- factor(as.numeric(combANOVA$`parallel item`))
biasANOVA$V9 <- factor(as.numeric(combANOVA$`correlated facets`))
biasANOVA$V10 <- factor(as.numeric(combANOVA$`maginute of sd`))
biasANOVA$V11 <- factor(as.numeric(combANOVA$`carry-over effects`))
biasANOVA$V12 <- factor(combANOVA$size)
names(biasANOVA) <- c("trad alpha", "trad lambda2", "trad lambda4", "item alpha", "item lambda2", "item lambda4",
                      "test length", "parallel item", "correlated facets",  "maginute of sd", "carry-over effects", "size")
save(biasANOVA, file="biasANOVAdata.RData")

anova(lm(biasANOVA$'trad alpha' ~ biasANOVA$'test length' + 
                                  biasANOVA$'parallel item' + 
                                  biasANOVA$'correlated facets' +
                                  biasANOVA$'maginute of sd' + 
                                  biasANOVA$'carry-over effects' + 
                                  biasANOVA$'size' +
                                  biasANOVA$'test length'*biasANOVA$'size' +
                                  biasANOVA$'parallel item'*biasANOVA$'size' +
                                  biasANOVA$'correlated facets'*biasANOVA$'size' +
                                  biasANOVA$'maginute of sd'*biasANOVA$'size'  +
                                  biasANOVA$'carry-over effects' *biasANOVA$'size'))

fit <- lm(biasANOVA$'trad alpha' ~ biasANOVA$'test length' + 
            biasANOVA$'parallel item' + 
            biasANOVA$'correlated facets' +
            biasANOVA$'maginute of sd' + 
            biasANOVA$'carry-over effects' + 
            biasANOVA$'size' +
            biasANOVA$'test length'*biasANOVA$'size' +
            biasANOVA$'parallel item'*biasANOVA$'size' +
            biasANOVA$'correlated facets'*biasANOVA$'size' +
            biasANOVA$'maginute of sd'*biasANOVA$'size'  +
            biasANOVA$'carry-over effects' *biasANOVA$'size')

summary(fit)
#The following are not used, but might be useful?
#groupindex <- data.frame(biasANOVA[, 7:12])
#groupindexNonDup <- groupindex[!duplicated(groupindex, margin=2),]
#dim(groupindexNonDup)
#splitANOVA <- split(biasANOVA, biasANOVA[, 7:12])
#length(splitANOVA)
