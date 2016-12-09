#---------------------------------------------------------------
# Note that this function has been added to the simulation.R.
# Thus, this function won't be updated any longer. 
# See the file simulation.R
# 
# 2016/12/09
#---------------------------------------------------------------


test_length <- c(9, 18, 27, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .3, .5, .8) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
magnitude_sd <- c(sqrt(.14), sqrt(.5))
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


df <- data.frame(matrix(unlist(conditions), nrow=192, byrow = T))
colnames(df) <- c('test length', 'parallel item', 'correlated facets',
                  'maginute of sd', 'carry-over effects')
library(xlsx)
library(rJava)
write.csv(df, 'conditions.csv', sep=',')
