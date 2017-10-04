#########################################################
#
# reliability of change scores at the population level
# via numerical integration
#
# Zhengguo Gu
#
##########################################################

set.seed(110)

#------------------------------------------------------------------------------ 
# parameteres (they are directly copied from "simularion.R)
#------------------------------------------------------------------------------

test_length <- c(9, 21, 36)
parallel_item <- c(1, 0) # 1== yes, 0 == no
correlated_facets <- c(1, .6, .8) #if == 1, then dimension of theta is 1, otherwise 3 dimensions
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


df <- data.frame(matrix(unlist(conditions), nrow=num_condition, byrow = T))
colnames(df) <- c('test length', 'parallel item', 'correlated facets',
                  'maginute of sd', 'carry-over effects')

load("D:/TilburgOffice/Dropbox/tilburg office/Research Individual change/Project 3 - item difference scores/20170929 itemPar/itemparLarge.RData")
load("/Users/zhengguogu/Dropbox/Tilburg office/Research Individual change/Project 3 - item difference scores/20170929 itemPar/itemparLarge.RData")  #macbook pro
n_person <- 1000  #! need to manually set n_person = 1000 or 100, which corresponds to itempar data.


rel <- array()  #reliability at the population level

for(sim in 1: num_condition){
  cond <- df[sim, ]

  ######## w.r.t. observed change scores ##########
  if(cond[3] == 1){
    # in this case, unidimensional
    if(cond[4] == sqrt(.14)){
      # in this case small change in theta: .14
      qbipoints <- Qbipoints(n_person=1000, mu=c(0,0), sigma=matrix(c(1, 0, 0, .14), 2, 2), bd=3)
      theta_pre <-  qbipoints[[1]][, 1]
      theta_post <- rowSums(qbipoints[[1]])  #because Qbipoints negerates quodrature points for pretest and change
      
      Expect_D <- array()
      Expect_D2 <- array()
      for(q in 1:1000){
        P_pre <- GRM_sim_1theta(theta = theta_pre[q], itempar =  ITEM_PAR[[5]])[[1]]
        w_pre <- Phi_X(P_pre)
        
        P_post <- GRM_sim_1theta(theta = theta_post[q], itempar =  ITEM_PAR[[5]])[[1]]
        w_post <- Phi_X(P_post)
        
        if(cond[5] == 10){
          car_eff <- "strong"
        }else if(cond[5] == 1){
          car_eff <- "weak"
        }else{
          car_eff <- "no"
        }
        
        observed_D <- Phi_D(w_pre, w_post, car_eff=car_eff)[[2]]
        Expect_D[q] <- sum(observed_D[1,]*observed_D[2,])
        observed_D2 <- observed_D[1,] ^ 2
        Expect_D2[q] <- sum(observed_D2*observed_D[2,])
        
      }
      
      Expect_Dpop_obs <- sum(Expect_D*qbipoints[[2]])
      Expect_D2pop_obs <- sum(Expect_D2*qbipoints[[2]])
      Var_Dpop_obs <- Expect_D2pop - Expect_Dpop^2
      
      Expect_Dpop_true <- sum(Expect_D*qbipoints[[2]])
      Expect_D2pop_true <- sum(Expect_D^2*qbipoints[[2]])
      Var_D2pop_true <- Expect_D2pop_true - Expect_Dpop_true^2
      
      rel[sim] <- Var_Dpop_obs/Var_D2pop_true
      
    }
  }
  
  
}

