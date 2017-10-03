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




for(sim in 1: num_condition){
  cond <- df[sim, ]

  ######## w.r.t. observed change scores ##########
  if(cond[3] == 1){
    # in this case, unidimensional
    if(cond[4] == sqrt(.14)){
      # in this case small change in theta: .14
      qpoint_pre <- Qpoints(n_person=1000, mu=0, sd=1, bd=3)
      qpoint_post <- Qpoints(n_person=1000, mu=0, sd=sqrt(1.14), bd=3)
      for(q in 1:1000){
        P_pre <- GRM_sim_1theta(theta = qpoint_pre[[1]][q], itempar =  ITEM_PAR[[5]])[[1]]
        w_pre <- Phi_X(P_pre)
        
        P_post <- GRM_sim_1theta(theta = qpoint_post[[1]][q], itempar =  ITEM_PAR[[5]])[[1]]
        w_post <- Phi_X(P_post)
        
        if(cond[5] == 10){
          aa <- Phi_D(w_pre, w_post, car_eff="strong")
        }
        
      }
      
      
      
    }
  }
  
  
}

