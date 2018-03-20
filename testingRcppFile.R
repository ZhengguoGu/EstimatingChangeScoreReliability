#######################################################
### Testing file
#######################################################

#### This is to test a few functions in populationRel.cpp. 
Rcpp::sourceCpp("C:\\Users\\Zhengguo\\Documents\\EstimatingChangeScoreReliability\\populationRel.cpp")

#### Test  GRMc_1theta
set.seed(1)
theta <- rnorm(100000, 0, 1)
theta <- sort(theta)
num_items <- 10
itempar <- matrix(NA,num_items,5)
itempar[,1] <- runif(num_items,1.5,2.5)  # discrimination
avg_beta <- runif(num_items, -1.5, 1.5)
itempar[,2] <- avg_beta - .75
itempar[,3] <- avg_beta - .25
itempar[,4] <- avg_beta + .25
itempar[,5] <- avg_beta + .75

result <- matrix(0, 100000, 10)
for (n in 1:100000){
  
  result[n,] <- GRMc_1theta(theta[n], itempar[,1], itempar[,2:5])
  
}

library(ltm)

parameteres<-ltm::grm(data.matrix(result))  #pretest scores

est_par <- matrix(unlist(parameteres$coefficients), ncol = 5, byrow = TRUE)
est_par <- cbind(est_par[, 5], est_par[, -5])
plot(itempar[, 1], est_par[, 1])  #slope
plot(itempar[, 2], est_par[, 2])  #item-difficulty 1
plot(itempar[, 3], est_par[, 3])  #item-difficulty 2
plot(itempar[, 4], est_par[, 4])  #item-difficulty 3
plot(itempar[, 5], est_par[, 5])  #item-difficulty 4 


#### Test a vector - matrix operation in GRMc_1thetaMD()
# GRMc_1thetaMD is not very different from GRMc_1theta: the only difference is that 
# GRMc_1thetaMD requires a vector of theta (i.e., each person has a vector of theta, because of multidimensionality)
# What is to be tested here is "theta_BMD(_, i) = thetaMD - Idiff(_, i);" (see below). The reason why I test it is because
# I notice that unlike R, C++ does not allow me to do Vector minus Matrix operation (so I have to do it per column of matrix)
library(Rcpp)
cppFunction('NumericMatrix V_M(NumericVector thetaMD, NumericMatrix Idiff) {
            int ncat = Idiff.ncol();
            NumericMatrix theta_BMD(thetaMD.size(), ncat);
            for(int i = 0; i < ncat; ++i){
            theta_BMD(_, i) = thetaMD - Idiff(_, i);
            }
            return(theta_BMD);
            }')

theta <- 1:5
diff <- matrix(1:30, 5, 6)
V_M(theta, diff)
