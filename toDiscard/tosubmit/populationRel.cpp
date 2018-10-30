#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector GRMc_1theta(double theta, NumericVector Islope, NumericMatrix Idiff){
  
  // Note that theta here is with respect to one person. 
  
  int nitems = Idiff.nrow();
  int ncat = Idiff.ncol();
  IntegerVector response(nitems);
  
  NumericMatrix numeritor(nitems, ncat);
  NumericMatrix Pstar(nitems, ncat);
  
  NumericVector randomnumber = runif(nitems, 0.0, 1.0); 
  NumericMatrix Idiff_theta = theta - Idiff;
  
  for(int j = 0; j < nitems; ++j){
    
    numeritor(j, _) = exp(Islope(j) * Idiff_theta(j,_));
    
    for(int k = 0; k < ncat; ++k){
      Pstar(j,k) = numeritor(j,k)/(1.0+numeritor(j,k));
      if(Pstar(j,k) > randomnumber(j)){
        response(j) = response(j) + 1;
      }
    }
  }
  return(response);
}



// [[Rcpp::export]]
IntegerVector GRMc_1thetaMD(NumericVector thetaMD, NumericVector Islope, NumericMatrix Idiff) {
  // this is for multidimentional theta (several sub-attributes subsumed under one factor)
  // Note: theta is a vector, each element in the vector correspondes to an item.
  
  int nitems = Idiff.nrow();
  int ncat = Idiff.ncol();
  NumericMatrix numeritorMD(nitems, ncat);
  NumericMatrix PstarMD(nitems, ncat);
  IntegerVector responseMD(nitems);
  
  
  // simulate pretest scores
  NumericVector randomnumberMD = runif(nitems, 0.0, 1.0); 
  NumericMatrix theta_BMD(nitems, ncat); 
  for(int i = 0; i < ncat; ++i){
    theta_BMD(_, i) = thetaMD - Idiff(_, i);  //this step has been tested, it works. In the future, I should use RcppArmadillo to handle matrix operations
  }
  for(int j = 0; j < nitems; ++j){
    
    numeritorMD(j, _) = exp(Islope(j) * theta_BMD(j,_));
    for(int k = 0; k <ncat; ++k){
      PstarMD(j,k) = numeritorMD(j,k)/(1.0+numeritorMD(j,k));
      if(PstarMD(j,k) > randomnumberMD(j)){
        responseMD(j) = responseMD(j) + 1;
      }
    }
  }
  
 return(responseMD);
}


// [[Rcpp::export]]

IntegerVector Carryover(IntegerVector pre, IntegerVector post, String eff){
  // eff == "S": strong effect
  // eff == "W": weak effect
  // eff == "N": no effect
  
  int n_items = pre.size();
  
  IntegerVector post_out(n_items);
  
  if(eff == "S"){
    for (int k = 1; k < n_items; ++k){
      if(pre(k) - post(k) < -1){
        post_out(k) = pre(k) + 1;
      }else if(pre(k) - post(k) > 1){
        post_out(k) = pre(k) - 1;
      }else{
        post_out(k) = pre(k);
      }
    }
    
  }else if(eff == "W"){
    for (int k = 1; k < n_items; ++k){
      if(pre(k) - post(k) < -1){
        post_out(k) = post(k) - 1;
      }else if(pre(k) - post(k) > 1){
        post_out(k) = post(k) + 1;
      }else{
        post_out(k) = post(k);
      }
    }
  }else if(eff == "N"){
    post_out = post;
  }
  
  return(post_out);
}

