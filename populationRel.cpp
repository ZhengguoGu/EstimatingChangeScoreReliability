#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

int Carryover(int pre, int post, String eff){
  // eff == "S": strong effect
  // eff == "W": weak effect
  // eff == "N": no effect
  int post_out = 0;
  
  if(eff == "S"){
    if(pre - post < -1){
      post_out = pre + 1;
    }else if(pre - post > 1){
      post_out = pre - 1;
    }else{
      post_out = pre;
    }
  }else if(eff == "W"){
    if(pre - post < -1){
      post_out = post - 1;
    }else if(pre - post > 1){
      post_out = post + 1;
    }else{
      post_out = post;
    }
  }else if(eff == "N"){
    post_out = post;
  }
  
  return(post_out);
}

// [[Rcpp::export]]
List GRMc_1theta(double theta_pre, double theta_post, NumericVector Islope, NumericMatrix Idiff, String eff) {
  
  // Note that theta here is with respect to one person. 
  
  int nitems = Idiff.nrow();
  int ncat = Idiff.ncol();
  IntegerVector response_pre(nitems);
  IntegerVector response_post(nitems);
  IntegerVector response_cary(nitems);
  
  NumericMatrix numeritor_pre(nitems, ncat);
  NumericMatrix numeritor_post(nitems, ncat);  
  NumericMatrix Pstar_pre(nitems, ncat);
  NumericMatrix Pstar_post(nitems, ncat);
  
  
  // simulate pretest scores
  NumericVector randomnumber_pre = runif(nitems, 0.0, 1.0); 
  NumericMatrix Idiff_pre = theta_pre - Idiff;
  
   for(int j = 0; j < nitems; ++j){
    
    numeritor_pre(j, _) = exp(Islope(j) * Idiff_pre(j,_));
    
    for(int k = 0; k < ncat; ++k){
      Pstar_pre(j,k) = numeritor_pre(j,k)/(1.0+numeritor_pre(j,k));
      if(Pstar_pre(j,k) > randomnumber_pre(j)){
        response_pre(j) = response_pre(j) + 1;
      }
    }
  }
  
  
  // simulate posttest scores
  NumericVector randomnumber_post = runif(nitems, 0.0, 1.0); //its necessary to re-generate random numbers, if not, the errors are correlated.
  NumericMatrix Idiff_post = theta_post - Idiff;
   
   for(int j = 0; j < nitems; ++j){
     
     numeritor_post(j, _) = exp(Islope(j) * Idiff_post(j,_));
     
     for(int k = 0; k < ncat; ++k){
       Pstar_post(j,k) = numeritor_post(j,k)/(1.0+numeritor_post(j,k));
       if(Pstar_post(j,k) > randomnumber_post(j)){
         response_post(j) = response_post(j) + 1;
       }
     }
   }

  // introducing carry-over
  for(int l = 0; l < nitems; ++l){
    response_cary(l) = Carryover(response_pre(l), response_post(l), eff);
  }
  
  List output;
  output["response_pre"] = response_pre;
  output["response_post"] = response_post;
  output["response_carry"] = response_cary;
  
  return(output);
}


// [[Rcpp::export]]
List TrueScore(double t_pre, double t_post, NumericVector Islope, NumericMatrix Idiff, String eff, int NResample) {

  IntegerVector change_S(NResample);
  IntegerMatrix Pretest_m(NResample, Islope.size());
  IntegerMatrix Posttest_m(NResample, Islope.size());
  
  for (int m = 0; m < NResample; ++m){ // this is to calculate the expectation i.e., true score
    List A_result = GRMc_1theta(t_pre, t_post, Islope, Idiff, eff);
    IntegerVector pretest_vec = A_result["response_pre"];
    IntegerVector posttest_vac = A_result["response_carry"];
    Pretest_m(m, _) = pretest_vec;
    Posttest_m(m, _) = posttest_vac;
    int preSum = sum(pretest_vec);
    int postSum = sum(posttest_vac); // if no carry-over effect, then "response_carry" contains the post scores without carryover effect
    change_S(m) = postSum - preSum;
  }
  double true_change = mean(change_S);
  
  List B_result = GRMc_1theta(t_pre, t_post, Islope, Idiff, eff); // just generate 1 dataset, and ust it as the observed score.
  IntegerVector pretest_Obs = B_result["response_pre"];
  IntegerVector posttest_Obs = B_result["response_carry"];
  int change_Obs = sum(posttest_Obs) - sum(pretest_Obs);
  
  List OUTPUT;
  OUTPUT["change score"] = change_S;
  OUTPUT["pretest matrix"] = Pretest_m;
  OUTPUT["posttest matrix"] = Posttest_m;
  OUTPUT["true change"] = true_change;
  OUTPUT["observed change"] = change_Obs;
  return(OUTPUT);
  
  
}