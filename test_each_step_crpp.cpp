#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix testMD(NumericVector thetaMD, NumericVector Islope, NumericMatrix Idiff) {
  
  int nitems = Idiff.nrow();
  int ncat = Idiff.ncol();
  NumericMatrix theta_BMD(nitems, ncat); 
  NumericMatrix numeritorMD(nitems, ncat);
  for(int i = 0; i < ncat; ++i){
    theta_BMD(_, i) = thetaMD - Idiff(_, i);  //this step has been tested, it works. In the future, I should use RcppArmadillo to handle matrix operations
  }
  
  print(theta_BMD);
  
  for(int j = 0; j < nitems; ++j){
    
    numeritorMD(j, _) = exp(Islope(j) * theta_BMD(j,_));

  }
  
  return numeritorMD;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
thetaMD <- 1:5
Idiff <- matrix(1:20, 5, 4)
Islope <- c(2, 4, 6, 8, 10)
testMD(thetaMD, Islope, Idiff)
*/
