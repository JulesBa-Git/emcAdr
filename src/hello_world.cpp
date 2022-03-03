// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <iostream>
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

//'Fonction test doc
//'
//'@return void
//[[Rcpp::export]]
void HelloWorld(){
    std::cout << "Hello, World !" << std::endl;
}