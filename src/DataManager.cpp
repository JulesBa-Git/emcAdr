#include "RcppArmadillo.h"
#include <vector>
#include <iostream>
using Rcpp::DataFrame;


// [[Rcpp::depends(RcppArmadillo)]]
//'Convert ATC Code for each patients to the corresponding DFS number of the ATC tree 
//'
//'@param tree : ATC tree (we assume that there is a column 'ATCCode' )
//'@param patients : patients observations, for each patient we got a string containing every medication he takes/took
//'@export
// [[Rcpp::export]]
void ATCtoNumeric(DataFrame& patients,const DataFrame& tree) {
  std::vector<std::string> cppTree= tree["ATCCode"];
  std::vector<std::string> patientATC = patients["patientATC"];
  std::vector<std::vector<int>> newPatientATC;
  
  newPatientATC.reserve(patientATC.size());
  
  std::string codeI;
  std::string delimiter = " ";
  int posOcc = 0,posTree = 0;
  std::vector<int> patientI(0);
  for(std::string& code : patientATC){
    patientI.clear();
    patientI.reserve(3);
    
    codeI = code;
    posOcc = 0;
    posTree = 0;
    
    while(posOcc > -1){
      posOcc = codeI.find(delimiter);
      
      //token is the ATC code being converted to numeric 
      std::string token = codeI.substr(0,posOcc);
      //new code without the token (because it is already converted)
      codeI = codeI.substr(posOcc+1);
      //first condition should be useless but could detect error -> when the ATC code is not found in the ATC tree
      while(posTree < cppTree.size() && token != cppTree[posTree]){
        ++posTree;
      }
      
      if(posTree == cppTree.size()){
        std::cerr<<"error : a patient take a medication that is not in the tree" << '\n';
        return;
      }
      //+1 because of the cpp indexes (starting at 0)
      patientI.push_back(posTree+1);
    }
    patientI.shrink_to_fit();
    newPatientATC.push_back(patientI);
    
  }
  patients["patientATC"] = newPatientATC;
  
}

//' test data transformation for histogram plot (in)
//' 
//' @param RRDistribution : the RR distribution Vector (given by the EMC algorithm)
//' @param f : the function used to plot, hist by default
//[[Rcpp::export]]
void frequencyHist(const Rcpp::IntegerVector& RRDistribution, Rcpp::Function f = Rcpp::Function("hist")){
  int size = 0;
  double RR = 0.9;
  for(auto tmp : RRDistribution)
    size+=tmp;
  
  std::vector<double> ret{};
  ret.reserve(size);
  
  for(auto tmp : RRDistribution){
    for(int i = 0 ; i < tmp ; ++i){
      ret.push_back(RR);
    }
    RR+=0.1;
  }
  
  Rcpp::NumericVector rRet;
  rRet = ret;
  f(rRet);
}

/*** R
#to test (there is hard coded path)
treeATC <- read.csv("your/path/to/ATCtree")
patientsATC <- read.csv("your/path/to/testPatientATCList")
ATCtoNumeric(patientsATC,treeATC)
*/


