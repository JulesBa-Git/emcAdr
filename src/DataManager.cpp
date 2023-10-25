#include "RcppArmadillo.h"
#include <vector>
#include <numeric>
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

//'Convert the histogram returned by the DistributionApproximation function, to a real number ditribution
//'(that can be used in a test for example) 
//'
//'@param vec : distribution returned by the DistributionAproximationFunction
//'@export
// [[Rcpp::export]]
Rcpp::NumericVector histogramToDitribution(const std::vector<int>& vec){
  std::vector<double> returnedVec;
  returnedVec.reserve(std::accumulate(vec.begin(),vec.end(),0));
  int count;
  for(int i = 0 ; i < vec.size(); ++i){
    count = vec[i];
    for(int j = 0 ; j < count ; ++j){
      returnedVec.push_back(static_cast<double>(i)/10.0);
    }
  }
  return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(returnedVec));
}

// [[Rcpp::export]]
Rcpp::NumericVector incorporateOustandingRRToDistribution(const std::vector<double>& outstandingRR, int RRmax){
  std::vector<double> returnedVec;
  returnedVec.resize((RRmax*10)+1);
  
  for(const auto& rr: outstandingRR){
    int index;
    if(rr < RRmax){
      index= rr*10;
    }
    else{
      index = returnedVec.size()-1;
    }
    ++returnedVec[index];
  }
  
  return Rcpp::wrap(returnedVec);
}
/*** R
#to test (there is hard coded path)
treeATC <- read.csv("your/path/to/ATCtree")
patientsATC <- read.csv("your/path/to/testPatientATCList")
ATCtoNumeric(patientsATC,treeATC)
*/


