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
      patientI.push_back(posTree);
    }
    patientI.shrink_to_fit();
    newPatientATC.push_back(patientI);
    
  }
  patients["patientATC"] = newPatientATC;
  
}


/*** R
#for test there is hard coded path
treeATC <- read.csv("your/path/to/ATCtree")
patientsATC <- read.csv("your/path/to/testPatientATCList")
ATCtoNumeric(patientsATC,treeATC)
*/

