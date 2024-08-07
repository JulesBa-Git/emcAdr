#include "RcppArmadillo.h"
#include <vector>
#include <numeric>
#include <iostream>
#include <sstream>
#include <algorithm>
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

bool hasExtension(const std::string& filename, const std::string& extension){
  size_t dotPos = filename.rfind('.');
  
  std::string fileExtension = dotPos == std::string::npos ? "" : filename.substr(dotPos);
  
  return fileExtension == extension;
}

std::vector<double> p_value_csv_file_size_k(const Rcpp::List& distribution_output,
                                            std::ifstream& input, int k){
  std::string line;
  std::vector<std::string> cellVec;
  cellVec.reserve(5);
  
  std::vector<double> solutions;
  solutions.reserve(100);
  
  Rcpp::Function compute_p_value = Rf_findFun(Rf_install("p_value_greater_than_empirical"),
                                              R_GlobalEnv);
  
  while(std::getline(input, line)){
    std::stringstream rowToAnalyze(line);
    std::string cell;
    
    int i = 0;
    //we only need the first 2 cells
    while(std::getline(rowToAnalyze, cell, ';')){
      cellVec.push_back(cell);
    }
    
    //if the cocktail is of the right size, we compute his p_value
    if(std::count(cellVec[1].begin(), cellVec[1].end(), ':') == k){
      solutions.push_back(Rcpp::as<double>(compute_p_value(distribution_output,
                                          std::stod(cellVec[0]),
                                          false)));
    }
    else{
      solutions.push_back(INT_MIN);
    }
    cellVec.clear();
  }
  return solutions;
}

// gérer le test sur l'extension fichier dans la fonction supérieur (celle 
// qui lance la boucle pour tous les k) dans le but d'ouvrir le buffer de fichier
// une seule fois 
//[[Rcpp::export]]
std::vector<double> p_value_of_genetic_size_k(const Rcpp::List& distribution_output, 
                                              const std::string& filename,int k
                                              ){
  std::vector<double> p_values;
  if(hasExtension(filename,".txt")){
    std::cout << ".txt file not supported for now \n";
  }
  else if(hasExtension(filename, ".csv")){
    std::ifstream ifstr(filename);
    if(!ifstr.is_open()){
      std::cerr << "the file " << filename << " has failed to open\n";
      return p_values;
    }
    p_values = p_value_csv_file_size_k(distribution_output, ifstr, k);
  }
  
  return p_values;
}

/*** R
library(emcAdr)
#to test (there is hard coded path)
treeATC <- read.csv("your/path/to/ATCtree")
patientsATC <- read.csv("your/path/to/testPatientATCList")
ATCtoNumeric(patientsATC,treeATC)
*/


