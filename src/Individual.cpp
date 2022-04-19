#include "Individual.h"
// [[Rcpp::depends(RcppArmadillo)]]

Individual::Individual(const std::vector<int>& medic, double temperature) :
  medications_{medic}, temperature_{temperature}
  {}

void Individual::printMedications() const{
  for(const int& med : medications_){
    std::cout << med<< ' ';
  }
  std::cout << '\n';
}
void Individual::printTemperature() const{
  std::cout << temperature_ << '\n';
} 

bool Individual::matches(const std::vector<int>& observation, const std::vector<int>& upperBound) const{
  int idx;
  bool inIt;

  for(const int& medIdx : medications_){
    inIt = false;
    idx = 0;
    
    while (idx < observation.size() && !inIt){
      //medix-1 because R index starts at 1 and the ATC tree is pretreated with R (will be changed ?)
      if(observation[idx] >= medIdx && observation[idx] <= upperBound[medIdx]-1){
        inIt = true;
      }else{
        ++idx;
      }
      
    }
    // if this condition is true, then no observation matched a medication of the individual, we return false
    if(idx == observation.size())
      return false;
  }

  return true;
  
}

double Individual::computeRR(const Rcpp::List& medications,const Rcpp::LogicalVector& ADR
                               , const Rcpp::DataFrame& ATCtree) const{
  int yesInDelt = 0, noInDelt = 0;
  int yesNotInDelt = 0, noNotInDelt = 0;
  double sumInDelt, sumNotInDelt;
  std::vector<int> upperBound = ATCtree["upperBound"];
  
  for(int i = 0; i < medications.size() ; ++i){
    //does the i-th observations is included in the individual ?
    bool isInDelta = this->matches(medications[i], upperBound);
    if(isInDelta){
      if(ADR[i]){
        ++yesInDelt;
      }else{
        ++noInDelt;
      }
    }
    else{
      if(ADR[i]){
        ++yesNotInDelt;
      }else{
        ++noNotInDelt;
      }
    }
    
  }
  sumInDelt = yesInDelt + noInDelt;
  sumNotInDelt = yesNotInDelt + noNotInDelt;
  //during the test the denominator was frequently 0 so I make it 1 if it is 0 
  //because if the denominator is 0 the numerator has to be 0 so the result would be 0 
  //no matter the denominator, we have 
  sumInDelt = (sumInDelt == 0) ? 1 : sumInDelt;
  
  double P_ADR_SEQ = yesInDelt / sumInDelt;
  double P_ADR_NotSEQ = yesNotInDelt / sumNotInDelt;
  
  //same as the sumInDelt
  P_ADR_NotSEQ = P_ADR_NotSEQ == 0 ? 0.00001 : P_ADR_NotSEQ;
  return P_ADR_SEQ / P_ADR_NotSEQ;
}