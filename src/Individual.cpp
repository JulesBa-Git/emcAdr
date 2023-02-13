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
      if(observation[idx] >= medIdx && observation[idx] <= upperBound[medIdx-1]){
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
  
  //if the cocktail is empty, the related risk is zero.
  if(this->medications_.size() == 0)
    return 0.0;
  
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
  //during the test the denominator was frequently 0 so I make it really small if it is 0 
  //because if the denominator is 0 the numerator has to be 0 so the result would be 0 
  //no matter the denominator, we have 
  sumInDelt = (sumInDelt == 0) ? 1 : sumInDelt;
  
  double P_ADR_SEQ = yesInDelt / sumInDelt;
  double P_ADR_NotSEQ = yesNotInDelt / sumNotInDelt;
  
  //same as the sumInDelt
  P_ADR_NotSEQ = P_ADR_NotSEQ == 0 ? 0.00001 : P_ADR_NotSEQ;

  return P_ADR_SEQ / P_ADR_NotSEQ;
}

std::vector<std::pair<int,int>> Individual::getVertexList(const Rcpp::DataFrame& ATCtree) const{
  std::vector<std::pair<int,int>> returnedVec{0};
  
  std::vector<int> upperBound = ATCtree["upperBound"];
  std::vector<int> depth = ATCtree["ATC_length"];
  
  int idx,depthMed,nextDepth,upperBMed;
  for(const auto& med : medications_){
    idx = med;
    depthMed = depth[med-1]; // minus 1 because per example the medication n° 317 in R numerotation will be at the index 316
    upperBMed = upperBound[med-1];
    //get the next depth
    if(depthMed == 1 || depthMed == 5)
      nextDepth = depthMed+2;
    else if(depthMed < 7)
      nextDepth = depthMed+1;
    
    //find the lower depth medications if we are not on a leaf
    if(depthMed != 7){
      
      while(idx < upperBMed){
        //if we are on the lower depth and the medication is not on the current medications vector
        if(depth[idx] == nextDepth && (std::find(medications_.begin(),medications_.end(),idx+1) == medications_.end())){
          returnedVec.emplace_back(med,idx+1);
        }
        ++idx;
      }
    
    }
    //we come back on the medication (so we got minus 1 because of the R <-> cpp index convertion)
    idx = med-1;
    //now find the upper depth medication (it should only have one because it is a tree), we do it if we are not on the first depth
    if(depthMed != 1){
      //the test >=0 may be deleted because we always should land on a lower depth if our current depth is not 1
      while(idx >= 0 && depthMed <= depth[idx]){
        idx--;
      }
      //we add it if it is not already on the medications vector
      if(std::find(medications_.begin(),medications_.end(),idx+1) == medications_.end())
        returnedVec.emplace_back(med,idx+1);
      
    }
    
  }

  return returnedVec;
}

bool Individual::operator==(const Individual& ind) const{
  if(ind.medications_.size() != medications_.size())
    return false;
  
  for(int i : medications_){
    int j =0;
    while(j < ind.medications_.size() && ind.medications_[j] != i){
      ++j;
    }
    if(j == ind.medications_.size())
      return false;
  }
  
  return true;
}

