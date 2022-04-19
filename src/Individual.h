#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <iostream>
#include "RcppArmadillo.h"
class Individual{
public:
  Individual() = default;
  Individual(const std::vector<int>& medic, double temperature=1);
  
  void printMedications() const;
  void printTemperature() const; 
  
  inline std::vector<int> getMedications() const{
    return medications_;
  }
  inline double getTemperature() const{
    return temperature_;
  }
  
  inline void setMedications(const std::vector<int>& newMed){
    medications_ = newMed;
  }
  inline void setTemperature(const double newTemp){
    temperature_ = newTemp;
  }
  
  bool matches(const std::vector<int>& observation, const std::vector<int>& upperBound) const;
  double computeRR(const Rcpp::List& medications,const Rcpp::LogicalVector& ADR
                     , const Rcpp::DataFrame& ATCtree) const;
private:
  std::vector<int> medications_;
  double temperature_;
};

#endif