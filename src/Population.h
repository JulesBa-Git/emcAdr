#ifndef POPULATION_H
#define POPULATION_H

#include "MCMC.h"
#include <queue>

class Population{
public:
  Population() = default;
  Population(int nbIndividuals);
  Population(int treeSize, int nbIndividuals,double meanMedic);
  Population(const Population& pop);
  
  inline std::vector<std::pair<double,Individual>> getIndividuals() const{
    return individuals_;
  }
  
  inline void setIndividuals(const std::vector<std::pair<double,Individual>>& newPopulation){
    individuals_ = newPopulation;
  }
  
  void addAnIndividualToPopulation(const std::pair<double,Individual>& ind);
  
  void evaluate(const Rcpp::List& medications, const Rcpp::LogicalVector& ADR
                  , const Rcpp::DataFrame& ATCtree);
  void keepElite(int nbElite, Population& matingPool) const;
  
  void tournamentSelection(int tournamentSize, Population& matingPool, int nbDrawing) const;
  
  void crossover(int nbElite, const std::vector<int>& ATClength, const std::vector<int>& upperBounds,
                 const Rcpp::DataFrame& ATCtree, double p_crossover);
  
  void mutate(int nbElite, double p_mutation, const Rcpp::DataFrame& ATCtree);
  
  void clear();
  
  double getMean() const;
  
  int bestIndividual() const;
  
  void printPopulation(std::ostream& ost) const;
  
  void printSummary(int epoch, double populationMean, int populationBestIndex) const;
  
  Population& operator=(const Population& pop);
  
  friend std::ostream& operator<<(std::ostream& ost, const Population& pop);
  
private:
  std::vector<std::pair<double,Individual>> individuals_;
};

#endif