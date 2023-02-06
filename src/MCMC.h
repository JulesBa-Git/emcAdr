#ifndef MCMC_H
#define MCMC_H


#include "Individual.h"
#include <memory>
#include <cmath>
#include <algorithm>
#include <set>

//' Get the average number of medications taken by a patient (made for the observations set)
//'
//'@param observations : a list of medications for each patient 
double meanMedications(const Rcpp::List& observations);

//' return the larger RR and the corresponding individual
//' 
//' @param a list of 4 RR to compare
std::pair<Individual,double> largerRR(const std::pair<Individual,double>& firstRR,const std::pair<Individual,double>& secRR,
                                      const std::pair<Individual,double>& thirdRR,const std::pair<Individual,double>& fourthRR);
 
 //' add the RR in the right box of the RR distributions array
 //' 
 //' @param an RR
 //' @param the list of RR distribution
void addRRtoDistribution(const double,std::vector<unsigned int>&);

//' Add the pair to the set of explored pair if needed
//' 
//' @param i : the individual that may be added
//' @param the current list of explored pairs
void addPairToSet(const Individual& i, std::set<std::pair<int,int>>&);
 
 //' Check if a Result is already in the result vector or not
 //' 
 //' @param bestResults : the current list of returned result (on the emc algorithm)
 //' @param bestResult : the result which needs to be tested
 //' 
 //' @return true if the result is already in the results vector, false otherwise
bool isNotInResultList(const std::vector<std::pair<Individual,double>>& bestResults,
                       const std::pair<Individual,double>& bestResult);
//'Create the individuals vector with starting individuals only
//'
//'@param startingInd : starting individuals given by the user
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOtemp(const Rcpp::List& startingInd);

//'Create the individuals vector with starting individuals and starting temperatures
//'
//'@param startingInd : starting individuals given by the user
//'@param startingTemp : starting temperatures given by the user
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_Wtemp(const Rcpp::List& startingInd,const Rcpp::NumericVector& startingTemp);

//'Create the individuals vector randomly
//'
//'@param treeSize : size of the ATC tree (to get the DFS index interval)
//'@param nbIndividuals : number of individuals in the population
//'@param meanMedic : the average number of medications took by the observations patients
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals,double meanMedic);

//'Return the Pairs cocktails causing the ADR (used to check which part of the space is explored)
//'
//'@param observationsMed : The List containing the medications taken by the real person
//'@param ADR : A logical vector containing the ADR associated for each patients of the first array (observationsMed)
//'
//'@return the Pairs cocktails causing the ADR 
std::set<std::pair<int,int>> getADRPairs(const Rcpp::List& observationsMed, const Rcpp::LogicalVector& ADR);

//' Return a Mutated version of the individual in parameter (using the 1st mutation)
//' 
//' @param indiv : the individual chosen to be mutated
//' @param treeSize : size of the ATC tree
//' @param alpha : a hyperparameter allowing us to manage to probability of adding a drug to the cocktail. The probability
//' to add a drug to the cocktail is the following : $$ \alpha / n$$ Where n is the original size of the cocktail. 
//' 
//' @return the mutated individual
Individual type1Mutation(const Individual& indiv, int treeSize, double alpha);

//' Return a mutated version of the individual in parameter (using the 2nd mutation)
//'
//'@param indiv : the individual chosen to be mutated
//'@param treeSize : size of the ATC tree
//'@param p : a pair of <int,int> representing the vertex being modified
//'
//'@return the mutated individual
Individual type2Mutation(const Individual& indiv, int treeSize, const std::pair<int,int>& p);

//' Return a Mutated version of the individual in parameter (using the crossover mutation)
//' 
//' @param indiv1 : the individual chosen to be mutated
//' @param indiv2 : the individual with which the indiv1 will be mutated
//' @param ATCtree : tree with every medication
//' @param selectedNode : represent the internal node of the tree on which we will perform the crossover
//' @param upperBound : the upper bound of the set to consider when performing a crossover. Note the interval 
//' to swap between indiv1 and indiv2 is [selectedNode ; upperBound[
//' 
//' @return the mutated individual
Individual crossoverMutation(const Individual& indiv,const Individual& indiv2,const Rcpp::DataFrame& ATCtree,
                             int selectedNode, int upperBound);

#endif