#ifndef MCMC_H
#define MCMC_H


#include "Individual.h"
#include <memory>
#include <cmath>
#include <algorithm>

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
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals);

//' Return a Mutated version of the individual in parameter (using the 1st mutation)
//' 
//' @param indiv : the individual chosen to be mutated
//' @param treeSize : size of the ATC tree
//' 
//' @return the mutated individual
Individual type1Mutation(const Individual& indiv, int treeSize);

//' Return a Mutated version of the individual in parameter (using the crossover mutation)
//' 
//' @param indiv1 : the individual chosen to be mutated
//' @param indiv2 : the individual with which the indiv1 will be mutated
//' @param ATCtree : tree with every medication
//' 
//' @return the mutated individual
Individual crossoverMutation(const Individual& indiv,const Individual& indiv2,const Rcpp::DataFrame& ATCtree);

#endif