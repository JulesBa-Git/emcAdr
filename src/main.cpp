// we only include MCMC.h which pulls RcppArmadillo.h and Rcpp.h in for us
#include "MCMC.h"
#include <iostream>
using Rcpp::DataFrame;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]


//'The Evolutionary MCMC method that runs the random walk
//'
//'@param n : number of step 
//'@param nbIndividuals : number on individuals in the population 5 by default
//'@param startingIndividuals : starting individuals, randomly initialized by default 
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'@param startingTemperatures : starting temperatures, randomly initialized by default
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observation : real observation of the ADR based on the medications of each real patients
//'(same form as the starting individuals)
//'@export
//[[Rcpp::export]]
void EMC(int n,const DataFrame& ATCtree,const DataFrame& observations, int nbIndividuals = 5,
         Rcpp::Nullable<Rcpp::List> startingIndividuals = R_NilValue,
         Rcpp::Nullable<Rcpp::NumericVector> startingTemperatures = R_NilValue ){
  
  Rcpp::NumericVector temperatures;
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  
  int chosenIndividual_k,chosenIndividual_l;
  double RRx_k,RRx_l, RRy_k,RRy_l, pMutation, pAcceptation, pDraw;
  
  std::vector<Individual> individualsCpp;
  Individual mutatedIndividual_k{};
  Individual mutatedIndividual_l{};
  
  //here if the user enter only temperatures the individuals would be initialized randomly 
  if(startingIndividuals.isNull()){
    //"completely" random initialization
    individualsCpp = DFtoCPP_WOtempAndIndividual(ATCtree.nrow(),nbIndividuals);
  }
  else{ 
    Rcpp::List tmpST = startingIndividuals.clone();
    
    if(startingTemperatures.isNull()){
      //initialized with random temperatures
      individualsCpp = DFtoCPP_WOtemp(tmpST);
    }
    else{
      temperatures = startingTemperatures.clone();
      if(temperatures.length() == tmpST.length()){
        individualsCpp = DFtoCPP_Wtemp(tmpST,temperatures);
        }
      else{
        std::cerr << "there must be the same number of starting point and temperatures" << '\n';
      }
      
    }
  }
  //test affichage individus /*
  for(const Individual& ind : individualsCpp){
    ind.printMedications();
    ind.printTemperature();
  }
  
  
  //for each step suggest a modification using one of the 4 mutations 
  //and compute the probability of acceptation,
  //1st step : compute the RR of the actual population (depending on the chosen mutation)
  //2nd step : apply a mutation 
  //3rd step : compute the RR of the "new population" 
  for(int i = 0; i < n ; ++i){
    pMutation = Rcpp::runif(1,0,1)[0];
    
    
    if(pMutation <= 0.25){     
      //type 1 mutation, we take nbIndividuals included because we trunc, there is almost no chance to
      //draw nbIndividuals
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);
std::cout << "mutation \n";
      mutatedIndividual_k = type1Mutation(individualsCpp[chosenIndividual_k],ATCtree.nrow());
      
      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);
      
      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature()));

      pDraw = Rcpp::runif(1,0,1)[0];
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
      }
      
    }
    else if(pMutation > 0.25 && pMutation <= 0.50){
      //type 2 mutation
      
      
    }
    else if(pMutation > 0.50 && pMutation <= 0.75){
      //crossover mutation
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;
      std::cout << "crossover \n";
      //we chose the second individual he has to be different as the first one so we 
      //do .. while equal to the first one 
      do{
        chosenIndividual_l = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
        chosenIndividual_l = chosenIndividual_l >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_l;
      } while (chosenIndividual_k == chosenIndividual_l);
      
      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication,observationsADR,ATCtree);
      RRx_l = individualsCpp[chosenIndividual_l].computeRR(observationsMedication,observationsADR,ATCtree);
      
      mutatedIndividual_k = crossoverMutation(individualsCpp[chosenIndividual_k],
                                              individualsCpp[chosenIndividual_l],ATCtree);
      mutatedIndividual_l = crossoverMutation(individualsCpp[chosenIndividual_l],
                                              individualsCpp[chosenIndividual_k],ATCtree);
      
      RRy_k = mutatedIndividual_k.computeRR(observationsMedication,observationsADR,ATCtree);
      RRy_l = mutatedIndividual_l.computeRR(observationsMedication,observationsADR,ATCtree);
      
      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature()) 
                           + ((RRy_l - RRx_l)/individualsCpp[chosenIndividual_l].getTemperature()));
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        individualsCpp[chosenIndividual_l] = mutatedIndividual_l;
      }
      
    }
    else{
      // Swap mutation
      
    }
    
  }
  
  for(const Individual& ind : individualsCpp){
    ind.printMedications();
    ind.printTemperature();
  }
}