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
//'@param nbIndividuals : number on individuals in the population, 5 by default
//'@param startingIndividuals : starting individuals, randomly initialized by default 
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'@param startingTemperatures : starting temperatures, randomly initialized by default
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observation : real observation of the ADR based on the medications of each real patients
//'(same form as the starting individuals)
//'@param nbResults : number of results returned (best RR individuals), 5 by default
//'
//'@return if no problem return an R List with : the distribution of RR we've met; bests individuals and the corresponding RRs. Otherwise the list is empty
//'@export
//[[Rcpp::export]]
Rcpp::List EMC(int n,const DataFrame& ATCtree,const DataFrame& observations, int nbIndividuals = 5, int nbResults = 5,
         Rcpp::Nullable<Rcpp::List> startingIndividuals = R_NilValue,
         Rcpp::Nullable<Rcpp::NumericVector> startingTemperatures = R_NilValue ){
  const int RRDistSize = 42;
  Rcpp::NumericVector temperatures;
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  
  int chosenIndividual_k,chosenIndividual_l;
  double RRx_k,RRx_l, RRy_k,RRy_l, pMutation, pAcceptation, pDraw;
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  
  std::vector<Individual> individualsCpp;
  Individual mutatedIndividual_k{};
  Individual mutatedIndividual_l{};
  
  int chosenVertexidx;
  int seqLen;
  
  std::vector<unsigned int> RRDistribution{};
  int smallerOneRR = 0;
  int largerFiveRR = 0;
  RRDistribution.resize(RRDistSize);
  std::vector<std::pair<Individual,double>> bestResults{};
  std::pair<Individual,double> bestResult{};
  bestResults.reserve(nbResults);
  double minRR = 0;
  
  if(nbIndividuals < 2){
    std::cerr << "There must be at least 2 individuals.\n";
    return Rcpp::List::create();
  }
  
  //here if the user enter only temperatures the individuals would be initialized randomly 
  if(startingIndividuals.isNull()){
    //"completely" random initialization using the mean medications per observations
    // we compute the mean minus 1 because we always want a medication so we will use the mean - 1 
    double meanMedicationPerObs = meanMedications(observationsMedication) - 1;
    individualsCpp = DFtoCPP_WOtempAndIndividual(ATCtree.nrow(),nbIndividuals, meanMedicationPerObs);
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
        return Rcpp::List::create();
      }
      
    }
  }
  //test affichage individus /*
  std::cout << "starting cocktails : \n";
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
    if( pMutation <=0.25){ // pMutation <=0.25    
      //type 1 mutation, we take nbIndividuals included because we trunc, there is almost no chance to
      //draw nbIndividuals
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);

      mutatedIndividual_k = type1Mutation(individualsCpp[chosenIndividual_k],ATCtree.nrow());
      seqLen = individualsCpp[chosenIndividual_k].getMedications().size();
      
      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);

      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature())) 
        * (((1/seqLen) * (1/(ATCtree.nrow()-seqLen))) / ((1/seqLen) * (1/(seqLen+1))));

      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
      bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
      }
      
    }
    else if(pMutation > 0.25 && pMutation <= 0.50){//pMutation > 0.25 && pMutation <= 0.50
      //type 2 mutation

      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);

      //get every vertex 0/1 for this patient
      vertexX = individualsCpp[chosenIndividual_k].getVertexList(ATCtree);
      chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
      chosenVertexidx = chosenVertexidx == vertexX.size() ? vertexX.size()-1 : chosenVertexidx;
      
      std::pair<int,int> chosenVertex = vertexX[chosenVertexidx];
      
      mutatedIndividual_k = type2Mutation(individualsCpp[chosenIndividual_k],ATCtree.nrow(),chosenVertex);

      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);
      vertexY = mutatedIndividual_k.getVertexList(ATCtree);

      pAcceptation = (exp(((RRy_k - RRx_k) / individualsCpp[chosenIndividual_k].getTemperature()))) * 
                      (vertexY.size()/vertexX.size());
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
      bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
      }
      
    }
    else if(pMutation > 0.50 && pMutation <= 0.75){//pMutation > 0.50 && pMutation <= 0.75

      //crossover mutation
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

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
      
      bestResult = largerRR({individualsCpp[chosenIndividual_k],RRx_k},{individualsCpp[chosenIndividual_l],RRx_l},
      {mutatedIndividual_k,RRy_k},{mutatedIndividual_l,RRy_l});
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        individualsCpp[chosenIndividual_l] = mutatedIndividual_l;
      }
      
    }
    else{

      // Swap mutation
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals-1)[0]);
      chosenIndividual_k = chosenIndividual_k == nbIndividuals-1 ? nbIndividuals-2 : chosenIndividual_k;
      
      chosenIndividual_l = chosenIndividual_k+1;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication,observationsADR,ATCtree);
      RRx_l = individualsCpp[chosenIndividual_l].computeRR(observationsMedication,observationsADR,ATCtree);

      //we just swap medications
      mutatedIndividual_k = {individualsCpp[chosenIndividual_l].getMedications(),
                             individualsCpp[chosenIndividual_k].getTemperature()};
      mutatedIndividual_l = {individualsCpp[chosenIndividual_k].getMedications(),
                             individualsCpp[chosenIndividual_l].getTemperature()};

      RRy_k = mutatedIndividual_k.computeRR(observationsMedication,observationsADR,ATCtree);
      RRy_l = mutatedIndividual_l.computeRR(observationsMedication,observationsADR,ATCtree);

      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature()) 
                           + ((RRy_l - RRx_l)/individualsCpp[chosenIndividual_l].getTemperature()));
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
     bestResult = largerRR({individualsCpp[chosenIndividual_k],RRx_k},{individualsCpp[chosenIndividual_l],RRx_l},
              {mutatedIndividual_k,RRy_k},{mutatedIndividual_l,RRy_l});
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        individualsCpp[chosenIndividual_l] = mutatedIndividual_l;
      }
    }
    //add the best results to the RR array
    if(bestResult.second >=1 && bestResult.second <= 5){
      addRRtoDistribution(bestResult.second,RRDistribution);
    }
    else{
      if(bestResult.second < 1 && bestResult.second > 0){
        ++smallerOneRR;
      }
      else if(bestResult.second > 5)
        ++largerFiveRR;
    }
    
    //adding the result to the best result if we need to 
    if(isNotInResultList(bestResults,bestResult)){
      if(bestResults.size() < nbResults){
        bestResults.emplace_back(bestResult);
        minRR = minRR < bestResult.second ? minRR : bestResult.second;
      }
      else if(minRR < bestResult.second){
        auto it = std::find_if(bestResults.begin(),bestResults.end(),
                     [minRR](const std::pair<Individual,double>& p){return p.second == minRR;});
        if(it != bestResults.end()){
          bestResults.erase(it);
        }
        bestResults.emplace_back(bestResult);
        auto tmpMin = *std::min_element(bestResults.begin(),bestResults.end(),
                                 [](const std::pair<Individual,double>& lp,const std::pair<Individual,double>& rp){
                                   return lp.second < rp.second; 
                                 });
        minRR = tmpMin.second;
      }
    }
  }
  
  //insert le smaller one and the larger five to the RR vector
  RRDistribution.insert(RRDistribution.begin(),smallerOneRR);
  RRDistribution[RRDistribution.size()-1] = largerFiveRR;
  
  //create the returned vector
  std::vector<std::vector<int>> returnedMed{};
  returnedMed.reserve(bestResults.size());
  std::vector<double>returnedRR{};
  returnedRR.reserve(bestResults.size());

  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returnedRR.push_back(pair.second);
  }
  
  return Rcpp::List::create(Rcpp::Named("bestCockatils") = returnedMed,Rcpp::Named("bestRR") = returnedRR,
                            Rcpp::Named("RRDistribution") = RRDistribution);
}