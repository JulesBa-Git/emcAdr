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
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observation : real observation of the ADR based on the medications of each real patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'
//'@param nbResults : number of results returned (best RR individuals), 5 by default
//'@param nbIndividuals : number on individuals in the population, 5 by default
//'@param startingIndividuals : starting individuals, randomly initialized by default 
//'(same form as the observations)
//'@param startingTemperatures : starting temperatures, randomly initialized by default
//'@param P_type1/P_type2/P_crossover : probability to operate respectively type1 mutation, type2 mutation and crossover. Note :
//'the probability to operate the swap is then 1 - sum(P_type1,P_type2,P_crossover). The sum must be less or equal to 1. 
//'@param alpha : a hyperparameter allowing us to manage to probability of adding a drug to the cocktail. The probability
//' to add a drug to the cocktail is the following : \eqn{\frac12}{\alpha/n} Where n is the original size of the cocktail. 1 is the default value.
//'
//'@return if no problem return an R List with : the distribution of RR we've met; bests individuals and the corresponding RRs. Otherwise the list is empty
//'@export
//[[Rcpp::export]]
Rcpp::List EMC(int n,const DataFrame& ATCtree,const DataFrame& observations, double P_type1 =.25,
               double P_type2 = .25, double P_crossover =.25, int nbIndividuals = 5, int nbResults = 5,
               double alpha = 1,
               Rcpp::Nullable<Rcpp::List> startingIndividuals = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> startingTemperatures = R_NilValue){
  const int RRDistSize = 42;
  Rcpp::NumericVector temperatures;
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  
  int chosenIndividual_k,chosenIndividual_l;
  double RRx_k,RRx_l, RRy_k,RRy_l, pMutation, pAcceptation, pDraw, q_y_given_x, q_x_given_y;
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  
  //crossover
  int selectedNode, upperBound;
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
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
  
  //acceptance rate
  int acceptedMove = 0;
  int nonZeroRR = 0;
  
  //ADR pairs 
  std::set<std::pair<int,int>> ADRPairs;
  std::set<std::pair<int,int>> exploredPairs;
  std::vector<std::pair<int,int>> interserction; 
  
  double type2Bound;
  double crossoverBound;
  if((P_type1 + P_type2 + P_crossover) > 1){
    std::cerr << "The sum of the 3 probabilities is greater than 1, it must be inferior or equal\n";
    return Rcpp::List::create();
  }else{
    type2Bound = P_type1 + P_type2;
    crossoverBound = type2Bound + P_crossover;
  }
  
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
    
    if(tmpST.size() != nbIndividuals){
      std::cerr << "give the same number of individuals as the value of the parameter nbIndividuals \n";
      return Rcpp::List::create();
    }
    
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
  

  ADRPairs = getADRPairs(observationsMedication, observationsADR);
  
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
    if( pMutation <= P_type1){ // pMutation <= P_type1    
      //type 1 mutation, we take nbIndividuals included because we trunc, there is almost no chance to
      //draw nbIndividuals
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);

      seqLen = individualsCpp[chosenIndividual_k].getMedications().size();
      bool emptySeq = (seqLen == 0);
      mutatedIndividual_k = type1Mutation(individualsCpp[chosenIndividual_k], ATCtree.nrow(), alpha, emptySeq);
      
      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);
      
      //for acceptance rate -> debug
      if(RRy_k > 0)
        ++nonZeroRR;
      
      //to have an overview of the explored space
      addPairToSet(mutatedIndividual_k, exploredPairs);
      
      
      //acceptation probability
      if(mutatedIndividual_k.getMedications().size() > seqLen){
        //if the Y cocktail size is larger than the X (so the mutation 1 added a medication)
        if(seqLen > 0){
          //if the cocktail was not empty
          double U_p = std::min(1.0, (alpha/static_cast<double>(seqLen)));
          q_y_given_x = U_p * ( 1.0 / (ATCtree.nrow() - seqLen) );
          
        }else if(seqLen == 0){
          // if it was empty U_p = 1 and seqLen = 0 
          q_y_given_x = 1.0 / static_cast<double>(ATCtree.nrow());
        }
        double U_p_plus = std::min(1.0, (alpha/static_cast<double>(seqLen + 1)));
        q_x_given_y = ((1 - U_p_plus) * (1.0 / static_cast<double>(seqLen+1)));
      }
      else{
        // if the X cocktail size is larger than the Y (so the mutation 1 removed a medication)
        double U_p = std::min(1.0, (alpha/static_cast<double>(seqLen)));
        q_y_given_x = ((1-U_p) * (1.0/static_cast<double>(seqLen)));
        
        if(seqLen > 1){
          // X cocktail had more than a single medication
          double U_p_minus = std::min(1.0, (alpha/static_cast<double>(seqLen - 1)));
          q_x_given_y = (U_p_minus * (1.0/static_cast<double>(ATCtree.nrow() - seqLen + 1)));
        }else if(seqLen == 1){
          q_x_given_y = (1.0 / static_cast<double>(ATCtree.nrow()));
        }
      }
      
      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature())) 
        * ( q_y_given_x / q_x_given_y);

      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
      bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        ++acceptedMove;
      }
      
    }
    else if(pMutation > P_type1 && pMutation <= type2Bound){//pMutation > 0.25 && pMutation <= 0.50 (default parameter)
      //type 2 mutation

      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;

      RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);

      //get every vertex 0/1 for this patient
      vertexX = individualsCpp[chosenIndividual_k].getVertexList(ATCtree);
      chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
      // TODO : /!\ quand le cocktail est vide, on peut avoir un indice de -1 -> interdire sur cocktail vide.
      chosenVertexidx = chosenVertexidx == vertexX.size() ? vertexX.size()-1 : chosenVertexidx;
      
      std::pair<int,int> chosenVertex = vertexX[chosenVertexidx];
      
      mutatedIndividual_k = type2Mutation(individualsCpp[chosenIndividual_k],ATCtree.nrow(),chosenVertex);

      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);
      vertexY = mutatedIndividual_k.getVertexList(ATCtree);
      
      //for acceptance rate -> debug
      if(RRy_k > 0)
        ++nonZeroRR;
      //to have an overview of the explored space
      addPairToSet(mutatedIndividual_k, exploredPairs);
      
      
      pAcceptation = (exp(((RRy_k - RRx_k) / individualsCpp[chosenIndividual_k].getTemperature()))) * 
                      (vertexY.size()/vertexX.size());
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
      bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        ++acceptedMove;
      }
      
    }
    else if(pMutation > type2Bound && pMutation <= crossoverBound){//pMutation > 0.50 && pMutation <= 0.75

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
      
      //selection of the node on which we will perform the crossover
      //while we are on a leaf 
      do{
        selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
        selectedNode = selectedNode == ATCtree.nrow() ? ATCtree.nrow()-1 : selectedNode;
      } while (ATClength[selectedNode] == 7);
      
      upperBound = upperBounds[selectedNode];

      mutatedIndividual_k = crossoverMutation(individualsCpp[chosenIndividual_k],
                                              individualsCpp[chosenIndividual_l], ATCtree,
                                              selectedNode, upperBound);
      mutatedIndividual_l = crossoverMutation(individualsCpp[chosenIndividual_l],
                                              individualsCpp[chosenIndividual_k], ATCtree,
                                              selectedNode, upperBound);

      RRy_k = mutatedIndividual_k.computeRR(observationsMedication,observationsADR,ATCtree);
      RRy_l = mutatedIndividual_l.computeRR(observationsMedication,observationsADR,ATCtree);
      
      //for acceptance rate -> debug
      if((RRy_k + RRy_l) > 0)
        ++nonZeroRR;
      
      //to have an overview of the explored space
      addPairToSet(mutatedIndividual_k, exploredPairs);
      addPairToSet(mutatedIndividual_l, exploredPairs);
      

      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature()) 
                           + ((RRy_l - RRx_l)/individualsCpp[chosenIndividual_l].getTemperature()));
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult = largerRR({individualsCpp[chosenIndividual_k],RRx_k},{individualsCpp[chosenIndividual_l],RRx_l},
      {mutatedIndividual_k,RRy_k},{mutatedIndividual_l,RRy_l});
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        individualsCpp[chosenIndividual_l] = mutatedIndividual_l;
        ++acceptedMove;
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
      
      //for acceptance rate
      if((RRy_k + RRy_l) > 0)
        ++nonZeroRR;
      
      //to have an overview of the explored space
      addPairToSet(mutatedIndividual_k, exploredPairs);
      addPairToSet(mutatedIndividual_l, exploredPairs);

      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature()) 
                           + ((RRy_l - RRx_l)/individualsCpp[chosenIndividual_l].getTemperature()));
      
      pDraw = Rcpp::runif(1,0,1)[0];
      
     bestResult = largerRR({individualsCpp[chosenIndividual_k],RRx_k},{individualsCpp[chosenIndividual_l],RRx_l},
              {mutatedIndividual_k,RRy_k},{mutatedIndividual_l,RRy_l});
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        individualsCpp[chosenIndividual_l] = mutatedIndividual_l;
        ++acceptedMove;
      }
    }
    //add the best results to the RR array -> do we consider the 0 RR in the distribution ? 
    if(bestResult.second >=1 && bestResult.second <= 5){
      addRRtoDistribution(bestResult.second,RRDistribution);
    }
    else{
      if(bestResult.second < 1&& bestResult.second > 0){
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
  
  double acceptanceRate = static_cast<double>(acceptedMove) / static_cast<double>(n);
  double acceptedNonZeroRR = static_cast<double>(acceptedMove) / static_cast<double>(nonZeroRR);
  std::cout << "acceptance rate : " << acceptanceRate << '\n';
  std::cout << "acceptance rate when a move is on a non zero RR direction : " << acceptedNonZeroRR<< "\n";
  
  std::vector<std::pair<int,int>> inter;
  //check the intersection between the eplored pairs and the pairs of the observations
  std::set_intersection(ADRPairs.begin(), ADRPairs.end(), exploredPairs.begin(), exploredPairs.end(),
                        std::back_inserter(inter));
  std::vector<std::vector<int>> pairs;
  pairs.reserve(inter.size());
  for(const auto& p: inter){
    pairs.push_back({p.first, p.second});
  }
  
  return Rcpp::List::create(Rcpp::Named("bestCockatils") = returnedMed, Rcpp::Named("bestRR") = returnedRR,
                            Rcpp::Named("RRDistribution") = RRDistribution, Rcpp::Named("ADRPairs") = pairs); //Rcpp::Named("ADRPairs") = ADRPairsSingleton
}