// we only include MCMC.h which pulls RcppArmadillo.h and Rcpp.h in for us
#include "MCMC.h"
#include "Population.h"
#include <iostream>
#include <string_view>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <set>

#ifdef _OPENMP
  #include <omp.h>
#endif

using Rcpp::DataFrame;
// [[Rcpp::plugins(openmp)]]
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
  double q_ratio;

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
  //best individual in the population each round
  std::vector<double> populationRR{};
  std::vector<unsigned int> RRDistributionPopulation{};
  populationRR.resize(nbIndividuals);
  RRDistributionPopulation.resize(RRDistSize);
  
  //acceptance rate
  int acceptedMove = 0;
  int nonZeroRR = 0;
  int acceptedSize = 0;
  int accepted_type1 =0;
  int accepted_type2 =0;
  int accepted_cross =0;
  int accepted_swap =0;
  int type1_move=0, type2_move=0, cross_move=0, swap_move=0;
  
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
    if(!startingTemperatures.isNull()){
      temperatures = startingTemperatures.clone();

      if(temperatures.size() != nbIndividuals){
        std::cerr << "there must be nbIndividuals temperatures given, here there is : "<< temperatures.size() <<"\n";
        return Rcpp::List::create();
      }
      else
        individualsCpp = DFtoCPP_WOIndividual(ATCtree.nrow(),nbIndividuals, meanMedicationPerObs, temperatures);
    }
    else
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
  
  //compute each RR before the loop
  for(int i = 0; i < individualsCpp.size(); ++i){
    populationRR[i] = individualsCpp[i].computeRR(observationsMedication, observationsADR, ATCtree);
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

      //RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);
      RRx_k = populationRR[chosenIndividual_k];
      
      seqLen = individualsCpp[chosenIndividual_k].getMedications().size();
      bool emptySeq = (seqLen == 0);
      
      //mutatedIndividual_k = type1Mutation(individualsCpp[chosenIndividual_k], ATCtree.nrow(), alpha, emptySeq); test new mutation
      mutatedIndividual_k = adjustedType1Mutation(individualsCpp[chosenIndividual_k], ATCtree.nrow(), alpha, emptySeq);
      
      RRy_k = mutatedIndividual_k.computeRR(observationsMedication, observationsADR, ATCtree);
      
      //for acceptance rate -> debug
      if(RRy_k > 0)
        ++nonZeroRR;
      //to have an overview of the explored space
      addPairToSet(mutatedIndividual_k, exploredPairs);
      
      
      //acceptation probability
      /*if(mutatedIndividual_k.getMedications().size() > seqLen){
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
        q_x_given_y = ((1 - U_p_plus) * (1.0 / static_cast<double>(seqLen + 1)));
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
        * ( q_x_given_y / q_y_given_x); //q_y_given_x / q_x_given_y*/
      //new acceptation probability (for the adjusted type1 mutation)! 
      if(mutatedIndividual_k.getMedications().size() > seqLen){
        //we added a medication
        if(emptySeq){
          q_ratio = 1-alpha;
        }
        else{
          q_ratio = ((1-alpha) * (1.0 / static_cast<double>(seqLen))) / static_cast<double>(alpha);
        }
      }
      else if(mutatedIndividual_k.getMedications().size() < seqLen){
        //if the modification led to a deletion
        if(seqLen == 1){
          q_ratio = 1.0 / (1-alpha);
        }
        else{
          q_ratio = alpha / ((1-alpha) * (1.0/static_cast<double>(seqLen)));
        }
      }else{
        //here we just modified the cocktail (no addition no deletion)
        q_ratio = (1-alpha) / alpha;
      }
      pAcceptation = exp(((RRy_k - RRx_k)/individualsCpp[chosenIndividual_k].getTemperature())) 
        * q_ratio;

      pDraw = Rcpp::runif(1,0,1)[0];
      
      bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
      bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
      
      if(pAcceptation > pDraw){
        individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
        populationRR[chosenIndividual_k] = RRy_k;
        ++acceptedMove;
        ++accepted_type1;
        acceptedSize += individualsCpp[chosenIndividual_k].getMedications().size();
      };
      ++type1_move;
    }
    else if(pMutation > P_type1 && pMutation <= type2Bound){//pMutation > 0.25 && pMutation <= 0.50 (default parameter)
      //type 2 mutation

      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals)[0]);
      chosenIndividual_k = chosenIndividual_k >= nbIndividuals ? nbIndividuals-1 : chosenIndividual_k;
      
      //if the selected indivudual is empty, the type 2 mutation cannot occur
      if(individualsCpp[chosenIndividual_k].getMedications().size() != 0){
        // RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication, observationsADR, ATCtree);
        RRx_k = populationRR[chosenIndividual_k];
        //get every vertex 0/1 for this patient
        vertexX = individualsCpp[chosenIndividual_k].getVertexList(ATCtree);
        chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
        
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
          (static_cast<double>(vertexX.size())/static_cast<double>(vertexY.size()));
        
        pDraw = Rcpp::runif(1,0,1)[0];
        
        bestResult.second = RRx_k > RRy_k ? RRx_k : RRy_k;
        bestResult.first = RRx_k > RRy_k ? individualsCpp[chosenIndividual_k] : mutatedIndividual_k;
        
        if(pAcceptation > pDraw){
          individualsCpp[chosenIndividual_k] = mutatedIndividual_k;
          populationRR[chosenIndividual_k] = RRy_k;
          ++acceptedMove;
          ++accepted_type2;
          acceptedSize += individualsCpp[chosenIndividual_k].getMedications().size();
        }
        ++type2_move;
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

      //RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication,observationsADR,ATCtree);
      //RRx_l = individualsCpp[chosenIndividual_l].computeRR(observationsMedication,observationsADR,ATCtree);
      RRx_k = populationRR[chosenIndividual_k];
      RRx_l = populationRR[chosenIndividual_l];
      
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
        populationRR[chosenIndividual_k] = RRy_k;
        populationRR[chosenIndividual_l] = RRy_l;
        ++acceptedMove;
        ++accepted_cross;
        acceptedSize += individualsCpp[chosenIndividual_k].getMedications().size();
      }
      ++cross_move;
    }
    else{

      // Swap mutation
      chosenIndividual_k = trunc(Rcpp::runif(1,0,nbIndividuals-1)[0]);
      chosenIndividual_k = chosenIndividual_k == nbIndividuals-1 ? nbIndividuals-2 : chosenIndividual_k;
      
      chosenIndividual_l = chosenIndividual_k+1;

      //RRx_k = individualsCpp[chosenIndividual_k].computeRR(observationsMedication,observationsADR,ATCtree);
      //RRx_l = individualsCpp[chosenIndividual_l].computeRR(observationsMedication,observationsADR,ATCtree);
      RRx_k = populationRR[chosenIndividual_k];
      RRx_l = populationRR[chosenIndividual_l];
      
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
        populationRR[chosenIndividual_k] = RRy_k;
        populationRR[chosenIndividual_l] = RRy_l;
        ++acceptedMove;
        ++accepted_swap;
        acceptedSize += individualsCpp[chosenIndividual_k].getMedications().size();
      }
      ++swap_move;
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
    
    /*if(pMutation < P_type1){
      std::cout << "epochs " << i <<'\n';
      std::cout << "proba mutation " << pMutation << "\n";
      std::cout << "proba d'acceptation " << pAcceptation << '\n';
      std::cout << "sequence mutée "<< chosenIndividual_k << '\n';
      std::cout << "RRx avant mutation " << RRx_k << " RRy apres mutation " << RRy_k << '\n';
      //std::cout << "q(y | x) " << q_y_given_x << " q(x|y) " << q_x_given_y << " rapport " << q_x_given_y / q_y_given_x << "\n";
      std::cout << "q ratio : " << q_ratio << '\n'; 
      for(const auto& ind : individualsCpp){
        std::cout << ind.getMedications().size() << '\n';
      }
      
    }*/
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
  double meanCocktailsSize = static_cast<double>(acceptedSize) / static_cast<double>(acceptedMove);
  std::cout << "global acceptance rate : " << acceptanceRate << '\n';
  std::cout << "type 1 mutation acceptance rate : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move) << "\n";
  std::cout << "type 2 mutation acceptance rate : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move) << "\n";
  std::cout << "crossover acceptance rate : " << static_cast<double>(accepted_cross) / static_cast<double>(cross_move) << "\n";
  std::cout << "swap mutation acceptance rate : " << static_cast<double>(accepted_swap) / static_cast<double>(swap_move) << "\n";
  
  std::cout << "acceptance rate when a move is on a non zero RR direction : " << acceptedNonZeroRR<< "\n";
  std::cout << "mean size of accepted cocktails : " << meanCocktailsSize << "\n";
  
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

//'The Evolutionary MCMC method that runs the random walk on a single cocktail
//'
//'@param epochs : number of step 
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observation : real observation of the ADR based on the medications of each real patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'
//'@param temperature : starting temperature, default = 1
//'@param nbResults : Number of returned solution (Cocktail of size Smax), 5 by default
//'@param Smax : Size of the cocktail we approximate the distribution from
//'@param p_type1: probability to operate type1 mutation. Note :
//'the probability to operate the type 2 mutation is then 1 - P_type1. P_type1 must be in [0;1]. 
//'@param beta : filter the minimum number of patients that must have taken the 
//'cocktail for it to be considered 'significant', default is 4
//'@param Upper bound of the considered Metric. 
//'@param num_thread : Number of thread to run in parallel
//'
//'@return if no problem return an array of the approximation of the RR distribution : the distribution of RR we've met; Otherwise the list is empty
//'@export
//[[Rcpp::export]]
Rcpp::List DistributionApproximation(int epochs, const DataFrame& ATCtree, const DataFrame& observations,
                                              int temperature = 1, int nbResults = 5, int Smax = 2,
                                              double p_type1 = .01, int beta = 4, int max_Metric = 100,
                                              int num_thread = 1){
  //arguments verification
  if(p_type1 > 1 || p_type1 < 0 || epochs < 1){
    std::cerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::List();
  }
  
  // OMP SET NUM THREAD = k, s.t. 1 <= k <= omp_get_num_procs()
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    std::cerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  double p_type2 = 1 - p_type1;
  Individual cocktail, mutatedIndividual;
  
  //for the moment the distribution is bounded by 0 and RRmax
  //const int distribSize = 300;
  std::vector<unsigned int> RRDistribution((max_Metric*10) +1); // +1 stand for every RR over the RRmax value
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstandingRR{};
  outstandingRR.reserve(10);
  std::pair<double, std::pair<int, int>> currentRR = std::make_pair(0.0,std::make_pair(0,0));
  std::pair<double, std::pair<int, int>> computeRROutput; // pair< RR, <N° of people taking the cocktail and having the ADR, N° of people taking the cocktail>>
  
  std::vector<unsigned int> RRDistributionGreaterBeta((max_Metric*10) +1);
  
  //used in the phyper function
  int ADRCount = 0;
  for(const auto& adr : observationsADR){
    if(adr)
      ++ADRCount;
  }
  int notADRCount = observationsMedication.size() - ADRCount;
  double ADRProp = static_cast<double>(ADRCount) / static_cast<double>(observationsMedication.size());
    
  std::pair<double, std::pair<int, int>> currentGeom = std::make_pair(0.0,std::make_pair(0,0));
  std::pair<double, std::pair<int, int>> computeGeomOutput; // pair< phypergeometric, <N° of people taking the cocktail and having the ADR, N° of people taking the cocktail>>
  double minGeom = 0;
  double minGeomBeta = 0;
  /*
  std::pair<double, std::pair<int, int>> currentBinom = std::make_pair(0.0,std::make_pair(0,0));
  std::pair<double, std::pair<int, int>> computePBinomOutput; // pair< phypergeometric, <N° of people taking the cocktail and having the ADR, N° of people taking the cocktail>>
  double minBinom = 0;
  double minBinomBeta = 0;*/
  

  //acceptance rate
  int acceptedMove = 0;
  int accepted_type1 =0;
  int accepted_type2 =0;
  int type1_move=0, type2_move=0;
  int type1_move_inF = 0, type2_move_inF = 0;
  int falseAcceptedCocktailCount = 0, falseSampledCocktailCount= 0;
  
  double RRx_k, RRy_k, q_y_given_x, q_x_given_y, pMutation, pAcceptation, pDraw;
  double q_ratio;
  int seqLen;
  
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  int chosenVertexidx;
  
  std::vector<std::pair<Individual,double>> bestResults;
  bestResults.reserve(nbResults);
  std::vector<std::pair<Individual,double>> bestResultsBeta;
  bestResultsBeta.reserve(nbResults);
  
  double minRR = 0;
  double minRRbeta = 0;
  std::pair<Individual, double> currentResult;
  
  //if p.second is greater than 0, it means that the cocktail correspond to at least one person
  auto belongToF = [](const std::pair<double,std::pair<int,int>>& p){
    return (p.second.second > 0);
  };

  //initialization (every cocktail has to contain Smax medication)
  do{
    cocktail = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
    //computeRROutput = cocktail.computeRR(observationsMedication, observationsADR, upperBounds, RRmax, num_thread);
    computeGeomOutput = cocktail.computePHypergeom(observationsMedication, observationsADR,
                                                   upperBounds, ADRCount, notADRCount,
                                                   max_Metric, num_thread);
    /*computePBinomOutput = cocktail.computePBinomial(observationsMedication, observationsADR,
                                                    upperBounds, ADRProp,
                                                    max_Metric, num_thread);*/
  } while (!belongToF(computeGeomOutput));
  //currentRR = computeRROutput;
  currentGeom = computeGeomOutput;
  
#ifdef _OPENMP
  std::cout<< "openMP available \n";
#endif
  
  //minRR = currentRR.first;
  minGeom = currentGeom.first;
  bestResults.emplace_back(cocktail, currentGeom.first);
  if(currentGeom.second.second > beta){
    minGeomBeta = currentGeom.first;
    bestResultsBeta.emplace_back(cocktail, currentGeom.first);
  }
  
  for(int i = 0; i < epochs; ++i){
      pMutation = Rcpp::runif(1,0,1)[0];
      
      if(pMutation < p_type1){
        //type 1 mutation
        RRx_k = currentGeom.first;
        
        //mutatedIndividual = type1Mutation(cocktail, ATCtree.nrow(), alpha, emptySeq);
        //mutatedIndividual = adjustedType1Mutation(cocktail, ATCtree.nrow(), alpha, emptySeq);
        //here the current type 1 mutation consist in drawing a new cocktail of the same size
        mutatedIndividual = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
        //computeRROutput = mutatedIndividual.computeRR(observationsMedication, observationsADR, upperBounds, max_Metric, num_thread);
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_Metric,
                                                                num_thread);
        //std::cout << "size : " << mutatedIndividual.getMedications().size() << '\n';
        RRy_k = computeGeomOutput.first;
        //std::cout << "prev RR : " << RRx_k << " new RR : " << RRy_k << '\n';

        //to have an overview of the explored space (not in this method for the moment)
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        

        
        if(belongToF(computeGeomOutput)){
          // with this mutation, our ration q(X|Y) / q(Y|X) = 1
          pAcceptation = exp(((RRy_k - RRx_k)/static_cast<double>(cocktail.getTemperature()))); 
          //std::cout << "proba d'acceptation : " << pAcceptation << '\n';
          pDraw = Rcpp::runif(1,0,1)[0];
          ++type1_move_inF;
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentGeom = computeGeomOutput;
            ++acceptedMove;
            ++accepted_type1;
          }
        }else{
          ++nbCocktailNotInPopulation;
        }
        
        ++type1_move;
      }
      else{
        //type 2 mutation
        //if the selected indivudual is empty, the type 2 mutation cannot occur

        RRx_k = currentGeom.first;
        //get every vertex 0/1 for this patient
        vertexX = cocktail.getVertexList(ATCtree);
        
        chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
        chosenVertexidx = chosenVertexidx == vertexX.size() ? vertexX.size()-1 : chosenVertexidx;
        
        std::pair<int,int> chosenVertex = vertexX[chosenVertexidx];
        
        mutatedIndividual = type2Mutation(cocktail, ATCtree.nrow(), chosenVertex);
        
        //computeRROutput = mutatedIndividual.computeRR(observationsMedication, observationsADR, upperBounds, max_Metric, num_thread);
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_Metric, 
                                                                num_thread);
        RRy_k = computeGeomOutput.first;
        vertexY = mutatedIndividual.getVertexList(ATCtree);
        
        //to have an overview of the explored space
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        
        if(belongToF(computeGeomOutput)){
          pAcceptation = (exp(((RRy_k - RRx_k) / cocktail.getTemperature()))) * 
            (static_cast<double>(vertexX.size())/static_cast<double>(vertexY.size()));
          
          pDraw = Rcpp::runif(1,0,1)[0];
          ++type2_move_inF;
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentGeom = computeGeomOutput;
            ++acceptedMove;
            ++accepted_type2;
          }
        }else{
          ++nbCocktailNotInPopulation;
        }
        ++type2_move;
      
    }
    if(currentGeom.first < max_Metric){
      int index = 10 * currentGeom.first;
      ++RRDistribution[index];
      // in every cases we add the RR to the "normal returned" distribution, and if
      // more than beta persons take it, we add the RR to the other ditribution named
      // RRDistributionGreaterBeta
      if(currentGeom.second.second > beta){ // second.first = N° of people taking the cocktail and having the ADR
        ++RRDistributionGreaterBeta[index];
      }
    }
    else{ // since we have an RR max, we just increment the last elements of the distribution
      // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
      outstandingRR.push_back(currentGeom.first); // could remove this line ?
      ++RRDistribution[RRDistribution.size()-1];
    }
    
    if(!cocktail.isTrueCocktail(upperBounds))
      falseAcceptedCocktailCount++;
    if(!mutatedIndividual.isTrueCocktail(upperBounds))
      falseSampledCocktailCount++;
        
    currentResult = std::make_pair(cocktail, currentGeom.first);
    //adding the result to the best result if we need to 
    minGeom = addToBestCocktails(bestResults, currentResult, nbResults, minGeom,
                               upperBounds);
    
    if(currentGeom.second.second > beta){
      minGeomBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                       minGeomBeta, upperBounds);
    }
  }
  
  //create the returned vector
  std::vector<std::vector<int>> returnedMed{};
  returnedMed.reserve(bestResults.size());
  std::vector<double>returnedRR{};
  returnedRR.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returnedRR.push_back(pair.second);
  }
  
  //create the returned vector with the cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returnedRRBeta{};
  returnedRRBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returnedRRBeta.push_back(pair.second);
  }
  
  
  std::cout << "acceptance rate : " << static_cast<double>(acceptedMove) / static_cast<double>(epochs)<< "\n";
  std::cout << "acceptance rate type1 mutation : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move)<< "\n";
  std::cout << "acceptance rate type2 mutation : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move)<< "\n";
  std::cout << "acceptance rate type1 mutation when the proposal is in F : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move_inF)<< "\n";
  std::cout << "acceptance rate type2 mutation when the proposal is in F : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move_inF)<< "\n";
  std::cout << "number of proposed cocktail that was taken by nobody in the population : " << nbCocktailNotInPopulation << '\n';
  std::cout << "number of false cocktail sampled : " << falseSampledCocktailCount << '\n';
  std::cout << "number of false cocktail concidered in the distribution during the run : " << falseAcceptedCocktailCount << '\n';
  
  return Rcpp::List::create(Rcpp::Named("Distribution") = RRDistribution, Rcpp::Named("OutstandingRR") = outstandingRR, 
                            Rcpp::Named("bestCockatils") = returnedMed, Rcpp::Named("bestRR") = returnedRR,
                            Rcpp::Named("FilteredDistribution") = RRDistributionGreaterBeta,
                            Rcpp::Named("bestCocktailsBeta") = returnedMedBeta,
                            Rcpp::Named("bestRRBeta") = returnedRRBeta,
                            Rcpp::Named("cocktailSize") = Smax);
}

//'Genetic algorithm, trying to reach the best cocktail (the one which maximize
//'the fitness function, Related Risk in our case)
//'
//'@param epochs : number of step 
//'@param nbIndividuals : size of the popultation
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observation : real observation of the ADR based on the medications of each real patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'@param diversity : enable the diversity mechanism of the algorithm
//' (favor the diversity of cocktail in the population)
//'@param p_crossover: probability to operate a crossover on the crossover phase.
//'@param p_mutation: probability to operate a mutation after the crossover phase.
//'@param nbElite : number of best individual we keep from generation to generation
//'@param tournamentSize : size of the tournament (select the best individual be tween tournamentSize individuals) 
//'
//'@return if no problem return the best cocktail we found (according to the fitness function which is the Relative Risk)
//'@export
//[[Rcpp::export]]
Rcpp::List GeneticAlgorithm(int epochs, int nbIndividuals, const DataFrame& ATCtree, 
                            const DataFrame& observations, int num_thread = 1, 
                            bool diversity = false, double p_crossover = .80,
                            double p_mutation = .01, int nbElite = 0, 
                            int tournamentSize = 2, double alpha = 1,
                            bool summary = true){ 
  //arguments verification
  if(p_crossover > 1 || p_crossover < 0 || nbIndividuals < 1 || p_mutation > 1 || p_mutation < 0 || epochs < 1){
    std::cerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::List();
  }
  if(tournamentSize > nbIndividuals || nbElite > nbIndividuals){
    std::cerr << "the tournament size and the nbElite parameters must be less equal than the number of individuals \n";
    return Rcpp::List();
  }
  // OMP SET NUM THREAD = k, s.t. 1 <= k <= omp_get_num_procs()
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    std::cerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  
  
  //since there is a diversity mechanism, we may not want to set an upper bound
  // to the metric (here Phypergeometric), so we put it equal to INT_MAX
  // we may want to put this as a parameter
  int max_metric = INT_MAX;
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
    
  int ADRCount = 0;
  for(const auto& adr : observationsADR){
    if(adr)
      ++ADRCount;
  }
  int notADRCount = observationsMedication.size() - ADRCount;
  
  //generate the initial population randomly (do we consider an Smax ?)
  double meanMedicationPerObs = meanMedications(observationsMedication) - 1;
  Population population(ATCtree.nrow(), nbIndividuals, meanMedicationPerObs,
                        upperBounds);
  Population matingPool(nbIndividuals);
  
  std::vector<double> meanScore;
  meanScore.reserve(epochs);
  std::vector<double> bestScore;
  bestScore.reserve(epochs);
  int bestIndividualIndex;
  
  int remainingDraw = 0;
  
  std::vector<double> score_before_penalization;
  score_before_penalization.reserve(population.getIndividuals().size());
  
  //here we may want to have a more sophisticated stopping condition (like, if the RR is 
  //significantly high given the previous calculated distribution)
  for(int i =0; i < epochs; ++i){
    
    //do we apply the diversity mechanism ?
    if(diversity){
      population.penalize(depth,father);
    }
    
    //1st : fit every individual
    population.evaluate(observationsMedication, observationsADR, upperBounds,
                        ADRCount, notADRCount, max_metric, num_thread);
    
    for(const auto& ind : population.getIndividuals()){
      score_before_penalization.push_back(ind.first);
    }
    
    
    //do we apply the diversity mechanism ?
    if(diversity){
      population.penalize(depth,father);
    }
    
    //2nd : make a ss and apply the crossover on the selected individual
    //keep the elite first
    population.keepElite(nbElite, matingPool);

    //select the individual according to the fitness
    remainingDraw = nbIndividuals - nbElite;
    population.tournamentSelection(tournamentSize, matingPool, remainingDraw);

    //operate a crossover over the mating pool 
    matingPool.crossover(nbElite, ATClength, upperBounds, ATCtree, p_crossover);
    
    //operate a mutation over the mating pool
    matingPool.mutate(nbElite, p_mutation, ATCtree, upperBounds, alpha);
    

    //3rd : replace the population
    population = matingPool;
    matingPool.clear();
    
    meanScore.push_back(population.getMean());
    bestIndividualIndex = population.bestIndividual();
    bestScore.push_back(population.getIndividuals()[bestIndividualIndex].first);
    
    if(summary)
      population.printSummary(i, meanScore[i], bestIndividualIndex);
    
    score_before_penalization.clear();
    
  }
  population.evaluate(observationsMedication, observationsADR, upperBounds,
                      ADRCount, notADRCount, max_metric, num_thread);

  //output the population 
  std::vector<std::vector<int>> medications;
  medications.reserve(nbIndividuals);
  std::vector<double> populationScore;
  populationScore.reserve(nbIndividuals);
  
  medications = population.getMedications();
  populationScore = population.getRR();

  
  Rcpp::List returnedPop = Rcpp::List::create(Rcpp::Named("cocktails") = medications,
                                              Rcpp::Named("score") = populationScore);
  
  return Rcpp::List::create(Rcpp::Named("meanFitnesses") = meanScore,
                            Rcpp::Named("BestFitnesses") = bestScore,
                            Rcpp::Named("FinalPopulation") = returnedPop);
}


//'The true distribution of drugs
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param obervations : population of patients on which we want to compute the risk distribution
//'@param beta : minimum number of cocktail takers 
//'@param max_risk : maximum risk, at which point the risk is considered exceptional (outliers)
//'
//'@return the risk distribution among drugs
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionDrugs(const DataFrame& ATCtree, const DataFrame& observations,
                                 int beta, int max_risk = 100, int num_thread = 1){
 
#ifdef _OPENMP
 if(num_thread < 1 || num_thread > omp_get_num_procs()){
   std::cerr << "Wrong thread number, it should be between 1 and " 
             << omp_get_num_procs() << " \n";
   return Rcpp::List();
 }
#endif
 
 Rcpp::List observationsMedicationTmp = observations["patientATC"];
 std::vector<std::vector<int>> observationsMedication;
 observationsMedication.reserve(observationsMedicationTmp.size());
 Rcpp::LogicalVector observationsADR = observations["patientADR"];
 std::vector<int> upperBounds = ATCtree["upperBound"];
 int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
 
 for(int i =0; i < observationsMedicationTmp.size(); ++i){
   observationsMedication.push_back(observationsMedicationTmp[i]);
 }
 
 int notADRCount = observationsMedication.size() - ADRCount;
 
 
 //for the moment the distribution is bounded by 0 and 100
 const int distribSize = max_risk *10;
 std::vector<unsigned int> RRDistribution(distribSize);
 std::vector<unsigned int> RRDistributionGreaterBeta(distribSize);
 
 unsigned int nbCocktailNotInPopulation = 0;
 std::vector<double> outstandingRR{};
 outstandingRR.reserve(10);
 std::vector<double> outstandingRRBeta{};
 outstandingRRBeta.reserve(10);
 
 std::pair<double, std::pair<int,int>> computeRROutput;
 std::pair<double, std::pair<int,int>> computePhyperOutput;
 
 int nbResults = 100; // we keep the 100 best individuals
 double minRR = 0, minRRbeta = 0;
 double minPhyper = 0, minPhyperBeta = 0;
 std::vector<std::pair<Individual,double>> bestResults;
 bestResults.reserve(nbResults);
 std::vector<std::pair<Individual,double>> bestResultsBeta;
 bestResultsBeta.reserve(nbResults);
 std::pair<Individual, double> currentResult;
 
 Individual indiv{};
 indiv.setTemperature(1);
 int med;
 for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
   med = i;

   indiv.setMedications({med});
   //computeRROutput = indiv.computeRR(observationsMedication, observationsADR, upperBounds, 8000, num_thread);
   computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                 upperBounds, ADRCount, notADRCount,
                                                 8000, num_thread);
   if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
     if(computePhyperOutput.first < max_risk){
       int index = 10 * computePhyperOutput.first;
       ++RRDistribution[index];
       if(computePhyperOutput.second.second > beta){
         ++RRDistributionGreaterBeta[index];
       }
     }
     else{
       if(computePhyperOutput.second.second > beta){
         outstandingRRBeta.push_back(computePhyperOutput.first);
       }
       outstandingRR.push_back(computePhyperOutput.first);
     }
     
     //keep the best results
     currentResult = std::make_pair(indiv, computePhyperOutput.first);
     minPhyper = addToBestCocktails(bestResults, currentResult, nbResults, minPhyper,
                                    upperBounds);
     
     if(computePhyperOutput.second.second > beta){
       minPhyperBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                          minPhyperBeta, upperBounds);
     }
     
   
   }
 }
 
 //create the vector of the best cocktail
 std::vector<std::vector<int>> returnedMed{};
 returnedMed.reserve(bestResults.size());
 std::vector<double>returnedRR{};
 returnedRR.reserve(bestResults.size());
 
 for(const auto &pair : bestResults){
   returnedMed.push_back(pair.first.getMedications());
   returnedRR.push_back(pair.second);
 }
 
 //create the returned vector of the best cocktail taken by more than beta person
 std::vector<std::vector<int>> returnedMedBeta{};
 returnedMedBeta.reserve(bestResultsBeta.size());
 std::vector<double>returnedRRBeta{};
 returnedRRBeta.reserve(bestResultsBeta.size());
 
 for(const auto &pair : bestResultsBeta){
   returnedMedBeta.push_back(pair.first.getMedications());
   returnedRRBeta.push_back(pair.second);
 }
 
 return Rcpp::List::create(Rcpp::Named("Distribution") = RRDistribution,
                           Rcpp::Named("Filtered_distribution") = RRDistributionGreaterBeta,
                           Rcpp::Named("Outstanding_risk") = outstandingRR,
                           Rcpp::Named("Best_cocktails") = returnedMed,
                           Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                           Rcpp::Named("Best_risks") = returnedRR,
                           Rcpp::Named("Best_risks_beta") = returnedRRBeta);
}

//'The true distribution of size 2 cocktails
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param obervations : population of patients on which we want to compute the risk distribution
//'@param beta : minimum number of cocktail takers 
//'@param max_risk : maximum risk, at which point the risk is considered exceptional (outliers)
//'
//'@return the risk distribution among size 2 cocktails
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations,
                                           int beta, int max_risk = 100, int num_thread = 1){
  
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    std::cerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
  
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  int notADRCount = observationsMedication.size() - ADRCount;

  
  //for the moment the distribution is bounded by 0 and 30
  const int distribSize = max_risk *10;
  std::vector<unsigned int> RRDistribution(distribSize);
  std::vector<unsigned int> RRDistributionGreaterBeta(distribSize);
  
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstandingRR{};
  outstandingRR.reserve(10);
  std::vector<double> outstandingRRBeta{};
  outstandingRRBeta.reserve(10);
  
  double currentRR = 0;
  std::pair<double, std::pair<int,int>> computeRROutput;
  std::pair<double, std::pair<int,int>> computePhyperOutput;
  
  int nbResults = 100; // we keep the 100 best individuals
  double minRR = 0, minRRbeta = 0;
  double minPhyper = 0, minPhyperBeta = 0;
  std::vector<std::pair<Individual,double>> bestResults;
  bestResults.reserve(nbResults);
  std::vector<std::pair<Individual,double>> bestResultsBeta;
  bestResultsBeta.reserve(nbResults);
  std::pair<Individual, double> currentResult;

  Individual indiv{};
  indiv.setTemperature(1);
  std::vector<int> med;
  med.resize(2);
  for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
    med[0] = i;
    for(int j = i+1 ; j < ATCtree.nrow(); ++j){
      med[1] = j;
      indiv.setMedications(med);
      //computeRROutput = indiv.computeRR(observationsMedication, observationsADR, upperBounds, 8000, num_thread);
      computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                    upperBounds, ADRCount, notADRCount,
                                                    8000, num_thread);
      if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
        if(computePhyperOutput.first < max_risk){
          int index = 10 * computePhyperOutput.first;
          ++RRDistribution[index];
          if(computePhyperOutput.second.second > beta){
            ++RRDistributionGreaterBeta[index];
          }
        }
        else{
          if(computePhyperOutput.second.second > beta){
            outstandingRRBeta.push_back(computePhyperOutput.first);
          }
          // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
          outstandingRR.push_back(computePhyperOutput.first);
        }
        
        //keep the best results
        currentResult = std::make_pair(indiv, computePhyperOutput.first);
        minPhyper = addToBestCocktails(bestResults, currentResult, nbResults, minPhyper,
                                   upperBounds);
        
        if(computePhyperOutput.second.second > beta){
          minPhyperBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                             minPhyperBeta, upperBounds);
        }
        
      }
    }
  }
  
  //create the vector of the best cocktail
  std::vector<std::vector<int>> returnedMed{};
  returnedMed.reserve(bestResults.size());
  std::vector<double>returnedRR{};
  returnedRR.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returnedRR.push_back(pair.second);
  }
  
  //create the returned vector of the best cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returnedRRBeta{};
  returnedRRBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returnedRRBeta.push_back(pair.second);
  }
  
  return Rcpp::List::create(Rcpp::Named("Distribution") = RRDistribution,
                            Rcpp::Named("FilteredDistribution") = RRDistributionGreaterBeta,
                            Rcpp::Named("OutstandingRR") = outstandingRR,
                            Rcpp::Named("BestCocktails") = returnedMed,
                            Rcpp::Named("BestCocktailsBeta") = returnedMedBeta,
                            Rcpp::Named("BestRR") = returnedRR,
                            Rcpp::Named("BestRRBeta") = returnedRRBeta);
}

//'The true RR distribution of cocktail of size 3
 //'
 //'@param ATCtree : ATC tree with upper bound of the DFS (without the root),
 //'@param observations : the set of patient we calculate the distribution from
 //'@return the RR distribution among size 3 cocktail
 //'@export
 //[[Rcpp::export]]
 Rcpp::List trueDistributionSizeThreeCocktail(const DataFrame& ATCtree, const DataFrame& observations,
                                            int beta, int num_thread = 1){
   
#ifdef _OPENMP
   if(num_thread < 1 || num_thread > omp_get_num_procs()){
     std::cerr << "Wrong thread number, it should be between 1 and " 
               << omp_get_num_procs() << " \n";
     return Rcpp::List();
   }
#endif
   
   Rcpp::List observationsMedicationTmp = observations["patientATC"];
   std::vector<std::vector<int>> observationsMedication;
   observationsMedication.reserve(observationsMedicationTmp.size());
   Rcpp::LogicalVector observationsADR = observations["patientADR"];
   std::vector<int> upperBounds = ATCtree["upperBound"];
   int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
   
   for(int i =0; i < observationsMedicationTmp.size(); ++i){
     observationsMedication.push_back(observationsMedicationTmp[i]);
   }
   
   int notADRCount = observationsMedication.size() - ADRCount;
   
   //for the moment the distribution is bounded by 0 and 30
   const int distribSize = 300;
   std::vector<unsigned int> RRDistribution(distribSize);
   std::vector<unsigned int> RRDistributionGreaterBeta(distribSize);
   
   unsigned int nbCocktailNotInPopulation = 0;
   std::vector<double> outstandingRR{};
   outstandingRR.reserve(10);
   std::vector<double> outstandingRRBeta{};
   outstandingRRBeta.reserve(10);
   
   double currentRR = 0;
   std::pair<double, std::pair<int,int>> computeRROutput;
   std::pair<double, std::pair<int,int>> computePhyperOutput;
   
   int compteur = 0;
   
   int nbResults = 100; // we keep the 100 best individuals
   double minRR = 0, minRRbeta = 0;
   double minPhyper = 0, minPhyperBeta = 0;
   std::vector<std::pair<Individual,double>> bestResults;
   bestResults.reserve(nbResults);
   std::vector<std::pair<Individual,double>> bestResultsBeta;
   bestResultsBeta.reserve(nbResults);
   std::pair<Individual, double> currentResult;
   
   Individual indiv{};
   indiv.setTemperature(1);
   std::vector<int> med;
   med.resize(3);
   for(int i = 0 ; i < ATCtree.nrow() - 2 ; ++i){
     med[0] = i;
     for(int j = i+1 ; j < ATCtree.nrow() - 1; ++j){
       med[1] = j;
       for(int k = j+1 ; k < ATCtree.nrow(); ++k){
         med[2] = k;
         indiv.setMedications(med);
         computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                       upperBounds, ADRCount, notADRCount,
                                                       8000, num_thread);
         if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
           if(computePhyperOutput.first < 30){
             int index = 10 * computePhyperOutput.first;
             ++RRDistribution[index];
             if(computePhyperOutput.second.second > beta){
               ++RRDistributionGreaterBeta[index];
             }
           }
           else{
             if(computePhyperOutput.second.second > beta){
               outstandingRRBeta.push_back(computePhyperOutput.first);
             }
             // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
             outstandingRR.push_back(computePhyperOutput.first);
           }
           
           //keep the best results
           currentResult = std::make_pair(indiv, computePhyperOutput.first);
           minPhyper = addToBestCocktails(bestResults, currentResult, nbResults, minPhyper,
                                          upperBounds);
           
           if(computePhyperOutput.second.second > beta){
             minPhyperBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                                minPhyperBeta, upperBounds);
           }
         }
       }
     }
   }
   
   //create the vector of the best cocktail
   std::vector<std::vector<int>> returnedMed{};
   returnedMed.reserve(bestResults.size());
   std::vector<double>returnedRR{};
   returnedRR.reserve(bestResults.size());
   
   for(const auto &pair : bestResults){
     returnedMed.push_back(pair.first.getMedications());
     returnedRR.push_back(pair.second);
   }
   
   //create the returned vector of the best cocktail taken by more than beta person
   std::vector<std::vector<int>> returnedMedBeta{};
   returnedMedBeta.reserve(bestResultsBeta.size());
   std::vector<double>returnedRRBeta{};
   returnedRRBeta.reserve(bestResultsBeta.size());
   
   for(const auto &pair : bestResultsBeta){
     returnedMedBeta.push_back(pair.first.getMedications());
     returnedRRBeta.push_back(pair.second);
   }
   
   return Rcpp::List::create(Rcpp::Named("Distribution") = RRDistribution,
                             Rcpp::Named("FilteredDistribution") = RRDistributionGreaterBeta,
                             Rcpp::Named("OutstandingPHyper") = outstandingRR,
                             Rcpp::Named("BestCocktails") = returnedMed,
                             Rcpp::Named("BestCocktailsBeta") = returnedMedBeta,
                             Rcpp::Named("BestPHyper") = returnedRR,
                             Rcpp::Named("BestPHyperBeta") = returnedRRBeta);
}


double computeCSS(const std::pair<double, std::pair<int,int>>& RR1_and_2,
                  const std::pair<double, std::pair<int,int>>& RR1,
                  const std::pair<double, std::pair<int,int>>& RR2,
                  int ADRCount, int popSize){
  // lower .025 interval for RR of concomittant use of D1 and D2
  double RR1_2_N11 = RR1_and_2.second.first;
  double RR1_2_N1plus = RR1_and_2.second.second;
  double RR1_2_N0plus = popSize - RR1_2_N1plus;
  double RR1_2_N01 = ADRCount - RR1_2_N11;
  double tmp;
  if(RR1_2_N11 != 0 && RR1_2_N1plus != 0 && RR1_2_N0plus != 0 && RR1_2_N01 != 0)
  {
     tmp = log(RR1_and_2.first) - (1.96 * sqrt(((1.0/RR1_2_N11) 
                                            - (1.0/RR1_2_N1plus)
                                            + (1.0/RR1_2_N01)
                                            - (1.0/RR1_2_N0plus))));
  }
  else{
    tmp = 0;
  }
  
  double RR_1_et_2_interval025 = exp(tmp);
  
  // upper .975 interval for RR of Drug 1 use 
  double RR1_N11 = RR1.second.first;
  double RR1_N1plus = RR1.second.second;
  double RR1_N0plus = popSize - RR1_N1plus;
  double RR1_N01 = ADRCount - RR1_N11;
  if(RR1_N11 != 0 && RR1_N1plus != 0 && RR1_N0plus != 0 && RR1_N01 != 0)
  {
    tmp = log(RR1.first) + (1.96 * sqrt(((1.0/RR1_N11) 
                                                - (1.0/RR1_N1plus)
                                                + (1.0/RR1_N01)
                                                - (1.0/RR1_N0plus))));
  }
  else{
    tmp = 0;
  }
                                                      
  double RR1_interval975 = exp(tmp);
                                                      
    // upper .975 interval for RR of Drug 1 use 
  double RR2_N11 = RR2.second.first;
  double RR2_N1plus = RR2.second.second;
  double RR2_N0plus = popSize - RR2_N1plus;
  double RR2_N01 = ADRCount - RR2_N11;
  if(RR2_N11 != 0 && RR2_N1plus != 0 && RR2_N0plus != 0 && RR2_N01 != 0)
  {
    tmp = log(RR2.first) + (1.96 * sqrt(((1.0/RR2_N11) 
                                        - (1.0/RR2_N1plus)
                                        + (1.0/RR2_N01)
                                        - (1.0/RR2_N0plus))));
  }else{
    tmp = 0;
  }    
  
  double RR2_interval975 = exp(tmp);

  return std::max(RR1_interval975,RR2_interval975) <= 1.0 || RR_1_et_2_interval025 <= 1.0? // meaning that tmp == 0
         0 : RR_1_et_2_interval025 / std::max(RR1_interval975,RR2_interval975);
}

//'Function used to compare diverse metrics used in Disproportionality analysis
 //'
 //'@return the RR distribution among size 2 cocktail
 //'@export
 //[[Rcpp::export]]
std::vector<double> MetricCalc(const std::vector<int> &cocktail, 
                               const std::vector<int>& ATClength,
                               const std::vector<int>& upperBounds,
                               const std::vector<std::vector<int>>& observationsMedication,
                               const Rcpp::LogicalVector& observationsADR,
                               int ADRCount, int num_thread = 1){
  
  std::vector<double> solution;
  solution.reserve(5);
  
  double propADR = static_cast<double>(ADRCount) /
    static_cast<double>(observationsMedication.size());
  
  Individual ind{cocktail};
  Individual ind_D1{{cocktail[0]}};
  Individual ind_D2{{cocktail[1]}};
  //10 000 stads for the upper limit of the metrics we don't want any upper limit
  // 10 000 should never be encoutered
  auto RR_1_et_2 = ind.computeRR(observationsMedication, observationsADR, upperBounds,
                10000, num_thread);
  solution.push_back(RR_1_et_2.first);

  double phyper = ind.computePHypergeom(observationsMedication, observationsADR,
                                        upperBounds, ADRCount, 
                                        observationsMedication.size() - ADRCount,
                                        10000, num_thread).first;
  solution.push_back(phyper);
  
  // -1 because we want P(X >= D1_D2_AE_count)
  Rcpp::IntegerVector D1_D2_AE_count{RR_1_et_2.second.first - 1};
  double D1_D2_count = RR_1_et_2.second.second;
  
  double pbinom = -(Rcpp::pbinom(D1_D2_AE_count, D1_D2_count, propADR, false, true)[0]);
  solution.push_back(pbinom);
  
  auto RR_1 = ind_D1.computeRR(observationsMedication, observationsADR, upperBounds,
                               10000, num_thread);
  
  auto RR_2 = ind_D2.computeRR(observationsMedication, observationsADR, upperBounds,
                               10000, num_thread);
  
  double CRR = std::max(RR_1.first, RR_2.first) == 0 ?
                0 : RR_1_et_2.first / std::max(RR_1.first, RR_2.first);
  solution.push_back(CRR);
  
  double CSS = computeCSS(RR_1_et_2, RR_1, RR_2, ADRCount,
                          observationsMedication.size());
  solution.push_back(CSS);
  
  return solution;
}

//[[Rcpp::export]]
Rcpp::DataFrame computeMetrics(const Rcpp::DataFrame& df, const DataFrame& ATCtree,
                    const DataFrame& observations, int num_thread = 1 ){
  
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  
  
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }

  std::unordered_map<std::string, std::vector<double>> metrics{{"CSS" , {}}, 
                                                    {"CRR" , {}},
                                                    {"pbinom" , {}},
                                                    {"phyper" , {}},
                                                    {"RR" , {}}};
  std::vector<double> metricsResults;
  metricsResults.reserve(5);
  
  
  Rcpp::List CocktailList =  df["Cocktail"];
  for(const std::vector<int>& cocktail : CocktailList){
    metricsResults = MetricCalc(cocktail, ATClength, upperBounds, 
                                observationsMedication, observationsADR, ADRCount,
                                num_thread);
    int i = 0;
    for(auto& pair : metrics){
      pair.second.push_back(metricsResults[i++]);
    }
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("Label") = df["Label"],
                                 Rcpp::Named("RR") = metrics["RR"],
                                 Rcpp::Named("phyper") = metrics["phyper"],
                                 Rcpp::Named("pbinom") = metrics["pbinom"],
                                 Rcpp::Named("CRR") = metrics["CRR"],
                                 Rcpp::Named("CSS") = metrics["CSS"]);            
}  


void print_list_in_file(const Rcpp::List& resultsGeneticAlgorithm,
                        const std::string& filename){
  std::ofstream ost(filename, std::ios::app);
  if(!ost.is_open()){
    std::cerr << "erreur ouverture fichier \n";
  }
  
  Rcpp::List final_population = resultsGeneticAlgorithm["FinalPopulation"];
  std::vector<std::vector<int>> cocktails = final_population["cocktails"];
  std::vector<double> score = final_population["score"];
  
  std::vector<int> mean_score_per_epoch = resultsGeneticAlgorithm["meanFitnesses"];
  std::vector<int> best_score_per_epoch = resultsGeneticAlgorithm["BestFitnesses"];
  
  //print cocktail + score 
  int i = 0;
  for(const auto& vec : cocktails){
    for(int med : vec){
      ost << med << ' ';
    }
    ost << score[i++] << '\n';
  }
  
  //for each epoch we print "mean_score best_score"
  for(int j = 0 ; j < mean_score_per_epoch.size() ; ++j){
    ost << mean_score_per_epoch[j] << ' ' << best_score_per_epoch[j] << '\n';
  }
  ost.close();
}
//[[Rcpp::export]]
void hyperparam_test_genetic_algorithm(int epochs, int nb_individuals, 
                                       const DataFrame& ATCtree, 
                                       const DataFrame& observations,
                                       int nb_test_desired, 
                                       const std::vector<double>& mutation_rate,
                                       const std::vector<int>& nb_elite,
                                       const std::vector<double>& alphas,
                                       const std::string& path = "./",
                                       int num_thread=1){
  for(const auto& mr : mutation_rate){
    for(const auto& ne : nb_elite){
      for(const auto& alpha : alphas){
        std::ostringstream filename;
        filename << path << epochs << "e_" << nb_individuals << "ind_" << mr << "mr_" <<
          ne << "ne_" << alpha << "alpha.txt";
        for(int i = 0; i < nb_test_desired; ++i){
          auto out = GeneticAlgorithm(epochs, nb_individuals, ATCtree, observations,
                                      num_thread, true, 0.8, mr, ne, 2, alpha,
                                      false);
          print_list_in_file(out, filename.str());
        }
      }
    }
  }
}

using answer_set = std::vector<std::vector<int>>;

std::vector<int> recup_cocktail(const std::string& line){
  std::istringstream stream(line.data());
  int medoc;
  std::vector<int> returned_vec;
  
  while(stream >> medoc){
    returned_vec.push_back(medoc);
  }
  returned_vec.pop_back();
  
  return returned_vec;
}

std::pair<double, std::vector<int>> recup_solution(const std::string& line){
  std::istringstream stream(line.data());
  
  std::vector<int> vec;
  double number;
  
  while(stream >> number){
    vec.push_back(number);
  }
  
  vec.pop_back();
  std::sort(vec.begin(),vec.end());
  return {number, vec};
}

void analyze(const std::deque<std::vector<int>>& cocktail_trouves,
             const std::vector<std::vector<int>>& vraies_reponses,
             const std::string& input_filename, const std::vector<int>& depth,
             const std::vector<int> father,
             const std::string& output_filename = "analytics.txt"){
  
  std::ofstream out(output_filename, std::ios::app);
  std::deque<std::vector<int>> population_cocktails = cocktail_trouves;
  for(const auto& solution : vraies_reponses){
    population_cocktails.push_front(solution);
    Population pop_tmp({population_cocktails.begin(),population_cocktails.end()});
    
    
    IntMatrix M;
    std::vector<int> indexM;
    
    std::tie(M, indexM) = pop_tmp.pretraitement(depth,father);
    RealMatrix S = pop_tmp.similarity(M, indexM);
    
    //now we have to find the cocktail that is the most similar to the 
    // good solution
    int i_max = 1;
    double max_sim = S[0][1]; // we start at (0,1) because (0,0) is 
    // a comparaison between the real solution and himself
    for(int i = 2 ; i < population_cocktails.size(); ++i){
      if(max_sim < S[0][i]){
        i_max = i;
        max_sim = S[0][i];
      }
    }
    out << input_filename << " | ";
    for(const auto& med : solution)
      out << med << " ";
    out << "| ";
    for(const auto& med : population_cocktails[i_max])
      out << med << " ";
    out << "| " << max_sim << '\n';
    population_cocktails.pop_front();
  }
}

//[[Rcpp::export]]
void analyse_resultats(const std::vector<std::vector<int>>& reponses,
                       const std::string& input_filename,
                       int repetition, const DataFrame& ATCtree){
  
  std::ifstream input(input_filename);
  if(!input.is_open()){
    std::cerr << "erreur ouverture du fichier " << input_filename << "\n";
    return;
  }
  int epochs = std::stoi(input_filename.substr(input_filename.find('/')+1,
                                               input_filename.find('e')).data());
  int nb_individuals = std::stoi(input_filename.substr(input_filename.find('_')+1,
                                                       input_filename.find('i')).data());
  
  std::vector<std::deque<std::vector<int>>> cocktail_par_essai;
  std::string line;
  cocktail_par_essai.reserve(repetition);
  
  for(int i = 0; i < repetition; ++i){
    std::deque<std::vector<int>> tmp;
    for(int j = 0; j < nb_individuals;++j){
      std::getline(input, line);
      tmp.push_back(recup_cocktail(line));
    }
    cocktail_par_essai.push_back(tmp);
    for(int j = 0; j < epochs; ++j){
      std::getline(input, line);
    }
  }
  input.close();
  
  std::vector<int> ATClength = ATCtree["ATC_length"];
  
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  for(int i = 0; i < repetition ; ++i){
    analyze(cocktail_par_essai[i], reponses, input_filename, depth, father);
  }
}

//[[Rcpp::export]]
Rcpp::DataFrame true_results_dissimilarity_and_class(std::deque<std::vector<int>> cocktails,
                                                     const std::deque<std::vector<int>>& solutions,
                                                     const DataFrame& ATCtree){
  std::vector<short> classes;
  std::vector<double> similarity;
  
  classes.reserve(cocktails.size());
  similarity.reserve(similarity.size());
  
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  for(const auto& sol : solutions){
      cocktails.push_front(sol);
  }
  
  Population population({cocktails.begin(),cocktails.end()});
  
  IntMatrix M;
  std::vector<int> indexM;
  
  std::tie(M, indexM) = population.pretraitement(depth,father);
  RealMatrix S = population.similarity(M, indexM);
  
  for(int i = solutions.size(); i < cocktails.size(); ++i){
    int jmax = 0;
    double max_sim = S[i][jmax];
    for(int j = 1 ; j < solutions.size(); ++j){
      if(S[i][j] > S[i][jmax]){
        max_sim = S[i][j];
        jmax = j;
      }
    }
    classes.push_back(jmax);
    similarity.push_back(max_sim);
  }
  return Rcpp::DataFrame::create(Rcpp::Named("class") = classes,
                         Rcpp::Named("similarity") = similarity);
}

//[[Rcpp::export]]
std::vector<std::vector<double>> test_func(const std::vector<int>& ATClength){
  
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  std::vector<std::vector<int>> cocktails = {{1903},{521, 904, 1696, 1903, 4731}};
  Population population({cocktails.begin(),cocktails.end()});
  
  IntMatrix M;
  std::vector<int> indexM;
  
  std::tie(M, indexM) = population.pretraitement(depth,father);
  
  std::cout << "M\n";
  for(const auto& row : M){
    for(const auto& i : row){
      std::cout<<i << " ";
    }
    std::cout<<"\n";
  }
  
  std::cout << std::count(M[0].begin(),M[0].end(), M[0][0]) << std::endl;
  std::cout << std::count(M[1].begin(),M[1].end(), M[1][0]) << std::endl;
  RealMatrix D = population.dissimilarity(M, indexM, false);
  
  return D;
}

void print_most_similar(const std::deque<std::pair<double, std::vector<int>>>& sol,
                        const std::vector<std::vector<int>>& reponses,
                        const std::vector<int>& depth, const std::vector<int> father,
                        const std::string& input_filename,
                        const std::string& output_filename = "solution_rank.txt"){
  std::ofstream out(output_filename, std::ios::app);
  out << input_filename << " : \n";
  
  auto find_sol = sol;
  
  for(int i = 0; i < reponses.size(); ++i){
    find_sol.push_front({0.0, reponses[i]});
  }
  
  std::vector<std::vector<int>> vec_sol_tmp;
  vec_sol_tmp.reserve(find_sol.size());
  for(const auto& pair : find_sol){
    vec_sol_tmp.push_back(pair.second);
  }
  
  Population population_cocktail({vec_sol_tmp.begin(),vec_sol_tmp.end()});
  
  IntMatrix M;
  std::vector<int> indexM;
  
  std::tie(M, indexM) = population_cocktail.pretraitement(depth,father);
  RealMatrix S = population_cocktail.similarity(M, indexM);
  
  //we go from reponses.size() to find_sol.size() because we added the real
  // solution in the deque
  for(int i = reponses.size(); i < find_sol.size(); ++i){
    int jmax = 0;
    double max_sim = S[i][jmax];
    
    for(int j = 1; j < reponses.size(); ++j){
      if(max_sim < S[i][j]){
        max_sim = S[i][j];
        jmax = j;
      }
    }
    
    out << "Rank " << i - reponses.size() + 1 << " | ";
    for(int a : find_sol[i].second){
      out << a << ' ';
    }
    out << "| ";
    for(int a : find_sol[jmax].second){
      out << a << ' ';
    }
    out << "| similarity : " << max_sim << " | score : " << find_sol[i].first << '\n';
  }
  
}

void print_without_solutions(const std::deque<std::pair<double, std::vector<int>>>& sol,
                             const std::vector<std::string>& actives_principles,
                             const std::string& input_filename,
                             const std::string& output_filename = "solution_rank.txt"){
  std::ofstream out(output_filename, std::ios::app);
  out << input_filename << " : \n";
  
  auto find_sol = sol;
  
  //we go from reponses.size() to find_sol.size() because we added the real
  // solution in the deque
  for(int i = 0; i < sol.size(); ++i){
    
    out << "Rank " << i+1 << " | ";
    for(int a : sol[i].second){
      out << actives_principles[a] << " ; ";
    }
    out << " | score : " << sol[i].first << '\n';
  }
}


//[[Rcpp::export]]
void analyse_resultats_2(const std::vector<std::vector<int>>& reponses,
                       const std::string& input_filename,
                       int repetition, const DataFrame& ATCtree,
                       bool have_solution){
  std::ifstream input(input_filename);
  if(!input.is_open()){
    std::cerr << "erreur ouverture du fichier " << input_filename << "\n";
    return;
  }
  int epochs = std::stoi(input_filename.substr(input_filename.find('/')+1,
                                               input_filename.find('e')).data());
  int nb_individuals = std::stoi(input_filename.substr(input_filename.find('_')+1,
                                                       input_filename.find('i')).data());
  
  std::vector<std::pair<double, std::vector<int>>> solutions;
  
  std::string line;
  solutions.reserve(repetition*nb_individuals);
  
  for(int i = 0; i < repetition; ++i){
    for(int j = 0; j < nb_individuals;++j){
      std::getline(input, line);
      solutions.push_back(recup_solution(line));
    } 
    
    for(int j = 0; j < epochs; ++j){
      std::getline(input, line);
    } 
  }
  input.close();
  
  //remove all duplicates by constructing a set
  std::set<
    std::pair<double, std::vector<int>>,
    std::greater<std::pair<double, std::vector<int>>>
  > set_sol(solutions.begin(), solutions.end());
  
  solutions.clear();
  solutions.reserve(set_sol.size());
  std::move(set_sol.begin(), set_sol.end(), std::back_inserter(solutions));
  
  
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);

  if(have_solution){
    print_most_similar({solutions.begin(),solutions.end()}, reponses, depth, father,
                       input_filename);
  }
  else{
    print_without_solutions({solutions.begin(),solutions.end()}, ATCtree["Name"],
                            input_filename);
  }
}


//'Print every cocktail found during the genetic algorithm 
//'
//'@param input_filenames : A List containing filename of hyperparam_test_genetic_algorithm output file
//'@param ATCtree : The ATC tree
//'@param csv_filename : Name of the output file
//'
//'@export
//[[Rcpp::export]]
void print_csv(const std::vector<std::string>& input_filenames,
               const DataFrame& observations,
               int repetition, const DataFrame& ATCtree,
               const std::string& csv_filename = "solutions.csv" ){
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];

  
  std::vector<std::pair<double, std::vector<int>>> solutions;
  
  for(const auto& filename : input_filenames){
    std::ifstream input(filename);
    if(!input.is_open()){
      std::cerr << "erreur ouverture du fichier " << filename << "\n";
      return;
    }
    
    int epochs = std::stoi(filename.substr(filename.find('/')+1,
                                           filename.find('e')).data());
    int nb_individuals = std::stoi(filename.substr(filename.find('_')+1,
                                                   filename.find('i')).data());
    
    std::string line;
    solutions.reserve(solutions.capacity() + repetition*nb_individuals);
    
    for(int i = 0; i < repetition; ++i){
      for(int j = 0; j < nb_individuals;++j){
        std::getline(input, line);
        auto tmp = recup_solution(line);
        if(tmp.first > 0)
          solutions.push_back(tmp);
      } 
      
      for(int j = 0; j < epochs; ++j){
        std::getline(input, line);
      } 
    }
    
    input.close();
  }
  
  std::set<
    std::pair<double, std::vector<int>>,
    std::greater<std::pair<double, std::vector<int>>>
  > set_sol(solutions.begin(), solutions.end());
  
  std::ofstream output(csv_filename);
  if(!output.is_open()){
    std::cerr << "erreur ouverture du fichier " << csv_filename << "\n";
    return;
  }
  std::vector<std::string> ATCName = ATCtree["Name"];
  
  output << "score ; Cocktail ; n patient taking C ; n patient taking C and having AE ; RR \n";
  
  for(const auto& sol : set_sol){
    Individual c{sol.second};
    auto pair = c.computePHypergeom(observationsMedication, observationsADR,
                                    ATCtree["upperBound"], 1,1,1,1).second;
    
    double RR = c.computeRR(observationsMedicationTmp, observationsADR, ATCtree);
    output << sol.first << ";";
    for(auto ite = sol.second.begin(); ite != sol.second.end()-1; ++ite){
      output << ATCName[*ite] << ":"; 
    }
    output << ATCName[*(sol.second.end()-1)] << ";";
    output << pair.second << ";" << pair.first << ";" << RR << "\n";
  }
  
  output.close();
}

std::vector<std::vector<double>> dissim(const Population& pop,
                                        const std::vector<int>& depth,
                                        const std::vector<int>& father,
                                        bool normalization){
  IntMatrix M;
  std::vector<int> indexM;
  
  std::tie(M, indexM) = pop.pretraitement(depth,father);
  RealMatrix D = pop.dissimilarity(M, indexM, normalization);
  
  return D;
}

//[[Rcpp::export]]
std::vector<std::vector<double>> get_dissimilarity_from_list(const Rcpp::List& genetic_results,
                                                   const DataFrame& ATCtree){
  Rcpp::List population_list = genetic_results["FinalPopulation"];
  std::vector<std::vector<int>> cocktails = population_list["cocktails"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  Population population(cocktails);

  return dissim(population, depth, father, true);
}

//[[Rcpp::export]]
std::vector<std::vector<double>> get_dissimilarity(
                                                   const std::string& filename,
                                                   const DataFrame& ATCtree,
                                                   bool normalization = true){
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<std::vector<int>> cocktails;
  std::vector<int> current_cocktail;
  
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  std::string line;
  std::ifstream input(filename);
  if(!input.is_open()){
    std::cerr << "erreur ouverture du fichier " << filename << "\n";
    return {};
  }
  
  std::getline(input,line); // we don't want the first line as it is a title
  
  while(std::getline(input, line)){
    int beg_cocktail = line.find('|');
    int end_cocktail = line.substr(beg_cocktail+1).find('|');
    std::string cocktail = line.substr(beg_cocktail + 1, end_cocktail);
    // 3.8 is the expected number of char in the id of a randomly picked 
    // drug in the ATC tree (2014tree)
    int approximated_cocktail_size = cocktail.size() / 3.8;
    current_cocktail.reserve(approximated_cocktail_size);
    
    int med;
    std::istringstream iss(cocktail);
    while(iss >> med){
      current_cocktail.push_back(med);
    }
    cocktails.push_back(current_cocktail);
    current_cocktail.clear();
    
  }
  
  Population population(cocktails);
  
  return dissim(population, depth, father, normalization);
}


//[[Rcpp::export]]
std::vector<std::vector<double>> get_dissimilarity_from_cocktail(const std::vector<std::vector<int>>& cocktails,
                                                                const Rcpp::DataFrame& ATCtree,
                                                                bool normalization = true){
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  Population population(cocktails);
  
  return dissim(population, depth, father, normalization);
}

//[[Rcpp::export]]
Rcpp::DataFrame get_answer_class(const std::string& filename,
                                  const std::vector<std::string>& answer){
  std::vector<std::string> cocktails;
  std::string current_cocktail;
  std::vector<int> classe;
  classe.reserve(100);
  std::vector<double> similarity;
  similarity.reserve(100);
  
  std::vector<double> scores;
  scores.reserve(100);
  
  std::string line;
  std::ifstream input(filename);
  if(!input.is_open()){
    std::cerr << "erreur ouverture du fichier " << filename << "\n";
    return {};
  }
  
  std::getline(input,line); // we don't want the first line as it is a title
  
  while(std::getline(input, line)){
    //find cocktail
    int beg_cocktail = line.find('|')+1;
    int end_cocktail = line.find('|', beg_cocktail);
    cocktails.push_back(line.substr(beg_cocktail , end_cocktail-beg_cocktail));
    
    //find solution
    int beg_solution = end_cocktail+2;
    int end_solution = line.find('|', beg_solution) -1;
    std::string current_solution = line.substr(beg_solution, end_solution-beg_solution);
    
    int i =0;
    while(i < answer.size() && current_solution.compare(answer[i]) != 0){
      ++i;
    }
    classe.push_back(i);
    
    //find similarity
    int beg_similarity = line.find(':') +1;
    int end_similarity = line.find('|', beg_similarity) -1;
    similarity.push_back(std::stod(line.substr(beg_similarity, end_similarity-beg_similarity)));
    
    int beg_score = line.find(':',end_similarity +2)+1;
    scores.push_back(std::stod(line.substr(beg_score, line.size() - beg_score)));
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("cocktail") = cocktails,
                                 Rcpp::Named("class") = classe,
                                 Rcpp::Named("similarity") = similarity,
                                 Rcpp::Named("score") = scores);
}

///////////////////////// FIXED WEIGHT EM WITH COMPONENT SELECTION -> REFACTOR THE CODE IF IT WORKS WELL ///////////////////////////

arma::vec initialize_uniform_mixture_pi(int K_high){
  return arma::vec(K_high,arma::fill::ones) / K_high;
}

arma::mat initialize_mu_using_kmeans(const arma::mat& X, int K_high){
  arma::mat centers;
  bool status = arma::kmeans(centers, X.t(), K_high,arma::random_subset, 20, false);
  if(!status){
    Rcpp::stop("K-means clustering failed");
  }
  return centers.t();
}

arma::cube initialize_cov(const arma::mat& X, int K_high){
  arma::mat cov_X = arma::cov(X);
  arma::cube sigma(X.n_cols, X.n_cols, K_high);
  for(int k = 0 ; k < K_high; ++k){
    sigma.slice(k) = cov_X;
  }
  return sigma;
}

//[[Rcpp::export]]
arma::vec dmvnrm_arma(const arma::mat& X, const arma::rowvec& mean,
                      const arma::mat& sigma_k, const arma::vec& w){
  int n = X.n_rows;
  int dim = X.n_cols;
  double log2pi = std::log(2.0 * M_PI);
  
  arma::vec densities(n);
  
  
  for(int i = 0 ; i < n ; ++i){
    arma::mat rooti = arma::inv(arma::trimatu(arma::chol((sigma_k/w(i)))));
    double logdet = arma::sum(arma::log(rooti.diag()));
    double constants = (-0.5 * (dim * log2pi)) + logdet;
    
    arma::rowvec z = (X.row(i) - mean) * rooti ;
    densities(i) = constants - 0.5 * arma::dot(z,z);
  }
  
  return arma::exp(densities);
  
}

arma::rowvec update_mean(const arma::mat& X, const arma::vec& eta_k,
                         const arma::vec& w_i){
  auto schur_w_eta = w_i % eta_k;
  arma::rowvec numerator = arma::sum((schur_w_eta) % X.each_col());
  double denom = sum(schur_w_eta);
  
  return numerator / denom;
}

arma::mat compute_posteriors(const arma::mat& X, const arma::vec& pi_k,
                             const arma::mat& mu_k, const arma::cube& Sigma_k,
                             const arma::vec& w_i) {
  
  int n = X.n_rows;
  int K = pi_k.n_elem;
  arma::mat posteriors(n,K, arma::fill::zeros);
  
  for(int k = 0; k < K; ++k){
    posteriors.col(k) = pi_k(k) * dmvnrm_arma(X, mu_k.row(k), Sigma_k.slice(k),
                   w_i);
  }
  
  arma::vec row_sums = arma::sum(posteriors, 1);
  
  posteriors.each_col() /= row_sums;
  
  return posteriors;
}

arma::mat update_covariance(const arma::mat& x, const arma::vec& eta_k,
                            const arma::vec& w, const arma::rowvec& mu_k){
  arma::mat new_sigma(x.n_cols, x.n_cols, arma::fill::zeros);
  
  for (int i = 0; i < x.n_rows; ++i) {
    arma::rowvec x_i = x.row(i);
    new_sigma += (w(i) * eta_k(i)) * ((x_i - mu_k).t() * (x_i - mu_k));
  }
  double total_weight = sum(eta_k);
  return new_sigma / total_weight;
}

double compute_q_function(const arma::mat& X, const arma::vec& pi_k,
                          const arma::mat& mu_k, const arma::cube& Sigma_k,
                          const arma::vec& w_i, const arma::mat& posteriors) {
  double q_function = 0.0;
  int N = X.n_rows;
  int K = pi_k.n_elem;
  int D = X.n_cols;
  double log2pi = std::log(2*M_PI);
  for (int k = 0; k < K; ++k) {
    arma::mat Sigma_k_weighted = Sigma_k.slice(k);
    for (int i = 0; i < N; ++i) {
      arma::rowvec x_i = X.row(i);
      arma::mat adjusted_sigma = Sigma_k_weighted / w_i(i);
      arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(adjusted_sigma))));
      double logdet = 2.0 * arma::sum(arma::log(rooti.diag()));
      double constants = -0.5 * (D * log2pi + logdet);
      arma::vec z = rooti * (x_i.t() - mu_k.row(k).t());
      double log_pdf = constants - 0.5 * arma::dot(z, z);
       
      q_function += posteriors(i, k) * (log(pi_k(k)) + log_pdf);
    }
  }
   
  return q_function;
} 

//[[Rcpp::export]]
Rcpp::List FWD_EM(const arma::mat& X, int K, double eps, const arma::vec& w_i,
                  int max_steps){
  int r = 0;
  int N = X.n_rows;
  int D = X.n_cols;

  double Q_r = std::numeric_limits<double>::infinity();
  double absolute_diff_q = std::numeric_limits<double>::infinity();
  
  arma::vec pi_k = initialize_uniform_mixture_pi(K);
  arma::mat mu_k = initialize_mu_using_kmeans(X, K);
  arma::cube sigma_k = initialize_cov(X, K); // voir si on utilise la moyenne
  //des clusters en X
  
  Rcpp::List theta_min;
  arma::mat posteriors(N,K, arma::fill::zeros);
  
  do{
    //E-step
    posteriors = compute_posteriors(X, pi_k, mu_k, sigma_k, w_i);
    
    //M-step
    arma::vec N_k = arma::sum(posteriors,0).t();
    pi_k = N_k / static_cast<double>(N);
    
    for(int k = 0; k < K; ++k){
      mu_k.row(k) = update_mean(X, posteriors.col(k), w_i);
      sigma_k.slice(k) = update_covariance(X, posteriors.col(k), w_i, mu_k.row(k));
    }
    
    double Q_r_1 = compute_q_function(X, pi_k,mu_k, sigma_k, w_i, posteriors);
    
    absolute_diff_q = std::abs(Q_r_1 - Q_r);
    ++r;
    Q_r = Q_r_1;
    
  } while (r < max_steps && absolute_diff_q > eps);
  
  return Rcpp::List::create(Rcpp::Named("posteriors") = posteriors, 
                            Rcpp::Named("pi") = pi_k, 
                            Rcpp::Named("mu") = mu_k,
                            Rcpp::Named("sigma") = sigma_k,
                            Rcpp::Named("Q_r") = Q_r,
                            Rcpp::Named("iter") = r);
}