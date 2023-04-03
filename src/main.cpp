// we only include MCMC.h which pulls RcppArmadillo.h and Rcpp.h in for us
#include "MCMC.h"
#include "Population.h"
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

//'The Evolutionary MCMC method that runs the random walk
 //'
 //'@param epochs : number of step 
 //'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
 //'@param observation : real observation of the ADR based on the medications of each real patients
 //'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
 //'
 //'@param temperature : starting temperature, default = randomly initialized
 //'@param P_type1: probability to operate type1 mutation. Note :
 //'the probability to operate the type 2 mutation is then 1 - P_type1. P_type1 must be in [0;1]. 
 //'@param alpha : a hyperparameter allowing us to manage to probability of adding a drug to the cocktail. The probability
 //' to add a drug to the cocktail is the following : \eqn{\frac12}{\alpha/n} Where n is the original size of the cocktail. 1 is the default value.
 //'
 //'@return if no problem return an array of the approximation of the RR distribution : the distribution of RR we've met; Otherwise the list is empty
 //'@export
 //[[Rcpp::export]]
Rcpp::IntegerVector DistributionApproximation(int epochs, const DataFrame& ATCtree, const DataFrame& observations,
                                              int temperature_M1 = 1, int temperature_M2 = 1, int Smax = 4, double p_type1 = .01){
  //arguments verification
  if(p_type1 > 1 || p_type1 < 0 || epochs < 1){
    std::cerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::IntegerVector();
  }
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
  double p_type2 = 1 - p_type1;
  Individual cocktail, mutatedIndividual;
  
  //for the moment the distribution is bounded by 0 and 30
  const int distribSize = 300;
  std::vector<unsigned int> RRDistribution(distribSize);
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstandingRR{};
  outstandingRR.reserve(10);
  double currentRR = 0;
  std::pair<double, bool> computeRROutput;
  
  //acceptance rate
  int acceptedMove = 0;
  int accepted_type1 =0;
  int accepted_type2 =0;
  int type1_move=0, type2_move=0;
  
  double RRx_k, RRy_k, q_y_given_x, q_x_given_y, pMutation, pAcceptation, pDraw;
  double q_ratio;
  int seqLen;
  
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  int chosenVertexidx;
  
  //if p.second is true, it means that the cocktail correspond to at least one person
  auto belongToF = [](const std::vector<int>& med, int acceptedSize, const std::pair<double,bool>& p){
    return ((med.size() == acceptedSize) && (p.second));
  };
  
  
  //initialization (every cocktail has to contain Smax medication)
  do{
    cocktail = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature_M1)[0];
    computeRROutput = cocktail.computeRR(observationsMedication, observationsADR, ATCtree, true);
  } while (!belongToF(cocktail.getMedications(), Smax, computeRROutput));
  currentRR = computeRROutput.first;
  
  for(int i = 0; i < epochs; ++i){
      pMutation = Rcpp::runif(1,0,1)[0];
      
      if(pMutation < p_type1){
        //type 1 mutation
        RRx_k = currentRR;
        
        //mutatedIndividual = type1Mutation(cocktail, ATCtree.nrow(), alpha, emptySeq);
        //mutatedIndividual = adjustedType1Mutation(cocktail, ATCtree.nrow(), alpha, emptySeq);
        //here the current type 1 mutation consist in drawing a new cocktail of the same size
        mutatedIndividual = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature_M1)[0];
        computeRROutput = mutatedIndividual.computeRR(observationsMedication, observationsADR, ATCtree, true);
        std::cout << "size : " << mutatedIndividual.getMedications().size() << '\n';
        RRy_k = computeRROutput.first;
        std::cout << "prev RR : " << RRx_k << " new RR : " << RRy_k << '\n';

        //to have an overview of the explored space (not in this method for the moment)
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        

        
        if(belongToF(mutatedIndividual.getMedications(), Smax, computeRROutput)){
          // with this mutation, our ration q(X|Y) / q(Y|X) = 1
          pAcceptation = exp(((RRy_k - RRx_k)/static_cast<double>(temperature_M1))); //cocktail.getTemperature()
          std::cout << "proba d'acceptation : " << pAcceptation << '\n';
          pDraw = Rcpp::runif(1,0,1)[0];
          
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentRR = RRy_k;
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

        RRx_k = currentRR;
        //get every vertex 0/1 for this patient
        vertexX = cocktail.getVertexList(ATCtree);
        
        chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
        chosenVertexidx = chosenVertexidx == vertexX.size() ? vertexX.size()-1 : chosenVertexidx;
        
        std::pair<int,int> chosenVertex = vertexX[chosenVertexidx];
        
        mutatedIndividual = type2Mutation(cocktail, ATCtree.nrow(), chosenVertex);
        
        computeRROutput = mutatedIndividual.computeRR(observationsMedication, observationsADR, ATCtree, true);
        RRy_k = computeRROutput.first;
        vertexY = mutatedIndividual.getVertexList(ATCtree);
        
        //to have an overview of the explored space
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        
        if(belongToF(mutatedIndividual.getMedications(), Smax, computeRROutput)){
          pAcceptation = (exp(((RRy_k - RRx_k) / temperature_M2))) * 
            (static_cast<double>(vertexX.size())/static_cast<double>(vertexY.size()));
          
          pDraw = Rcpp::runif(1,0,1)[0];
          
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentRR = RRy_k;
            ++acceptedMove;
            ++accepted_type2;
          }
        }else{
          ++nbCocktailNotInPopulation;
        }
        ++type2_move;
      
    }
    if(currentRR <= 30){
      int index = 10 * currentRR;
      ++RRDistribution[index];
    }
    else{
      // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
      outstandingRR.push_back(currentRR);
    }
  }
  std::cout << "acceptance rate : " << static_cast<double>(acceptedMove) / static_cast<double>(epochs)<< "\n";
  std::cout << "acceptance rate type1 mutation : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move)<< "\n";
  std::cout << "acceptance rate type2 mutation : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move)<< "\n";
  std::cout << "number of proposed cocktail that was taken by nobody in the population : " << nbCocktailNotInPopulation << '\n';
  Rcpp::IntegerVector result = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(RRDistribution));
  return result;
}

//'The Evolutionary MCMC method that runs the random walk
 //'
 //'@param epochs : number of step 
 //'@param nbIndividuals : size of the popultation
 //'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
 //'@param observation : real observation of the ADR based on the medications of each real patients
 //'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
 //'@param p_crossover: probability to operate a crossover on the crossover phase.
 //'@param p_mutation: probability to operate a mutation after the crossover phase.
 //'@param nbElite : number of best individual we keep from generation to generation
 //'@param tournamentSize : size of the tournament (select the best individual be tween tournamentSize individuals) 
 //'
 //'@return if no problem return the best cocktail we found (according to the fitness function which is the Relative Risk)
 //'@export
 //[[Rcpp::export]]
Rcpp::List GeneticAlgorithm(int epochs, int nbIndividuals, const DataFrame& ATCtree, const DataFrame& observations,
                            double p_crossover = .80, double p_mutation = .01, int nbElite = 0, int tournamentSize = 2){
  //arguments verification
  if(p_crossover > 1 || p_crossover < 0 || nbIndividuals < 1 || p_mutation > 1 || p_mutation < 0 || epochs < 1){
    std::cerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::List();
  }
  if(tournamentSize > nbIndividuals || nbElite > nbIndividuals){
    std::cerr << "the tournament size and the nbElite parameters must be less equal than the number of individuals \n";
    return Rcpp::List();
  }
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];

  //generate the initial population randomly (do we consider an Smax ?)
  double meanMedicationPerObs = meanMedications(observationsMedication) - 1;
  Population population(ATCtree.nrow(), nbIndividuals, meanMedicationPerObs);
  Population matingPool(nbIndividuals);
  
  std::vector<double> meanRR;
  meanRR.reserve(epochs);
  std::vector<double> bestRR;
  bestRR.reserve(epochs);
  int bestIndividualIndex;
  
  int remainingDraw = 0;
  
  //here we may want to have a more sophisticated stopping condition (like, if the RR is 
  //significantly high given the previous calculated ditribution)
  for(int i =0; i < epochs; ++i){
    //1st : fit every individual
    population.evaluate(observationsMedication, observationsADR, ATCtree);
 
    //2nd : make a selection given every individual fitness and apply the crossover on the selected individual
    //keep the elite first
    population.keepElite(nbElite, matingPool);

        //select the individual according to the fitness
    remainingDraw = nbIndividuals - nbElite;
    population.tournamentSelection(tournamentSize, matingPool, remainingDraw);

        //operate a crossover over the mating pool 
    matingPool.crossover(nbElite, ATClength, upperBounds, ATCtree, p_crossover);

    //operate a mutation over the mating pool
    matingPool.mutate(nbElite, p_mutation, ATCtree);
    //3rd : replace the population
    population = matingPool;
    matingPool.clear();
    
    meanRR.push_back(population.getMean());
    bestIndividualIndex = population.bestIndividual();
    bestRR.push_back(population.getIndividuals()[bestIndividualIndex].first);
    population.printSummary(i, meanRR[i], bestIndividualIndex);
    
  }
  population.evaluate(observationsMedication, observationsADR, ATCtree);
  //do we keep the results if there is only more than a single medication ?
  
  return Rcpp::List::create(Rcpp::Named("meanFitnesses") = meanRR, Rcpp::Named("BestFitnesses") = bestRR);
}

//'The true RR distribution of cocktail of size 2
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@return the RR distribution among size 2 cocktail
//'@export
 //[[Rcpp::export]]
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations){
  Rcpp::List observationsMedication = observations["patientATC"];
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  
  //for the moment the distribution is bounded by 0 and 30
  const int distribSize = 300;
  std::vector<unsigned int> RRDistribution(distribSize);
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstandingRR{};
  outstandingRR.reserve(10);
  double currentRR = 0;
  std::pair<double, bool> computeRROutput;
  int compteur = 0;
  
  Individual indiv{};
  indiv.setTemperature(1);
  std::vector<int> med;
  med.resize(2);
  for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
    med[0] = i;
    for(int j = i+1 ; j < ATCtree.nrow(); ++j){
      med[1] = j;
      indiv.setMedications(med);
      std::cout << ++compteur << '\n';
      computeRROutput = indiv.computeRR(observationsMedication, observationsADR, ATCtree, true);
      if(computeRROutput.second){
        if(computeRROutput.first <= 30){
          int index = 10 * computeRROutput.first;
          ++RRDistribution[index];
        }
        else{
          // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
          outstandingRR.push_back(computeRROutput.first);
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("Distribution") = RRDistribution, Rcpp::Named("OutstandingRR") = outstandingRR);
}
  