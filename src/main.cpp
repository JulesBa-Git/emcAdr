// we only include MCMC.h which pulls RcppArmadillo.h and Rcpp.h in for us
#include "MCMC.h"
#include "Population.h"
#include <iostream>
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
                            Rcpp::Named("bestRRBeta") = returnedRRBeta);
}

//'Genetic algorithm, trying to reach the best cocktail (the one which maximize
//'the fitness function, Related Risk in our case)
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

  //output the population 
  std::vector<std::vector<int>> medications;
  medications.reserve(nbIndividuals);
  std::vector<double> populationRR;
  populationRR.reserve(nbIndividuals);
  
  medications = population.getMedications();
  populationRR = population.getRR();

  
  Rcpp::List returnedPop = Rcpp::List::create(Rcpp::Named("cocktails") = medications,
                                              Rcpp::Named("RR") = populationRR);
  
  return Rcpp::List::create(Rcpp::Named("meanFitnesses") = meanRR,
                            Rcpp::Named("BestFitnesses") = bestRR,
                            Rcpp::Named("FinalPopulation") = returnedPop);
}

//'The true RR distribution of cocktail of size 2
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@return the RR distribution among size 2 cocktail
//'@export
 //[[Rcpp::export]]
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations,
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
  med.resize(2);
  for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
    med[0] = i;
    for(int j = i+1 ; j < ATCtree.nrow(); ++j){
      med[1] = j;
      indiv.setMedications(med);
      std::cout << ++compteur << '\n';
      //computeRROutput = indiv.computeRR(observationsMedication, observationsADR, upperBounds, 8000, num_thread);
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
                    const DataFrame& observations, int num_thread ){
  
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
