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


//'The MCMC method that runs the random walk on a single cocktail in order to estimate the distribution of score among cocktails of size Smax.
//'
//'@param epochs : number of steps for the MCMC algorithm
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root, also see on the github repo for an example)
//'@param observations : real observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'
//'@param temperature : starting temperature, default = 1 (denoted T in the article)
//'@param nbResults : Number of returned solution (Cocktail of size Smax with the best oberved score during the run), 5 by default
//'@param Smax : Size of the cocktail we approximate the distribution from
//'@param p_type1: probability to operate type1 mutation. Note :
//'the probability to operate the type 2 mutation is then 1 - P_type1. P_type1 must be in [0;1]. Default is .01
//'@param beta : filter the minimum number of patients that must have taken the 
//'cocktail for his risk to be taken into account in the DistributionScoreBeta default is 4
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 500
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return I no problem, return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - OutstandingScore : An array of the score greater than max_score,
//' - bestCockatils : the nbResults bests cocktails encountered during the run.
//' - bestScore : Score corresponding to the bestCocktails.
//' - FilteredDistribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - bestCocktailsBeta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - bestScoreBeta : Score corresponding to the bestCocktailsBeta.
//' - cocktailSize : Smax parameter used during the run.
//'; Otherwise the list is empty
//'@export
//[[Rcpp::export]]
Rcpp::List DistributionApproximation(int epochs, const DataFrame& ATCtree, const DataFrame& observations,
                                              int temperature = 1, int nbResults = 5, int Smax = 2,
                                              double p_type1 = .01, int beta = 4, int max_score = 500,
                                              int num_thread = 1, bool verbose = false){
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
  
  Individual cocktail, mutatedIndividual;
  
  //for the moment the distribution is bounded by 0 and RRmax
  //const int distribSize = 300;
  std::vector<unsigned int> score_distribution((max_score*10) +1); // +1 stand for every RR over the RRmax value
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstanding_score{};
  outstanding_score.reserve(10);
  
  std::vector<unsigned int> score_distribution_beta((max_score*10) +1);
  
  //used in the phyper function
  int ADRCount = 0;
  for(const auto& adr : observationsADR){
    if(adr)
      ++ADRCount;
  }
  int notADRCount = observationsMedication.size() - ADRCount;
    
  std::pair<double, std::pair<int, int>> currentGeom = std::make_pair(0.0,std::make_pair(0,0));
  std::pair<double, std::pair<int, int>> computeGeomOutput; // pair< phypergeometric, <N° of people taking the cocktail and having the ADR, N° of people taking the cocktail>>
  double minGeom = 0;
  double minGeomBeta = 0;

  //acceptance rate
  int acceptedMove = 0;
  int accepted_type1 =0;
  int accepted_type2 =0;
  int type1_move=0, type2_move=0;
  int type1_move_inF = 0, type2_move_inF = 0;
  int falseAcceptedCocktailCount = 0, falseSampledCocktailCount= 0;
  
  double RRx_k, RRy_k, pMutation, pAcceptation, pDraw;
  
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  int chosenVertexidx;
  
  std::vector<std::pair<Individual,double>> bestResults;
  bestResults.reserve(nbResults);
  std::vector<std::pair<Individual,double>> bestResultsBeta;
  bestResultsBeta.reserve(nbResults);
  
  std::pair<Individual, double> currentResult;
  
  //if p.second is greater than 0, it means that the cocktail correspond to at least one person
  auto belongToF = [](const std::pair<double,std::pair<int,int>>& p){
    return (p.second.second > 0);
  };

  //initialization (every cocktail has to contain Smax medication)
  do{
    cocktail = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
    computeGeomOutput = cocktail.computePHypergeom(observationsMedication, observationsADR,
                                                   upperBounds, ADRCount, notADRCount,
                                                   max_score, num_thread);
  } while (!belongToF(computeGeomOutput));
  currentGeom = computeGeomOutput;
  
#ifdef _OPENMP
  std::cout<< "openMP available \n";
#endif
  
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
        
        //here the current type 1 mutation consist in drawing a new cocktail of the same size
        mutatedIndividual = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
        
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_score,
                                                                num_thread);
        
        RRy_k = computeGeomOutput.first;

        //to have an overview of the explored space (not in this method for the moment)
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        
        if(belongToF(computeGeomOutput)){
          // with this mutation, our ration q(X|Y) / q(Y|X) = 1
          pAcceptation = exp(((RRy_k - RRx_k)/static_cast<double>(cocktail.getTemperature()))); 
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
        
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_score, 
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
    if(currentGeom.first < max_score){
      int index = 10 * currentGeom.first;
      ++score_distribution[index];
      // in every cases we add the RR to the "normal returned" distribution, and if
      // more than beta persons take it, we add the RR to the other ditribution named
      // score_distribution_beta
      if(currentGeom.second.second > beta){ // second.first = N° of people taking the cocktail and having the ADR
        ++score_distribution_beta[index];
      }
    }
    else{ // since we have an RR max, we just increment the last elements of the distribution
      // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
      outstanding_score.push_back(currentGeom.first); // could remove this line ?
      ++score_distribution[score_distribution.size()-1];
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
  std::vector<double>returned_score{};
  returned_score.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returned_score.push_back(pair.second);
  }
  
  //create the returned vector with the cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returned_scoreBeta{};
  returned_scoreBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returned_scoreBeta.push_back(pair.second);
  }
  
  if(verbose){
    std::cout << "acceptance rate : " << static_cast<double>(acceptedMove) / static_cast<double>(epochs)<< "\n";
    std::cout << "acceptance rate type1 mutation : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move)<< "\n";
    std::cout << "acceptance rate type2 mutation : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move)<< "\n";
    std::cout << "acceptance rate type1 mutation when the proposal is in F : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move_inF)<< "\n";
    std::cout << "acceptance rate type2 mutation when the proposal is in F : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move_inF)<< "\n";
    std::cout << "number of proposed cocktail that was taken by nobody in the population : " << nbCocktailNotInPopulation << '\n';
    std::cout << "number of false cocktail sampled : " << falseSampledCocktailCount << '\n';
    std::cout << "number of false cocktail concidered in the distribution during the run : " << falseAcceptedCocktailCount << '\n';
  }
  
  return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution, Rcpp::Named("OutstandingScore") = outstanding_score, 
                            Rcpp::Named("bestCockatils") = returnedMed, Rcpp::Named("bestScore") = returned_score,
                            Rcpp::Named("FilteredDistribution") = score_distribution_beta,
                            Rcpp::Named("bestCocktailsBeta") = returnedMedBeta,
                            Rcpp::Named("bestScoreBeta") = returned_scoreBeta,
                            Rcpp::Named("cocktailSize") = Smax);
}

//'Genetic algorithm, trying to reach riskiest cocktails (the ones which maximize
//'the fitness function, Hypergeometric score in our case)
//'
//'@param epochs : number of step or the algorithm 
//'@param nbIndividuals : size of the population
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : real observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'@param diversity : enable the diversity mechanism of the algorithm
//' (favor the diversity of cocktail in the population),  default is false
//'@param p_crossover: probability to operate a crossover on the crossover phase. Default is 80%
//'@param p_mutation: probability to operate a mutation after the crossover phase. Default is 1%
//'@param nbElite : number of best individual we keep from generation to generation. Default is 0
//'@param tournamentSize : size of the tournament (select the best individual 
//'between tournamentSize sampled individuals) 
//'
//'@return If no problem, return a List :
//' - meanFitnesses : The mean score of the population at each epochs of the algorithm.
//' - BestFitnesses : The best score of the population at each epochs of the algorithm.
//' - FinalPopulation : The final population of the algorithm when finished (medications
//' and corresponding scores)
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
    
    //2nd : make a selection and apply the crossover on the selected individual
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


//'The true distribution of the score among every single nodes of the ATC
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
 //'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
 //' on which we want to compute the risk distribution
//'@param beta : minimum number of person taking the cocktails in order to consider it
//'in the beta score distribution 
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 1000
//'@param nbResults : Number of returned solution (Cocktail with the
//' best oberved score during the run), 100 by default
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return Return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - Filtered_score_distribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - Outstanding_score : An array of the score greater than max_score,
//' - Best_cocktails : the nbResults bests cocktails encountered during the run.
//' - Best_cocktails_beta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - Best_scores : Score corresponding to the Best_cocktails.
//' - Best_scores_beta : Score corresponding to the Best_cocktails_beta.
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionDrugs(const DataFrame& ATCtree, const DataFrame& observations,
                                 int beta, int max_score = 1000, int nbResults = 100,
                                 int num_thread = 1){
 
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
 
 
 //for the moment the distribution is bounded by 0 and 1000
 const int distribSize = max_score *10;
 std::vector<unsigned int> score_distribution(distribSize);
 std::vector<unsigned int> score_distributionGreaterBeta(distribSize);
 
 std::vector<double> outstanding_score{};
 outstanding_score.reserve(10);
 std::vector<double> outstanding_scoreBeta{};
 outstanding_scoreBeta.reserve(10);
 
 std::pair<double, std::pair<int,int>> computePhyperOutput;
 
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
   computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                 upperBounds, ADRCount, notADRCount,
                                                 8000, num_thread);
   if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
     if(computePhyperOutput.first < max_score){
       int index = 10 * computePhyperOutput.first;
       ++score_distribution[index];
       if(computePhyperOutput.second.second > beta){
         ++score_distributionGreaterBeta[index];
       }
     }
     else{
       if(computePhyperOutput.second.second > beta){
         outstanding_scoreBeta.push_back(computePhyperOutput.first);
       }
       outstanding_score.push_back(computePhyperOutput.first);
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
 std::vector<double>returned_score{};
 returned_score.reserve(bestResults.size());
 
 for(const auto &pair : bestResults){
   returnedMed.push_back(pair.first.getMedications());
   returned_score.push_back(pair.second);
 }
 
 //create the returned vector of the best cocktail taken by more than beta person
 std::vector<std::vector<int>> returnedMedBeta{};
 returnedMedBeta.reserve(bestResultsBeta.size());
 std::vector<double>returned_scoreBeta{};
 returned_scoreBeta.reserve(bestResultsBeta.size());
 
 for(const auto &pair : bestResultsBeta){
   returnedMedBeta.push_back(pair.first.getMedications());
   returned_scoreBeta.push_back(pair.second);
 }
 
 return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution,
                           Rcpp::Named("Filtered_score_distribution") = score_distributionGreaterBeta,
                           Rcpp::Named("Outstanding_score") = outstanding_score,
                           Rcpp::Named("Best_cocktails") = returnedMed,
                           Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                           Rcpp::Named("Best_scores") = returned_score,
                           Rcpp::Named("Best_scores_beta") = returned_scoreBeta);
}

//'The true distribution of the score among every size-two cocktails
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param beta : minimum number of person taking the cocktails in order to consider it
//'in the beta score distribution 
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 1000
//'@param nbResults : Number of returned solution (Cocktail with the
//' best oberved score during the run), 100 by default
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return Return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - Filtered_score_distribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - Outstanding_score : An array of the score greater than max_score,
//' - Best_cocktails : the nbResults bests cocktails encountered during the run.
//' - Best_cocktails_beta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - Best_scores : Score corresponding to the Best_cocktails.
//' - Best_scores_beta : Score corresponding to the Best_cocktails_beta.
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations,
                                           int beta, int max_score = 100, int nbResults = 100,
                                           int num_thread = 1){
  
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
  const int distribSize = max_score *10;
  std::vector<unsigned int> score_distribution(distribSize);
  std::vector<unsigned int> score_distributionGreaterBeta(distribSize);
  
  std::vector<double> outstanding_score{};
  outstanding_score.reserve(10);
  std::vector<double> outstanding_scoreBeta{};
  outstanding_scoreBeta.reserve(10);
  
  std::pair<double, std::pair<int,int>> computePhyperOutput;
  
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
        if(computePhyperOutput.first < max_score){
          int index = 10 * computePhyperOutput.first;
          ++score_distribution[index];
          if(computePhyperOutput.second.second > beta){
            ++score_distributionGreaterBeta[index];
          }
        }
        else{
          if(computePhyperOutput.second.second > beta){
            outstanding_scoreBeta.push_back(computePhyperOutput.first);
          }
          // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
          outstanding_score.push_back(computePhyperOutput.first);
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
  std::vector<double>returned_score{};
  returned_score.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returned_score.push_back(pair.second);
  }
  
  //create the returned vector of the best cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returned_scoreBeta{};
  returned_scoreBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returned_scoreBeta.push_back(pair.second);
  }
  
  return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution,
                            Rcpp::Named("Filtered_score_distribution") = score_distributionGreaterBeta,
                            Rcpp::Named("outstanding_score") = outstanding_score,
                            Rcpp::Named("Best_cocktails") = returnedMed,
                            Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                            Rcpp::Named("Best_scores") = returned_score,
                            Rcpp::Named("Best_scores_beta") = returned_scoreBeta);
}


//'Function used to compare diverse metrics used in Disproportionality analysis
//'
//'@return RR and hypergeometric score among size 3 cocktail in "cocktail"
//'@export
//[[Rcpp::export]]
std::vector<double> MetricCalc_size3(const std::vector<int> &cocktail, 
                                const std::vector<int>& ATClength,
                                const std::vector<int>& upperBounds,
                                const std::vector<std::vector<int>>& observationsMedication,
                                const Rcpp::LogicalVector& observationsADR,
                                int ADRCount, int num_thread = 1){
 
 std::vector<double> solution;
 solution.reserve(2);
 
 Individual ind{cocktail};
 

 double RR_cocktail = 0.0;

 RR_cocktail = ind.computeRR(observationsMedication, observationsADR, upperBounds,
                              1000000, num_thread).first;
 
 
 //phyper
 double phyper = ind.computePHypergeom(observationsMedication, observationsADR,
                                       upperBounds, ADRCount, 
                                       observationsMedication.size() - ADRCount,
                                       10000, num_thread).first;

 solution.push_back(phyper);
 solution.push_back(RR_cocktail);
 
 return solution;
}

//'Function used to compare diverse metrics used in Disproportionality analysis
 //'
 //'@return multiple risk among size 2 cocktail
 //'@export
 //[[Rcpp::export]]
 std::vector<double> MetricCalc_2(const std::vector<int> &cocktail, 
                                  const std::vector<int>& ATClength,
                                  const std::vector<int>& upperBounds,
                                  const std::vector<std::vector<int>>& observationsMedication,
                                  const Rcpp::LogicalVector& observationsADR,
                                  int ADRCount, int num_thread = 1){
   
   std::vector<double> solution;
   solution.reserve(5);
   
   Individual ind{cocktail};
   Individual ind_D1{{cocktail[0]}};
   Individual ind_D2{{cocktail[1]}};
   
   
   int n000 = 0,n001 = 0;
   int n100 = 0, n101 = 0;
   int n011 = 0, n010 = 0;
   int n111 = 0, n110 = 0;
   bool have_D1 = false, have_D2 = false;
   for(int i = 0 ; i < observationsMedication.size(); ++i){
     have_D1 = ind_D1.matches(observationsMedication[i], upperBounds);
     have_D2 = ind_D2.matches(observationsMedication[i], upperBounds);
     if(!have_D1 && !have_D2){
       if(observationsADR[i]){
         ++n001;
       }
       else{
         ++n000;
       }
     } 
     else if(have_D1){
       if(have_D2){
         if(observationsADR[i])
           ++n111;
         else
           ++n110;
       }else{
         if(observationsADR[i])
           ++n101;
         else
           ++n100;
       }
     }
     else{
       if(observationsADR[i]){
         ++n011;
       } 
       else{
         ++n010;
       } 
     }
   } 
   
   
   double D1_rate = static_cast<double>(n101) / (static_cast<double>(n100 + n101));
   double D2_rate = static_cast<double>(n011) / (static_cast<double>(n010 + n011));
   double background_rate = static_cast<double>(n001) / (static_cast<double>(n000 + n001));
   double D1_D2_rate = static_cast<double>(n111) / (static_cast<double>(n110 + n111));
   int n00 = n001 + n000;
   int n10 = n101 + n100;
   int n01 = n011 + n010;
   int n11 = n111 + n110;
   
   // RR
   double RR_cocktail = 0.0;
   if(n111 > 0){
     RR_cocktail = (D1_D2_rate) / 
       ((n001 + n011 + n101) / (static_cast<double>(n00 + n10 + n01)));
   }
   
   
   // PRR
   double signal_PRR = 0;
   
   double RR1 = ((n101 + n111) / (static_cast<double>(n11 + n10))) /
     ((n011 + n001) / (static_cast<double>(n01 + n00)));
   
   double RR2 = ((n011 + n111) / (static_cast<double>(n11 + n01))) /
     ((n101 + n001) / (static_cast<double>(n10 + n00)));
   
   double SD_D1;
   if((n101 != 0 || n111 !=0) && (n001 !=0 || n011 !=0)){
     SD_D1 = std::sqrt((1.0/(static_cast<double>(n101+n111))) -
       (1.0 / (static_cast<double>(n10+n11))) + (1.0/(static_cast<double>(n001+n011))) -
       (1.0/(static_cast<double>(n01+n00))));
   }else{
     SD_D1 = 0;
   }
   
   double SD_D2;
   if((n101 != 0 || n001 !=0) && (n111 !=0 || n011 !=0)){
     SD_D2 = std::sqrt((1.0/(static_cast<double>(n011+n111))) -
       (1.0 / (static_cast<double>(n01+n11))) + (1.0/(static_cast<double>(n001+n101))) -
       (1.0/(static_cast<double>(n10+n00))));
   }else{
     SD_D2 = 0;
   }
   
   double SD_D1D2;
   if((n101 + n001 + n011 !=0) && (n111 !=0)){
     SD_D1D2 = std::sqrt((1.0/(static_cast<double>(n111))) -
       (1.0 / (static_cast<double>(n11))) + (1.0/(static_cast<double>(n001+n101+n011))) -
       (1.0/(static_cast<double>(n10+n00+n01))));
   }else{
     SD_D1D2 = 0;
   }
   
   double PRR_D1_025 = std::exp((std::log(RR1) - 1.96*SD_D1));
   double PRR_D2_025 = std::exp((std::log(RR2) - 1.96*SD_D2));
   double PRR_D1D2_025 = std::exp((std::log(RR_cocktail) - 1.96*SD_D1D2));
   if(n111 > 0)
     signal_PRR = PRR_D1D2_025 > std::max(PRR_D1_025, PRR_D2_025) ? 1 : 0;
   else
     signal_PRR = 0;
   
   
   //CSS
   double PRR_D1_975 = std::exp((std::log(RR1) + 1.96*SD_D1));
   double PRR_D2_975 = std::exp((std::log(RR2) + 1.96*SD_D2));
   double CSS = 0.0;
   if(n111 > 0)
     CSS = PRR_D1D2_025 / std::max(PRR_D1_975,PRR_D2_975);
   
   
   
   //omega
   double odds_ratio_r00 = background_rate / (1 - background_rate);
   double odds_ratio_r10 = D1_rate / (1 - D1_rate);
   double odds_ratio_r01 = D2_rate / (1 - D2_rate);
   
   double s11 = 1 - (1.0 / (std::max(odds_ratio_r00, odds_ratio_r10) + 
                     std::max(odds_ratio_r00, odds_ratio_r01) - 
                     odds_ratio_r00 + 1));

   double omega = std::log2((n111 + 0.5) / (s11*n11 + 0.5));
   Rcpp::NumericVector tmp{0.975};
   double lower_bound_omegaIC = INT_MIN;
   if(n111 > 0)
     lower_bound_omegaIC = omega - (Rcpp::qnorm(tmp)[0] / (std::log(2) * std::sqrt(n111)));
   
   
   //phyper
   double phyper = ind.computePHypergeom(observationsMedication, observationsADR,
                                         upperBounds, ADRCount, 
                                         observationsMedication.size() - ADRCount,
                                         10000, num_thread).first;
   
   solution.push_back(lower_bound_omegaIC);
   solution.push_back(n110);
   solution.push_back(signal_PRR);
   solution.push_back(RR_cocktail);
   solution.push_back(n111);
   solution.push_back(phyper);
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

  std::unordered_map<std::string, std::vector<double>> metrics{
    {"n110", {}},
    {"n111", {}},
    {"RR" , {}}, 
    {"PRR" , {}},
    {"CSS" , {}},
    {"omega_025" , {}},
    {"phyper" , {}}};
  
  std::vector<double> metricsResults;
  metricsResults.reserve(metrics.size());
  
  Rcpp::List CocktailList =  df["Cocktail"];
  for(auto& cocktail : CocktailList){
    metricsResults = MetricCalc_2(cocktail, ATClength, upperBounds, 
                                observationsMedication, observationsADR, ADRCount,
                                num_thread);
   
    int i = 0;
    for(auto& pair : metrics){
      pair.second.push_back(metricsResults[i++]);
    }
  }
  auto label = df["Label"];
  return Rcpp::DataFrame::create(Rcpp::Named("Label") = label,
                                 Rcpp::Named("RR") = metrics["RR"],
                                 Rcpp::Named("phyper") = metrics["phyper"],
                                 Rcpp::Named("PRR") = metrics["PRR"],
                                 Rcpp::Named("CSS") = metrics["CSS"],
                                 Rcpp::Named("omega_025") = metrics["omega_025"],
                                 Rcpp::Named("n110") = metrics["n110"],
                                 Rcpp::Named("n111") = metrics["n111"]);            
}  


//[[Rcpp::export]]
Rcpp::DataFrame computeMetrics_size3(const Rcpp::DataFrame& df, const DataFrame& ATCtree,
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
  
  std::unordered_map<std::string, std::vector<double>> metrics{{"RR" , {}}, 
                                                               {"phyper" , {}}};
  std::vector<double> metricsResults;
  metricsResults.reserve(metrics.size());
  
  
  Rcpp::List CocktailList =  df["Cocktail"];
  for(auto& cocktail : CocktailList){
    metricsResults = MetricCalc_size3(cocktail, ATClength, upperBounds, 
                                  observationsMedication, observationsADR, ADRCount,
                                  num_thread);
    int i = 0;
    for(auto& pair : metrics){
      pair.second.push_back(metricsResults[i++]);
    }
    
  } 
  auto label = df["Label"];
  return Rcpp::DataFrame::create(Rcpp::Named("Label") = label,
                                 Rcpp::Named("RR") = metrics["RR"],
                                 Rcpp::Named("phyper") = metrics["phyper"]);            
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

///////////////////// 

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
    
    auto RR = c.computeRR(observationsMedication, observationsADR, 
                            ATCtree["upperBound"], 100000);
    output << sol.first << ";";
    for(auto ite = sol.second.begin(); ite != sol.second.end()-1; ++ite){
      output << ATCName[*ite] << ":"; 
    }
    output << ATCName[*(sol.second.end()-1)] << ";";
    output << pair.second << ";" << pair.first << ";" << RR.first << "\n";
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

///////////////////////// FIXED WEIGHT EM WITH COMPONENT SELECTION ///////////////////////////

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