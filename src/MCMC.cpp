#include "MCMC.h"

// [[Rcpp::depends(RcppArmadillo)]]

double meanMedications(const Rcpp::List& observations){
  double sum = 0;
  for(const std::vector<int> ind : observations){
    sum += ind.size();
  }
  
  return (sum / observations.size());
}

std::pair<Individual,double> largerRR(const std::pair<Individual,double>& firstRR, const std::pair<Individual,double>& secRR,
                                      const std::pair<Individual,double>& thirdRR, const std::pair<Individual,double>& fourthRR){
  if(firstRR.second >= secRR.second && firstRR.second >= thirdRR.second && firstRR.second >= fourthRR.second)
    return firstRR;
  else if(secRR.second >= firstRR.second && secRR.second >= thirdRR.second && secRR.second >= fourthRR.second)
    return secRR;
  else if(thirdRR.second >= firstRR.second && thirdRR.second >= secRR.second &&thirdRR.second >= fourthRR.second)
    return thirdRR;
  else
    return fourthRR;
}

void addRRtoDistribution(const double RR,std::vector<unsigned int>& vec){
  double dIndex = (RR-1) * 10;
  int index = static_cast<int>(dIndex);
  vec[index]++;
}

void addPairToSet(const Individual& i, std::set<std::pair<int,int>>& p){
  int minPair,maxPair;
  
  std::vector<int> mutIndivTmp = i.getMedications();
  if(mutIndivTmp.size() == 2){
    if(mutIndivTmp[0] > mutIndivTmp[1]){
      maxPair = mutIndivTmp[0];
      minPair = mutIndivTmp[1];
    }else{
      maxPair = mutIndivTmp[1];
      minPair = mutIndivTmp[0];
    }
    p.insert(std::make_pair(minPair, maxPair));
  }
}

bool isNotInResultList(const std::vector<std::pair<Individual,double>>& bestResults,
                       const std::pair<Individual,double>& bestResult){
  return std::find(bestResults.begin(),bestResults.end(),bestResult) == bestResults.end();
}


std::vector<Individual> DFtoCPP_WOtemp(const Rcpp::List& startingInd){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(startingInd.size());
  Rcpp::NumericVector temp;
  int itemp = 1;
  
  for(const std::vector<int> ind : startingInd){
    temp = Rcpp::runif(1,itemp-1,itemp);
    returnedVec.push_back(Individual{ind,temp[0]});
    //because the temperatures have to be increasing
    ++itemp;
  }
  return returnedVec;
}

std::vector<Individual> DFtoCPP_Wtemp(const Rcpp::List& startingInd,const Rcpp::NumericVector& startingTemp){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(startingInd.size());

  for(int i = 0 ; i < startingInd.length(); ++i){
    returnedVec.push_back(Individual{startingInd[i],startingTemp[i]});
  }
  return returnedVec;
  
}

std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals,double meanMedic){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int nbMedic,num, itemp = 1;
  double temp;
  std::vector<int> medicVec(0);
  
  for(int i = 0; i < nbIndividuals; ++i){
    nbMedic = 1 + Rcpp::rpois(1,meanMedic)[0];
    medicVec.reserve(nbMedic);
    
    for(int j = 0; j < nbMedic ; ++j){
      //Draw every medications for a patient -> consider using sample ?
      num = trunc(Rcpp::runif(1,0,treeSize)[0]);
      num = num == treeSize ? treeSize-1 : num;
      medicVec.push_back(num);
    }
    temp = Rcpp::runif(1,itemp-1,itemp)[0];
    
    returnedVec.push_back(Individual{medicVec,temp});
    ++itemp;  
    medicVec.clear();
  }
  return returnedVec;
}

std::set<std::pair<int,int>> getADRPairs(const Rcpp::List& observationsMed, const Rcpp::LogicalVector& ADR){
  std::set<std::pair<int,int>> retSet;
  int i = 0;
  int max, min;
  for(const std::vector<int> tab : observationsMed){
    
    if(tab.size() == 2 && ADR[i]){
      if(tab[0] > tab[1]){
        max = tab[0];
        min = tab[1];
      }
      else{
        max = tab[1];
        min = tab[0];
      }
      retSet.insert(std::make_pair(min,max));
    }
    ++i;
  }

  return retSet;
}

Individual type1Mutation(const Individual& indiv, int treeSize, double alpha, bool emptyCocktail){
  //peut optimiser les appels à getMedication()
  int mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
  mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
  std::vector<int> newMed = indiv.getMedications();
  
  if(!emptyCocktail){
    
    double addAcceptation = (alpha/indiv.getMedications().size());
    double draw = Rcpp::runif(1,0,1)[0];
    
    if(addAcceptation >= draw){
      while (std::find(newMed.begin(), newMed.end(), mutateMed) != std::end(newMed)){
        // we draw a medication that is not inside the cocktail for the moment
        mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
        mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
      }
      
      newMed.push_back(mutateMed);
    }
    else{ // here we remove an element from the medications of the individuals
      int removeIndex = trunc(Rcpp::runif(1,0,indiv.getMedications().size())[0]);
      removeIndex = removeIndex == indiv.getMedications().size() ? removeIndex-1 : removeIndex;
      newMed.erase(newMed.begin() + removeIndex);
    }
    
  }
  else{ // here the cocktail is empty so we add a drug with probability 1 
    newMed.push_back(mutateMed);
  }
  
  return {newMed,indiv.getTemperature()};
}

Individual type2Mutation(const Individual& indiv, int treeSize, const std::pair<int,int>& p){
  // if a patient got every medication
  if(indiv.getMedications().size() == treeSize)
    return indiv;
  
  std::vector<int> prevMedic = indiv.getMedications();
  auto rmEnd = std::remove(prevMedic.begin(),prevMedic.end(),p.first);
  prevMedic.erase(rmEnd,prevMedic.end());
  prevMedic.push_back(p.second);
  
  return {prevMedic,indiv.getTemperature()};
}

Individual crossoverMutation(const Individual& indiv1, const Individual& indiv2,const Rcpp::DataFrame& ATCtree,
                             int selectedNode, int upperBound){
  std::vector<int> newMedi{};
  newMedi.reserve(indiv1.getMedications().size() + indiv2.getMedications().size());
  
  for(int med : indiv1.getMedications()){
    if(med < selectedNode || med >= upperBound){
      newMedi.push_back(med);
    }
  }
  for(int med : indiv2.getMedications()){
    if(med >= selectedNode && med < upperBound){
      newMedi.push_back(med);
    }
  }
  newMedi.shrink_to_fit();
  
  return {newMedi, indiv1.getTemperature()};
}

