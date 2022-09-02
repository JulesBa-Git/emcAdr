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

Individual type1Mutation(const Individual& indiv, int treeSize){
  double addAcceptation = (1/indiv.getMedications().size());
  double draw = Rcpp::runif(1,0,1)[0];
  std::vector<int> newMed = indiv.getMedications();
  
  if(addAcceptation > draw){
    int mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
    mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
    
    //note when changing the mutation 1 : do we have to keep that ? 
    std::vector<int>::iterator find = std::find(newMed.begin(), newMed.end(), mutateMed);
    //the medication is in the vector so we erase it 
    if( find != newMed.end()){
      newMed.erase(find);
    } else{ // the medication is not in the vector so we add it
      newMed.push_back(mutateMed);
    }
  }
  else{ // here we remove an element from the medications of the individuals
    int removeIndex = trunc(Rcpp::runif(1,0,indiv.getMedications().size())[0]);
    removeIndex = removeIndex == indiv.getMedications().size() ? removeIndex-1 : removeIndex;
    newMed.erase(newMed.begin() + removeIndex);
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

Individual crossoverMutation(const Individual& indiv1, const Individual& indiv2,const Rcpp::DataFrame& ATCtree){
  int selectedNode, upperBound;
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  std::vector<int> newMedi{};
  newMedi.reserve(indiv1.getMedications().size() + indiv2.getMedications().size());
  //while we are on a leaf 
  do{
    selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
    selectedNode = selectedNode == ATCtree.nrow() ? ATCtree.nrow()-1 : selectedNode;
  } while (ATClength[selectedNode] == 7);
  
  upperBound = upperBounds[selectedNode];
  for(int med : indiv1.getMedications()){
    if(med <= selectedNode || med > upperBound){
      newMedi.push_back(med);
    }
  }
  for(int med : indiv2.getMedications()){
    if(med > selectedNode && med <= upperBound){
      newMedi.push_back(med);
    }
  }
  newMedi.shrink_to_fit();
  
  return {newMedi, indiv1.getTemperature()};
  
}

