#include "MCMC.h"

// [[Rcpp::depends(RcppArmadillo)]]


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

std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int nbMedic,num, itemp = 1;
  double temp;
  std::vector<int> medicVec(0);
  
  for(int i = 0; i < nbIndividuals; ++i){
    nbMedic = trunc(Rcpp::rnorm(1,3)[0]);
    // to avoid 0 medications patient
    nbMedic = nbMedic == 0 ? 1 : nbMedic;
    medicVec.reserve(nbMedic);
    
    for(int j = 0; j < nbMedic ; ++j){
      //Draw every medications for a patient
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
  int mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
  mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
  std::vector<int> newMed = indiv.getMedications();
  
  std::vector<int>::iterator find = std::find(newMed.begin(), newMed.end(), mutateMed);
  //the medication is in the vector so we erase it 
  if( find != newMed.end()){
    newMed.erase(find);
  } else{ // the medication is not in the vector so we add it
    newMed.push_back(mutateMed);
  }
  return {newMed,indiv.getTemperature()};
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
  
  upperBound = upperBounds[selectedNode] - 1;
  for(int med : indiv1.getMedications()){
    if(med < selectedNode || med > upperBound){
      newMedi.push_back(med);
    }
  }
  for(int med : indiv2.getMedications()){
    if(med >= selectedNode && med <= upperBound){
      newMedi.push_back(med);
    }
  }
  newMedi.shrink_to_fit();
  
  return {newMedi, indiv1.getTemperature()};
  
}

