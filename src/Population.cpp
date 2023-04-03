#include "Population.h"

// [[Rcpp::depends(RcppArmadillo)]]

Population::Population(int nbIndividuals) : individuals_{}
{
  individuals_.reserve(nbIndividuals);
}

Population::Population(int treeSize, int nbIndividuals,double meanMedic) : individuals_{}{
  std::vector<Individual> individuals = DFtoCPP_WOIndividual(treeSize, nbIndividuals, meanMedic);
  std::vector<double> fitnessVector(nbIndividuals);
  
  individuals_.reserve(individuals.size());
  for(int i = 0 ; i < individuals.size(); ++i){
    individuals_.push_back({fitnessVector[i], individuals[i]});
  }
}


Population::Population(const Population& pop) : individuals_{pop.individuals_}
{}

void Population::addAnIndividualToPopulation(const std::pair<double,Individual>& ind){
  individuals_.push_back(ind);
}

void Population::evaluate(const Rcpp::List& medications, const Rcpp::LogicalVector& ADR
                , const Rcpp::DataFrame& ATCtree){
  for(auto& indiv : individuals_){
    indiv.first = indiv.second.computeRR(medications, ADR, ATCtree);
  }
}

void Population::keepElite(int nbElite, Population& matingPool) const{
  if (nbElite > 0){
    std::priority_queue<std::pair<double,Individual>, std::vector<std::pair<double,Individual>>,
                        std::greater<std::pair<double,Individual>>> pq;
    
    for (const auto & i : individuals_) {
      if (pq.size() < nbElite) {
        pq.push(i);
      } else if (i > pq.top()) {
        pq.pop();
        pq.push(i);
      }
    }
    
    while (!pq.empty()) {
      matingPool.addAnIndividualToPopulation(pq.top());
      pq.pop();
    }
  }
}

void Population::tournamentSelection(int tournamentSize, Population& matingPool, int nbDrawing) const{
  int popSize = individuals_.size();
  int indexDrawing;
  double indexMax;
  std::vector<int> selectedIndividuals;
  selectedIndividuals.reserve(tournamentSize);
  //we have to add nbDrawing in the mating pool
  for(int i = 0 ; i < nbDrawing; ++i){
    //select the index of the players of the tournament
    selectedIndividuals.clear();
    for(int j = 0; j < tournamentSize ; ++j){
      do{
        indexDrawing = Rcpp::runif(1, 0, popSize)[0];
      } while (std::find(selectedIndividuals.begin(), selectedIndividuals.end(), indexDrawing) != selectedIndividuals.end());
      selectedIndividuals.push_back(indexDrawing);
    }
    //now make them compete again each other
    indexMax = selectedIndividuals[0];
    for(int j = 1; j < tournamentSize ; ++j){
      indexMax = individuals_[indexMax] > individuals_[selectedIndividuals[j]] ? indexMax : selectedIndividuals[j];
    }
    matingPool.addAnIndividualToPopulation(individuals_[indexMax]);
  }
}

void Population::crossover(int nbElite, const std::vector<int>& ATClength, const std::vector<int>& upperBounds,
                           const Rcpp::DataFrame& ATCtree, double p_crossover){
  int remainingIndividuals = individuals_.size() - nbElite;
  int selectedNode;
  int upperBound;
  Individual tmp;
  double draw;
  

  for(int i = nbElite ; i < individuals_.size()-2; i+=2){
    draw = Rcpp::runif(1,0,1)[0];
    // With probability p_crossover we opperate a crossover, otherwise we let the individuals as it is
    if(draw <= p_crossover){
      do{
        selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
      } while (ATClength[selectedNode] == 7);
      
      upperBound = upperBounds[selectedNode];
      
      tmp = individuals_[i].second;
      individuals_[i].second = crossoverMutation(individuals_[i].second, individuals_[i+1].second, ATCtree,
                                                 selectedNode, upperBound);
      individuals_[i+1].second = crossoverMutation(individuals_[i+1].second, tmp, ATCtree,
                                                   selectedNode, upperBound);
    }
  }
  //if the number of remaining individuals is even, we have to operate one more crossover on the last 2 individuals
  //otherwise we just have one individual so we just copy it (equivalent of doing nothing)
  if(remainingIndividuals %2 == 0){
    draw = Rcpp::runif(1,0,1)[0];
    if(draw <= p_crossover){
      do{
        selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
      } while (ATClength[selectedNode] == 7);
      upperBound = upperBounds[selectedNode];
      
      tmp = individuals_[individuals_.size()-2].second;
      individuals_[individuals_.size()-2].second = crossoverMutation(individuals_[individuals_.size()-2].second, 
                                                                    individuals_[individuals_.size()-1].second, ATCtree,
                                                                    selectedNode, upperBound);
      individuals_[individuals_.size()-1].second = crossoverMutation(individuals_[individuals_.size()-1].second, tmp,
                                                                     ATCtree, selectedNode, upperBound);
    }
  }
}

void Population::mutate(int nbElite, double p_mutation, const Rcpp::DataFrame& ATCtree){
  double draw, drawMutation, chosenVertexIdx;
  double alpha = 1; // TODO : put this as a parameter
  bool emptyCocktail;
  std::vector<std::pair<int,int>> vertex;
  
  for(int i = nbElite ;  i < individuals_.size() ; ++i){
    draw = Rcpp::runif(1,0,1)[0];
    //for each individual in the population we draw a number uniformly in [0;1]
    // with probability p_mutation we mutate the individual i
    if(draw <= p_mutation){
      //we mutate with type 1 mutation or type 2 equiprobably
      emptyCocktail = (individuals_[i].second.getMedications().size() == 0);
      drawMutation = Rcpp::runif(1,0,1)[0];
      if(drawMutation <= 0.5){
        //type1 mutation
        individuals_[i].second = type1Mutation(individuals_[i].second, ATCtree.nrow(), alpha, emptyCocktail);
      }
      else{
        //type2 mutation
        if(!emptyCocktail){
          vertex = individuals_[i].second.getVertexList(ATCtree);
          chosenVertexIdx = trunc(Rcpp::runif(1,0,vertex.size())[0]);
          individuals_[i].second = type2Mutation(individuals_[i].second, ATCtree.nrow(), vertex[chosenVertexIdx]);
        
        }
      }
      
    }
  }
}

void Population::clear(){
  individuals_.clear();
}

double Population::getMean() const{
  const auto size = static_cast<double>(individuals_.size());
  return std::accumulate(individuals_.begin(), individuals_.end(),
                         0.0, [](double& a, const std::pair<double,Individual>& b){
                           return a + b.first;
                           }) / size;
}

int Population::bestIndividual() const{
  int i_max = 0;
  for(int i = 1; i < individuals_.size() ; ++i){
    i_max = individuals_[i].first > individuals_[i_max].first ? i : i_max;
  }
  return i_max;
}

void Population::printPopulation(std::ostream& ost) const{
  for(const auto& indiv : individuals_){
    ost << "RR : " << indiv.first << "\n medication : ";
    indiv.second.printMedications();
  }
}

void Population::printSummary(int epoch, double populationMean, int populationBestIndex) const{
  std::cout << "epoch : " << epoch << " | mean : " << populationMean << " | best RR : ";
  std::cout << individuals_[populationBestIndex].first << " | best cocktail : ";
  individuals_[populationBestIndex].second.printMedications();
}

Population& Population::operator=(const Population& pop){
  individuals_ = pop.individuals_;
  return *this;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop){
  pop.printPopulation(ost);
  return ost;
}

