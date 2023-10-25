medic = c()
med = c()
for(i in 1:nrow(ATC_Tree_UpperBound_2014)){
  if(ATC_Tree_UpperBound_2014[i, "ATC_length"] == 7){
    med <- append(med, (i-1))
  }
}
#bounds of the drug set we want to try (minus 1 because the data frame index of the tree will be shifted by
#1 because the algorithms are written in C++)
firstLowerBound <- 48
firstUpperBound <- 53
secondLowerBound <- 1394
secondUpperBound <- 1396

simulPatient <- list()
simulPatient$patientADR <- vector("logical",10000)
simulPatient$patientATC <- vector("list",10000)

med_chance <- .20 #probability to have the cocktail 
ADR_chance_med <- .80 # probability to have the ADR given that the cocktail was taken
ADR_chance_no_med <- .11 # probability to have the ADR given that the cocktail was not taken

# for 1 to the number of patient we want
for(i in 1:10000){
  #number of drugs drawn for a patient 
  nbMedic <- 1 + rpois(1,1)
  chance <- runif(1,0,1)
  # with a probability of 20% the current patient will take a medication under the bounds  
  if(chance <= med_chance){
    medic = append(medic,sample(firstLowerBound:firstUpperBound,1,replace = T))
    medic = append(medic,sample(secondLowerBound:secondUpperBound,1,replace = T))
    #given that the patient takes a medication under the bounds his probability to get the ADR will be 80%
    if(runif(1) < (1-ADR_chance_med)){
      simulPatient$patientADR[i] = F
    }
    else{
      simulPatient$patientADR[i] = T
    }
  }# otherwise we choose the drugs randomly
  else{
    for(j in 1:nbMedic){
      medic = append(medic,sample(med,1,replace = TRUE))
    }# given that the patient takes random drugs his probability to get the ADR will be 11% (may be changed) 
    if(runif(1) > ADR_chance_no_med){
      simulPatient$patientADR[i] = F
    }else{
      simulPatient$patientADR[i] = T
    }
  }
  #put the current drugs list into the patient list then reset the current drugs list
  simulPatient$patientATC[[i]] = c(medic)
  medic = c()
}
#transform into a dataFrame -> ! THIS is the dataFrame to use in the algorithm !
simulPatient_df_2sizeCocktail = data.frame(patientATC = c(1:10000),patientADR = c(FALSE))
simulPatient_df_2sizeCocktail$patientATC = c(simulPatient$patientATC)
simulPatient_df_2sizeCocktail$patientADR = simulPatient$patientADR
