medic = c()
med = c()
for(i in 1:nrow(ATC_Tree_UpperBound_2014)){
  if(ATC_Tree_UpperBound_2014[i, "ATC_length"] == 7){
    med <- append(med, (i-1))
  }
}
#bounds of the drug set we want to try (minus 1 because the data frame index of the tree will be shifted by
#1 because the algorithms are written in C++)
lowerBound <- 47
upperBound <- 53
simulPatient <- list()
simulPatient$patientADR <- vector("logical",10000)
simulPatient$patientATC <- vector("list",10000)
# for 1 to the number of patient we want
for(i in 1:10000){
  #number of drugs drawn for a patient 
  nbMedic <- 1 + rpois(1,1)
  chance <- runif(1,0,1)
  # with a probability of 20% the current patient will take a medication under the bounds  
  if(chance <= 0.20){
    medic = append(medic,sample(lowerBound:upperBound,1,replace = T))
    #given that the patient takes a medication under the bounds he's probability to get the ADR will be 80%
    if(runif(1) < 0.2){
      simulPatient$patientADR[i] = F
    }
    else{
      simulPatient$patientADR[i] = T
    }
  }# otherwise we choose the drugs randomly
  else{
    for(j in 1:nbMedic){
      medic = append(medic,sample(med,1,replace = TRUE))
    }# given that the patient takes random drugs he's probability to get the ADR will be 20% (may be changed) 
    if(runif(1) > 0.05){
      simulPatient$patientADR[i] = F
    }else{
      simulPatient$patientADR[i] = T
    }
  }
  #put the drugs list into the list then reset the drugs list
  simulPatient$patientATC[[i]] = c(medic)
  medic = c()
}
#transform into a dataFrame
simulPatient_df_2 = data.frame(patientATC = c(1:10000),patientADR = c(FALSE))
simulPatient_df_2$patientATC = c(simulPatient$patientATC)
simulPatient_df_2$patientADR = simulPatient$patientADR
