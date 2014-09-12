ccmtest <-
function(Arm, Arsd, Brm, Brsd, iter=100, cuts=10) {
  ###Test output from the CCM function for significant causal link
  #Input are the mean and standard deviation output from multispatial CCM
  #for two processes (A and B). Arm and Brm are the means of rho and
  #Arsd and Brsd are the standard deviations of rho.
  #Iter describes how many iterations should be used in the nonparametric test
  #(100 in the paper), and "cuts" describes the number of bins that data are
  #put into (10 in the paper).
  
  Arm<-Arm[!is.na(Arm)]
  Arsd<-Arsd[!is.na(Arsd)]
  
  Brm<-Brm[!is.na(Brm)]
  Brsd<-Brsd[!is.na(Brsd)]
  
  HOA<-numeric(iter) #Hypothesis "O" - counts how frequently rho decreases with L
  HAA<-numeric(iter) #Hypothesis "A" - counts how frequently rho increases with L
  diffvalA<-0 #Counts number of times that decreases outnumber increases
  rhoA<-numeric(iter) #Finds mean value of bins for A to test that rho is significantly greater than 0
  
  HOB<-numeric(iter) #"O" distribution for B
  HAB<-numeric(iter) #"A" distribution for B
  diffvalB<-0
  rhoB<-numeric(iter)
  
  for(ii in 1:iter) {
    #A
    smp<-rnorm(length(Arm), Arm, Arsd) #Create a sample observation for A based on the rho and sd from the CCM output
    smp_cut<-unlist(tapply(smp, cut(1:length(smp), cuts), mean)) #Bin the data
    
    #For neighboring bins, test how frequently a bin with a shorter L has a higher rho than a bin with a shorter L
    HAA[ii]<-sum(smp_cut[-1]>smp_cut[-length(smp_cut)]) #Counts number of times that rho decreases with L
    HOA[ii]<-sum(smp_cut[-1]<smp_cut[-length(smp_cut)]) #Counts number of times that rho increases with L
    diffvalA<-diffvalA+as.numeric(HOA[ii]>HAA[ii]) #Count total number of times that rho decreases rather than increases
    rhoA[ii]<-sum(smp_cut<0)/cuts #Records frequency of bins for which mean rho is less than 0
    
    #B
    smp<-rnorm(length(Brm), Brm, Brsd)
    smp_cut<-unlist(tapply(smp, cut(1:length(smp), cuts), mean))
    
    HAB[ii]<-sum(smp_cut[-1]>smp_cut[-length(smp_cut)])
    HOB[ii]<-sum(smp_cut[-1]<smp_cut[-length(smp_cut)])
    diffvalB<-diffvalB+as.numeric(HOB[ii]>HAB[ii])
    rhoB[ii]<-sum(smp_cut<0)/cuts
    
  }
  
  #Returns:
  #slopeA, slopeB, the frequency where rho_A or rho_B decreases rather than increases with increasing L across bins, averaged across all iterations
  #rhoA, rhoB, the frequency where rho_A or rho_B is less than 0, averaged across all iterations 
  return(c(slopeA=diffvalA/iter, rhoA=mean(rhoA), slopeB=diffvalB/iter, rhoB=mean(rhoB)))
}
