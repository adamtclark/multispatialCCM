\name{multispatialCCM-package}
\alias{multispatialCCM-package}
\alias{multispatialCCM}
\docType{package}
\title{
multispatial convergent cross mapping analysis
}
\description{
A package for conducting convergent cross mapping (CCM) tests for causal relations for short, replicated time series. This package contains functions for implementing CCM on replicated data, as well as testing for significance of causal relationships that are detected. This package can also be used to implement state space reconstruction on replicated data (dew-drop regression, Hsieh et al.)
}
\details{
\tabular{ll}{
Package: \tab multispatialCCM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-19\cr
License: \tab GPL (>=2)\cr
}
The functions in this package can analyze data from multiple replicated time series. To enter the time series into functions, they should be concatenated back-to-back into a single vector, with gaps between time series indicated by "NA".

The package follows the descriptions in Clark et al. There are two classes of functions: analytical, and diagnostic. Both of these should be used for conducting multispatial CCM.

Analytical: CCM_boot is used to conduct multispatial CCM. SSR_pred_boot is used to find the optimal embedding dimension E for the CCM analysis, and can also be used to implement a state space reconstruction algorithm to predict process dynamics based on another process, or on its own historical dynamics.

Diagnostic: SSR_check_signal is used to test whether dynamics from a time series are sufficiently nonlinear for application of CCM, and tests that they are not too dominated by noise. The ccmtest function is used to test the output from boot_CCM to determine whether it indicates a significant causal link between variables.

Finally, the package also includes the function make_ccm_data, which can be used to make fake data from a simulation of two competing species. This function is used in the examples to show how to utilize the multispatial CCM functions.
}
\author{
Adam Clark

Maintainer: Adam Clark <atclark@umn.edu>
~~ The author and/or maintainer of the package ~~
}
\references{
Sugihara, G., R. May, H. Ye, C. Hsieh, E. Deyle, M. Fogarty, and S. Munch. 2012. Detecting Causality in Complex Ecosystems. Science 338.

Hsieh, C., C. Anderson, and G. Sugihara. 2008. Extending Nonlinear analysis to short ecological time series. American Naturalist 171:71–80.

Clark, A.T., Yi, H., Cowles, J., Deyle, E., Isbell, F., Sugihara, G., Tilman, D. 2014. Spatial 'convergent cross mapping' to detect causal relationships from short time-series. In review.
}
\seealso{
CCM_boot, SSR_pred_boot, SSR_check_signal, ccmtest
}
\examples{
#Simulate data to use for multispatial CCM test
#See function for details - A is causally forced by B,
#but the reverse is not true.
ccm_data_out<-make_ccm_data()
Accm<-ccm_data_out$Accm
Bccm<-ccm_data_out$Bccm

#Calculate optimal E
maxE<-5 #Maximum E to test
#Matrix for storing output
Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")

#Loop over potential E values and calculate predictive ability
#of each process for its own dynamics
for(E in 2:maxE) {
  #Uses defaults of looking forward one prediction step (predstep)
  #And using time lag intervals of one time step (tau)
  Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
  Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
}

#Look at plots to find E for each process at which
#predictive ability rho is maximized
matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
          xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("A", "B"), lty=1:2, col=1:2, lwd=2, bty="n")

#Results will vary depending on simulation.
#Using the seed we provide,
#maximum E for A should be 2, and maximum E for B should be 3.
#For the analyses in the paper, we use E=2 for all simulations.

E_A<-2
E_B<-3

#Check data for nonlinear signal that is not dominated by noise
#Checks whether predictive ability of processes declines with
#increasing time distance
#See manuscript and R code for details
signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
  predsteplist=1:10)
signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
  predsteplist=1:10)

#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
# Does A "cause" B?
#Note - increase iterations to 100 for consistant results
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=10)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=10)

#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
                    CCM_boot_B))

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))

#Plot "A causes B"
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
     xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
     xlab="L", ylab="rho")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)

#Plot "B causes A"
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)

legend("topleft",
c("A causes B", "B causes A"),
lty=c(1,2), col=c(1,2), lwd=2, bty="n")
}
\keyword{ts}
