\name{ccmtest}
\alias{ccmtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for significant causal signal
}
\description{
Tests output from CCM_boot for significant causal signal. This function bins the mean and standard errors into discrete categories, and then simulates values from them to determine whether rho significantly increases with library length.  The function requires two sets of inputs - rho mean and standard error estimates for both process A as predicted by process B (which tests whether B is causally influenced by A), and B as predicted by A (which tests the reverse).
}
\usage{
ccmtest(Arm, Arsd, Brm, Brsd, iter=100, cuts=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Arm}{
Mean estimates from rho based on output from CCM_boot for process A as predicted by process B
}
  \item{Arsd}{
Standard error estimates from rho based on output from CCM_boot for process A as predicted by process B
}
  \item{Brm}{
Mean estimates from rho based on output from CCM_boot for process B as predicted by process A
}
  \item{Brsd}{
Standard error estimates from rho based on output from CCM_boot for process B as predicted by process A
}
  \item{iter}{
Number of random draws from the binned distribution that will be used to determine significant increase of rho with library length
}
  \item{cuts}{
Number of bins for the test
}
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
\item{slopeA}{The frequency with which rho of process A decreases rather than increases with increasing L across bins, averaged across all iterations. A value of <0.05 indicates significant increase of rho with library length.}
\item{slopeB}{The frequency with which rho of process B decreases rather than increases with increasing L across bins, averaged across all iterations. A value of <0.05 indicates significant increase of rho with library length.}
\item{rhoA}{The frequency where rho of process A is less than 0, averaged across all iterations. A value of <0.05 indicates that rho is significantly greater than zero.}
\item{rhoB}{The frequency where rho of process B is less than 0, averaged across all iterations. A value of <0.05 indicates that rho is significantly greater than zero.}

%% ...
}
\references{
Sugihara, G., R. May, H. Ye, C. Hsieh, E. Deyle, M. Fogarty, and S. Munch. 2012. Detecting Causality in Complex Ecosystems. Science 338.

Clark, A.T., Yi, H., Cowles, J., Deyle, E., Isbell, F., Sugihara, G., Tilman, D. 2014. Spatial 'convergent cross mapping' to detect causal relationships from short time-series. In prep.
}
\author{
Adam Clark
}
\seealso{
CCM_boot, SSR_pred_boot, SSR_check_signal
}
\examples{
#Simulate data to use for multispatial CCM test
#See function for details - A is causally forced by B,
#but the reverse is not true.
ccm_data_out<-make_ccm_data()
Accm<-ccm_data_out$Accm
Bccm<-ccm_data_out$Bccm

#Set optimal E - see multispatialCCM for details
E_A<-2
E_B<-3

#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
# Does A "cause" B?
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=100)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=100)

(CCM_significance_test<-ccmtest(Arm=CCM_boot_A$rho,
                    Arsd=CCM_boot_A$sdevrho,
                    Brm=CCM_boot_B$rho,
                    Brsd=CCM_boot_B$sdevrho,
                    iter=100, cuts=10))
#Significant causal relationship is detected if both slope and rho
#tests show p<0.05.
#Based on our starting random number seed, the test should
#show a non-significant forcing of B by A (slopeA==0.330, rhoA==0.001)
#and show a significant forcing of A by B (slopeB==0.000, rhoB==0.000)
}
\keyword{ts}
