\name{ccmtest}
\alias{ccmtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for significant causal signal
}
\description{
Tests output from CCM_boot for significant causal signal. This function compares the 95\% confidence intervals for esimated rho for the shortest and longest libraries calculated, and uses this to determine whether predictive power has significantly increased.
}
\usage{
ccmtest(CCM_boot_A, CCM_boot_B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CCM_boot_A}{
Output structure from a CCM test using CCM_boot
}
  \item{CCM_boot_B}{
Output structure from a CCM test using CCM_boot
}
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
\item{res}{Structure containing the p-values for both tests.}

%% ...
}
\references{
Sugihara, G., R. May, H. Ye, C. Hsieh, E. Deyle, M. Fogarty, and S. Munch. 2012. Detecting Causality in Complex Ecosystems. Science 338.

Clark, A.T., Yi, H., Cowles, J., Deyle, E., Isbell, F., Sugihara, G., Tilman, D. 2014. Spatial 'convergent cross mapping' to detect causal relationships from short time-series. In review.
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
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=10)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=10)

(CCM_significance_test<-ccmtest(CCM_boot_A,
                   CCM_boot_B))
}
\keyword{ts}
