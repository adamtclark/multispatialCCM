\name{CCM_boot}
\alias{CCM_boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run multispatial CCM algorithm on two time series
}
\description{
Runs the multispatial convergent cross mapping algorithm on two time series, A and B, to determine whether process A is a forcing process (i.e., causally affects) or process B. A and B do not need to be from single, long time series, but rather can be combinations of many (e.g. spatially-replicated) time series. See "Arguments" for details.
}
\usage{
CCM_boot(A, B, E, tau=1,
DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2),
iterations=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Time series that is being tested as a forcing process (i.e., the causal process). Should be a single vector. If data come from multiple time series, gaps between these should be marked with an "NA". E.g., c(1,2,3, NA, 1,2,3) implies two time series, each of length 3. Order of plots does not matter (because they will be shuffled during bootstrapping), but should match the order used in B.
}
  \item{B}{
Time series that is being tested as a response process (i.e., the process being affected by the causal process). Should be a single vector. If data come from multiple time series, gaps between these should be marked with an "NA". E.g., c(1,2,3, NA, 1,2,3) implies two time series, each of length 3. Order of plots does not matter (because they will be shuffled during bootstrapping), but should match the order used in A.
}
  \item{E}{
Embedding dimension to use for the analysis. Should be based on dimension that provides the best prediction of process A against itself using function "SSR_pred_boot" (state space reconstruction).
}
  \item{tau}{
Number of time steps to use for lagged components in the attractor space. Defaults to 1.
}
  \item{DesiredL}{
Desired library lengths for which to compute CCM. Defaults to the maximum possible length ((tau * (E - 1) + (E + 1)):length(A) - E + 2) (though number of resulting predictions may be smaller because of gaps in the time series). Shortening this list (e.g., only predicting every nth element) will reduce run-time for the algorithm, but may also reduce ability to detect causal relations.
}
  \item{iterations}{
Number of iterations for bootstrapping. Defaults to 100.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{A}{Input process A}
\item{Aest}{Estimated values of A from B using CCM, for the last (longest) library length considered}
\item{B}{Input process B}
\item{rho}{Pearson correlation coefficient of estimates of A from B for each library length tested}
\item{varrho}{Variance of rho based on bootstrapping for each library length}
\item{sdevrho}{Standard error (i.e., sqrt(var/n)) of the mean of rho for each library length}
\item{Lobs}{Library lengths for which rho was calculated}
\item{E}{Embedding dimension used for the analysis}
\item{tau}{Time lag length used for the analysis}
\item{FULLinfo}{Includes all output from .C call. Much of this is redundant with the variables listed above. Best not to look at this.}
%% ...
}
\references{
Sugihara, G., R. May, H. Ye, C. Hsieh, E. Deyle, M. Fogarty, and S. Munch. 2012. Detecting Causality in Complex Ecosystems. Science 338.

Adam T. Clark, H. Ye, Forest Isbell, Ethan R. Deyle, Jane Cowles, David Tilman, and George Sugihara. 2015. Spatial ’convergent cross mapping’ to detect causal relationships from short time-series. Ecology, 96(6):1174–1181.
}
\author{
Adam Clark
}
\section{Warning}{
If you do not separate distinct time series with "NA" as described in "Arguments", they will not be treated as such!
}

\seealso{
SSR_pred_boot, SSR_check_signal, ccmtest
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
}
\keyword{ts}
