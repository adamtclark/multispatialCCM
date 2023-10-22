\name{SSR_check_signal}
\alias{SSR_check_signal}
\title{
Test process for auto-predictability.
}
\description{
Predict elements of a process based historical observations of that process using cross-validation. Tests whether past observations are able to make good estimates of future elements of the time series.}
\usage{
SSR_check_signal(A, E, tau = 1,
predsteplist = 1:10, matchSugi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Process to be predicted. Should be a single vector. If data come from multiple time series, gaps between these should be marked with an "NA".
}
  \item{E}{
Embedding dimension to use for the analysis. Should be based on dimension that provides the best prediction of process A against itself using function "SSR_pred_boot" (state space reconstruction).
}
  \item{tau}{
Number of time steps to use for lagged components in the attractor space. Defaults to 1.
}
  \item{predsteplist}{
Vector of time step lengths for prediction.
}
  \item{matchSugi}{
Set to 1 to match results in Sugihara et al. publication described below, which removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
predatout=predatout,
rho_pre_slope=rho_pre_slope,
rho_predmaxCI=rho_predmaxCI

\item{predatout}{Vector of rho values describing predictive ability of process against itself for each prediction time step length}
\item{rho_pre_slope}{Slope of rho values as a function of prediction distance}
\item{rho_predmaxCI}{95\% confidence interval for rho value corresponding to the longest prediction interval tested}
%% ...
}
\references{
Sugihara, G., R. May, H. Ye, C. Hsieh, E. Deyle, M. Fogarty, and S. Munch. 2012. Detecting Causality in Complex Ecosystems. Science 338.

Adam T. Clark, H. Ye, Forest Isbell, Ethan R. Deyle, Jane Cowles, David Tilman, and George Sugihara. 2015. Spatial ’convergent cross mapping’ to detect causal relationships from short time-series. Ecology, 96(6):1174–1181.
}
\author{
Adam Clark
}
\seealso{
CCM_boot, SSR_pred_boot, ccmtest
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

#Check data for nonlinear signal that is not dominated by noise
#Checks whether predictive ability of processes declines with
#increasing time distance
#See manuscript and R code for details
signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
  predsteplist=1:10)
signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
  predsteplist=1:10)
}
\keyword{ts}
