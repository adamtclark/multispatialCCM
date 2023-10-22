\name{SSR_pred_boot}
\alias{SSR_pred_boot}
\title{
State space reconstruction function
}
\description{
Predict elements of A using B using state space reconstruction. If A=B, then the algorithm uses cross validation to assess the ability of historical portions of the A time series to predict future components of the time series. This function can be used to find the embedding dimension E that maximizes predictive ability.}
\usage{
SSR_pred_boot(A, B = A, E, tau = 1, predstep = 1, matchSugi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Process to be compared to B, or to itself. Should be a single vector. If data come from multiple time series, gaps between these should be marked with an "NA".
}
  \item{B}{
Process to be compared to A. If left empty, algorithm defaults to A=B.
}
  \item{E}{
Embedding dimension to use for the analysis. Should be based on dimension that provides the best prediction of process A against itself using function "SSR_pred_boot" (state space reconstruction).
}
  \item{tau}{
Number of time steps to use for lagged components in the attractor space. Defaults to 1.
}
  \item{predstep}{
Number of time steps into the future to make predictions from past observations.
}
  \item{matchSugi}{
Set to 1 to match results in Sugihara et al. publication described below, which removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{A}{Returns variable from input}
\item{Aest}{Estimated values for A}
\item{B}{Returns variable from input}
\item{E}{Returns variable from input}
\item{tau}{Returns variable from input}
\item{pAlength}{Length of A from input}
\item{pBlength}{Length of B from input}
\item{predstep}{Returns variable from input}
\item{prepvec}{Returns 1 if A and B were treated as same process}
\item{pmatchSugi}{Returns variable from input}
\item{acceptablelib}{List of library lengths that were used for the analysis, adjusting for ends and gaps in the library}
\item{plengthacceptablelib}{Length of library that was used for the analysis}
\item{rho}{Pearson correlation coefficient describing predictive ability of A against B or against itself}
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
CCM_boot, SSR_check_signal, ccmtest
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
}
\keyword{ts}
