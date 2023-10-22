\name{make_ccm_data}
\alias{make_ccm_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes fake data for other functions
}
\description{
Builds a fake data set of two interacting processes, based on the model in the Sugihara et al. publication below, and based on a two-species discrete-time competition model. In the model, process A is causally affected by process B, but process B is not influenced by process A.
}
\usage{
make_ccm_data(sp_sd=0.125, obs_sd=0.025,
Sstr=0.375, times=10, burnin=100,
number_of_chains=20, seednum=2718)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp_sd}{
Standard deviation used to add process noise. If you are simulating multiple plots, this adds normally distributed noise with mean=0 to the growth rates of each species in different plots.
}
  \item{obs_sd}{
Standard deviation used to add observation error. Observation error is added to process X as a lognormal variable (X*rlnorm), with mean=0 on the lognormal scale.
}
  \item{Sstr}{
Forcing strength defining the effect of process B on process A.
}
  \item{times}{
Number of sequential observations desired for the time series in each plot (i.e., length of each independent time series)?
}
  \item{burnin}{
Burnin time before starting the experiment. This can be used to remove correlation among plots that occurs because of starting conditions.
}
  \item{number_of_chains}{
Total number of time series (i.e., how many replicates will be assembled into a single long time series?)
}
  \item{seednum}{
Random seed used for simulation.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Accm}{Time series for process A. Gaps between time series are indicated by a "NA" entry.}
\item{Bccm}{Time series for process B. Gaps between time series are indicated by a "NA" entry.}
\item{time_ccm}{Time indices corresponding to process A and B.}

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
CCM_boot, SSR_pred_boot, SSR_check_signal, ccmtest
}
\examples{
#Simulate data to use for multispatial CCM test
#See function for details - A is causally forced by B,
#but the reverse is not true.
ccm_data_out<-make_ccm_data()
Accm<-ccm_data_out$Accm
Bccm<-ccm_data_out$Bccm
}
\keyword{ts}
