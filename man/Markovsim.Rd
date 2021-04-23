\name{Markovsim}
\alias{Markovsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Progression and monitoring simulation of a process with random shift size, random repair and random sampling time.
}
\description{
Wrapper for simulation of processes with a Markov chain-based control chart setup. Includes methods for different shift size distributions.
}
\usage{
Markovsim(shiftfun = c("exp", "exp-geo"), num = 100, h, k, sigma,
          s, delta, probmix = 1, probnbin = 0.5, disj=1,
          RanRep = FALSE, alpha = NULL, beta = NULL, RanSam = FALSE,
          StateDep = FALSE, a = NULL, b = NULL, q = NULL,
          z = NULL, detail = 100, Vd = 50, V, burnin = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shiftfun}{
A string defining the shift size distribution to be used. Must be either "exp", "exp-geo".
}
  \item{num}{
Integer. The number of sampling intervals simulated. This means that the time elapsed in the simulation is \code{num*h}.
}
  \item{h}{
The time between samplings. Must be a positive value.
}
  \item{k}{
The control limit (critical value). Must be a positive value. Only one sided shifts are allowed, thus there is only one control limit.
}
  \item{sigma}{
Process standard deviation (the distribution is normal).
}
  \item{s}{
Expected number of shifts in an unit time interval.
}
  \item{delta}{
Expected shift size.
}
  \item{probmix}{
The weight of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1.
}
  \item{probnbin}{
The probability parameter of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1.
}
  \item{disj}{
The size of a discrete jump in case of exponential-geometric mixture shift distribution, must be a positive number.
}
  \item{RanRep}{
Logical. Should the repair be random? Default is FALSE (no).
}
  \item{alpha}{
First shape parameter for the random repair beta distribution.
}
  \item{beta}{
Second shape parameter for the random repair beta distribution.
}
  \item{RanSam}{
Logical. Should the sampling be random? Default is FALSE (no).
}
  \item{StateDep}{
Logical. Should the sampling probability also depend on the distance from the target value (state dependency)? (If TRUE, a beta distribution is used for the sampling probability, if FALSE then a logistic function.)
}
  \item{a}{
First parameter*\code{h}  for the random sampling time beta distribution. The first shape parameter is \code{a/h} to create dependency on the time between samplings.
}
  \item{b}{
Second shape parameter for the random sampling time beta distribution.
}
  \item{q}{
The steepness of the curve of the random sampling time logistic function.
}
  \item{z}{
The logistic sigmoid's midpoint of the random sampling time logistic function.
}
  \item{detail}{
The detail of the simulation, i.e. how many data points should be simulated within a sampling interval (including the sampling itself). Should be a positive integer greater than 1.
}
  \item{Vd}{
Integer discretisation parameter: the number of states after the equidistant discretisation of the state space. Should be an integer value greater than 2. This parameter is needed to calculate a stationary distibution that can be compared to results of the \code{Markovchart} function.
}
  \item{V}{
Numeric discretisation parameter: the maximum (positive) distance from the target value taken into account. This parameter is needed to calculate a stationary distibution that can be compared to results of the \code{Markovchart} function and for the calculation of sampling probabilities in the case of random sampling.
}
  \item{burnin}{
Numeric burn-in parameter: the number of samplings deemed as a burn-in period. Should be an integer greater than one.
}

}
\details{
The simulation only includes the more complicated process and control chart cases and is meant for model checking and for situations when the exact calculation is problematic (such as low probabilities in the stationary distribution leading to rounding errors).
}
\value{
A \code{Markov_sim} object which is a list of length 4.
\item{Value_at_samplings }{The process value at sampling.}
\item{Sampling_event }{The event at sampling, each can either be success (there was a sampling but no alarm), alarm (sampling with alarm) or failure (no sampling occurred).}
\item{Simulation_data }{The simulated data (distances from the target value).}
\item{Stationary_distribution }{The stationary distribution of the Markov chain, created by discretising the simulated data. See the documentaion of the \code{\link[Markovchart]{Markovchart}} function.}
}
\references{
Zempleni A, Veber M, Duarte B and Saraiva P. (2004) Control charts: a cost-optimization approach for processes with random shifts. \emph{Applied Stochastic Models in Business and Industry}, 20(3), 185-200.

Dobi B and Zempleni A. (2019) Markov chain-based cost-optimal
control charts for health care data. \emph{Quality and Reliability Engineering
International}, 35(5), 1379-1395.

Dobi B and Zempleni A. (2019) Markov chain-based cost-optimal
control charts with different shift size distributions. \emph{Annales Univ. Sci.
Budapest., Sect. Comp.}, 49, 129-146.
}
\author{
Balazs Dobi and Andras Zempleni
}

\seealso{
\code{\link{Markovstat}}
}

\examples{
#Simulation using exponential shifts, random repair and random samling.
simres1	<-	Markovsim(shiftfun="exp", num=500, h=1, k=1, sigma=1, s=0.2, delta=2,
                          RanRep=TRUE, alpha=1, beta=3, RanSam=TRUE, StateDep=TRUE,
                          a=0.1, b=1, V=10)
simres1
hist(simres1[[1]], 20, freq=FALSE)

#Simulation using exponential-geometric mixture shifts, random repair and random samling.
simres2	<-	Markovsim(shiftfun="exp-geo", num=500, h=1, k=1, sigma=1, s=0.2, delta=2,
                          probmix=0.9, probnbin=0.6, RanRep=TRUE, alpha=1, beta=3, RanSam=TRUE,
                          StateDep=TRUE, a=0.1, b=1, V=10)
simres2
hist(simres2[[1]], 20, freq=FALSE)
}
\keyword{Markov chain}
\keyword{control chart}
\keyword{simulation}
