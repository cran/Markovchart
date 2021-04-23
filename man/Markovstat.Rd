\name{Markovstat}
\alias{Markovstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stationary distribution calculation for processes with fixed/random shift size, random repair and random sampling time.
}
\description{
Calculates the stationary distribution of a process described by a discrete state, discrete time Markov chain. The process is described by a degradation-repair cycle type model. The user must give parameters describing both the degradation and the repair. The process is not repaired until the problem is discovered by sampling, hence the control chart setup. The same, single element is monitored (i.e. the sample size is always 1).
}
\usage{
Markovstat(shiftfun = c("exp", "exp-geo", "deg"), h, k, sigma,
           s, delta, probmix = 1, probnbin = 0.5, disj = 1,
           RanRep = FALSE, alpha = NULL, beta = NULL,
           RanSam = FALSE, StateDep = FALSE, a = NULL,
           b = NULL, q = NULL, z = NULL, Vd = 100, V,
           Qparam = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shiftfun}{
A string defining the shift size distribution to be used. Must be either \code{"exp"} (exponential), \code{"exp-geo"} (exponential-geometric mixture) or \code{"deg"} (degenerate). Use \code{"deg"} for fixed shift size with perfect repair and guaranteed sampling, i.e. Duncan"s traditional cycle model.
}
  \item{h}{
The time between samplings. Must be a positive value.
}
  \item{k}{
The control limit (critical value). Must be a positive value. Only one sided shifts are allowed, thus there is only one control limit.
}
  \item{sigma}{
Process standard deviation (the distribution is assumed to be normal).
}
  \item{s}{
Expected number of shifts in an unit time interval.
}
  \item{delta}{
Expected shift size. Used as the parameter of the exponential distribution (\code{shiftfun="exp"} or \code{"exp-geo"}), or simply as the size of the shift (\code{shiftfun="deg"}).
}
  \item{probmix}{
The weight of the geometric distribution in case of exponential-geometric mixture shift distribution; should be between 0 and 1.
}
  \item{probnbin}{
The probability parameter of the geometric distribution in case of exponential-geometric mixture shift distribution; should be between 0 and 1.
}
  \item{disj}{
The size of a discrete jump in case of exponential-geometric mixture shift distribution, must be a positive number.
}
  \item{RanRep}{
Logical. Should the repair be random? Default is \code{FALSE} (the repair is perfect, the process is always repaired to the target value). The repair is always perfect (non-random) for \code{shiftfun="deg"}.
}
  \item{alpha}{
First shape parameter for the random repair beta distribution.
}
  \item{beta}{
Second shape parameter for the random repair beta distribution.
}
  \item{RanSam}{
Logical. Should the sampling be random? Default is \code{FALSE} (no). The sampling is never random for \code{shiftfun="deg"}.
}
  \item{StateDep}{
Logical. Should the sampling probability also depend on the distance from the target value (state dependency)? (If TRUE, a beta distribution is used for the sampling probability, if \code{FALSE} then a logistic function.)
}
  \item{a}{
First parameter*\code{h} for the random sampling time beta distribution. The first shape parameter is \code{a/h} to create dependency on the time between samplings as described at the \code{StateDep} parameter.
}
  \item{b}{
Second shape parameter for the random sampling time beta distribution.
}
  \item{q}{
The steepness of the curve of the random sampling time logistic function.
}
  \item{z}{
The logistic sigmoid"s midpoint of the random sampling time logistic function.
}
  \item{Vd}{
Integer discretisation parameter: the number of states in the equidistant discretisation of the state space. Should be an integer value greater than 2.
}
  \item{V}{
Numeric discretisation parameter: the maximum (positive) distance from the target value taken into account.
}
  \item{Qparam}{
Integer discretisation parameter: the number of maximum events taken into account within a sampling interval.
}
}

\value{
The function return a list object of class \code{Markov_stationary}. The list is of length 3:
\item{Stationary_distribution }{Stationary distribution of the Markov chain. The probabilities in the stationary distribution are labeled. If \code{shiftfun} is "deg" then the stationary distribution is always of length 4. If \code{shiftfun} is not "deg" then there are multiple out-of-control and true alarm states. These are labeled with an index and the value the state represents. If \code{shiftfun} is "deg" then the out-of-control and true alarm states are at a distance \code{delta} from the target value, and the in-control and the false alarm state are always at the target value.}
\item{Transition_matrix }{The transition matrix of the Markov chain. Not printed.}
\item{Param_list }{Parameters given to the function and various technical results used by the \code{\link{Markovchart}} function. Not printed.}
}
\references{
Zempleni A, Veber M, Duarte B and Saraiva P. (2004) Control charts: a cost-optimization approach for processes with random shifts. \emph{Applied Stochastic Models in Business and Industry}, 20(3), 185-200.

Dobi B and Zempleni A. (2019) Markov chain-based cost-optimal
control charts for health care data. \emph{Quality and Reliability Engineering
International}, 35(5), 1379-1395.

Dobi B and Zempleni A. (2019) Markov chain-based cost-optimal
control charts with different shift size distributions. \emph{Annales Univ. Sci.
Budapest., Sect. Comp.}, 49, 129-146.}
\author{
Balazs Dobi and Andras Zempleni
}


\seealso{
\code{\link{Markovchart}}
}
\examples{
#Fixed shift size (essentially Duncan's cycle model).
res1 <- Markovstat(shiftfun="deg", h=1, k=1, sigma=1, s=0.2, delta=2.5)
res1

#Exponential shift - perfect repair - deterministic sampling
res2 <- Markovstat(shiftfun="exp", h=1, k=1, sigma=1, s=0.2, delta=2, Vd=30, V=18)
res2
#Notice how the In-control and the False-alarm states have non-zero probabilities.
#If the repair would be random (RanRep=TRUE), then these states would have zero probability.

#Exponential-geometric mixture shift - random repair - random sampling.
res3 <- Markovstat(shiftfun='exp-geo', h=1.5, k=2, sigma=1, s=0.2,
                   delta=1.2, probmix=0.7, probnbin=0.8, disj=2,
                   RanRep=TRUE, alpha=1, beta=3, RanSam=TRUE,
                   StateDep=TRUE, a=1, b=15, Vd=40, V=8)
res3
}
\keyword{Markov chain}
\keyword{control chart}
\keyword{stationary distribution}
