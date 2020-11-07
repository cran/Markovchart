\name{Markovchart}
\alias{Markovchart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Economic design for the X-bar control chart with random shift size, random repair and random sampling time.
}
\description{
Wrapper for Markov chain-based cost optimal control charts. Includes methods for different shift size distributions and optimisation with respect to the average cost and cost standard deviation where the free parameters are the sampling interval (\code{h}) and the control limit (\code{k}).
}
\usage{
Markovchart(shiftfun = c("exp", "exp-geo", "deg"), h, k, sigma, s,
            delta, probmix = 0.5, probnbin = 0.5, disj = 1, RanRep
            = FALSE, alpha = NULL, beta = NULL, RanSam = FALSE,
            StateDep = FALSE, a = NULL, b = NULL, q = NULL,
            z = NULL, p = 1, Vd = 50, V, Qparam = 25,
            COST = c("no", "yes", "optim"), constantr = FALSE,
            ooc_rep = 0, cs = NULL, cofun = cofun_default,
            coparams = NULL, crfun = crfun_default, crparams = NULL,
            cf = crparams, vcofun = vcofun_default,
            vcoparams = c(0, 0), vcrfun = vcrfun_default,
            vcrparams = c(0, 0), method = c("L-BFGS-B", "Nelder-Mead",
            "BFGS", "CG", "SANN", "Brent"), parallel_opt = NULL,
            silent = TRUE, ...)
}
\arguments{
  \item{shiftfun}{
A string defining the shift size distribution to be used. Must be either \code{"exp"} (exponential), \code{"exp-geo"} (exponential-geometric mixture) or \code{"deg"} (degenerate). Use \code{"deg"} for fixed shift size with perfect repair and guaranteed sampling, i.e. Duncan"s traditional cycle model.
}
  \item{h}{
The time between samplings. Must be a positive value, can be a numeric vector. For optimisation, this is the initial value.
}
  \item{k}{
The control limit (critical value). Must be a positive value, can be a numeric vector. For optimisation, this is the initial value. Only one sided shifts are allowed, thus there is only one control limit.
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
  \item{p}{
The weight of the cost expectation in the calculation of the G-value; should be between 0 and 1.
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
  \item{COST}{
A string of value of either \code{"no"}, \code{"yes"} or \code{"optim"}. Default is no, in which case only the stationary distribution is returned. If yes, statistics related to costs are calculated and cost parameters must be provided. Optimisation (when the value is optim) is only possible with respect to costs.
}
  \item{constantr}{
Logical. Should the repair cost be assumed to constantly occur over time (\code{TRUE}) or assumed to only occur when there is a repair due to an alarm (\code{FALSE}, default)? If \code{TRUE}, then the repair cost should be given per unit time.
}
  \item{ooc_rep}{
Numeric value between 0 and 1. The percentage of repair cost ocurring during out-of-control operation. Default is 0. If a value greater than 0 is set, then \code{constantr} should be \code{TRUE}, but it is not forced.
}
  \item{cs}{
Sampling cost per sampling.
}
  \item{cofun}{
A function describing the relationship between the distance from the target value and the resulting out-of-control costs. Default is calculated using a base and a distance-scaling out-of-control parameter. See "Details".
}
  \item{coparams}{
Numeric vector. Parameters of \code{cofun}.
}
  \item{crfun}{
A function describing the relationship between the distance from the target value and the resulting repair costs. The default function assumes a linear relationship between the repair cost and the distance, and uses a base and a distance-scaling repair cost parameter. See "Details".
}
  \item{crparams}{
Numeric vector. Parameters of \code{crfun}.
}
  \item{cf}{
Numeric. The false alarm cost. Only relevant when \code{shiftfun} is \code{"deg"}.
}
  \item{vcofun}{
A function describing the relationship between the distance from the target value and the resulting out-of-control cost variance. For the default function see "Details".
}
  \item{vcoparams}{
Numeric vector. Parameters of \code{vcofun}.
}
  \item{vcrfun}{
A function describing the relationship between the distance from the target value and the resulting repair cost variance. For the default function see "Details".
}
  \item{vcrparams}{
Numeric vector. Parameters of \code{vcrfun}.
}
  \item{method}{
Method used for optimisation. Same as the argument of \code{optim}, but the default here is \code{"L-BFGS-B"}, because it turned out to be more robust in testing.
}
  \item{parallel_opt}{
A list of parallel options. See e.g. the argument \code{parallel} in the documentation of \code{\link[optimParallel]{optimParallel}}. Can be left empty, in this case the number of cores (threads) is automatically detected and all but one is used. (Single-core computers simply use one core.)
}
  \item{silent}{
Should the call be returned? Default is \code{FALSE}.
}
  \item{\dots}{
Further arguments to be passed down to \code{optimParallel}.
}
}
\details{
The \code{constantr} parameter is used for different repair assumptions. In traditional control chart theory, repair cost only occurs in case of an alarm signal. This is represented by \code{constantr=FALSE}, which is the default. In this case the repair is just a momentary cost, occurring at the time of the sampling. However this model is inappropriate in several cases in healthcare. For example there are chronic diseases that require constant medication (repair in the sense of the model). In this approach (\code{constantr=TRUE}) the repair cost still depends on the state of the process during sampling, but occurs even if there is no alarm and is divided by \code{h} to represent the constant repair through the whole sampling interval. Thus the repair cost should be given in a way which corresponds to the model used.

The default \code{cofun} calculates the out-of-control (OOC) cost using a base and a distance-scaling OOC parameter:

\deqn{c_{o} = c_{ob} + c_{os} A^2(v),}{c_o = c_ob + c_os A^2(v),}

where \eqn{c_{o}}{c_o} is the total OOC cost, \eqn{ c_{ob}}{c_ob} is the base OOC cost (even without shift), \eqn{c_{os}}{c_os} is the shift-scaling cost and \eqn{A^2(v)}{A^2(v)} is the squared distance from the target value. This latter part is defined like this because a Taguchi-type loss function is used. This \eqn{A^2(v)}{A^2(v)} incorporates the distances (the base of the losses) incurred not just at the time of the sampling, but also between samplings (hence it dependens on h). Even if the user defines a custom cost function for the OOC cost, this \eqn{A^2(v)}{A^2(v)} term must be included, as a closed form solution has been developed for the calculation of the squared distances in case of exponential shifts, considerably decreasing run times. Thus the arguments of the OOC cost function should look like this: function(\eqn{A^2(v)}{A^2(v)}, other parameters contained in a vector). \eqn{A^2(v)}{A^2(v)} is fed to the cost function as a vector, thus the function should vectorised with respect to this argument. The default function looks like this:

\preformatted{
cofun_default	<-	function(sqmudist,coparams)
{
  sqmudist=sqmudist
  cob=coparams[1]
  cos=coparams[2]
  co		<-	cob + cos*sqmudist
  return(co)
}
}

The default \code{vcofun} also uses a Taguchi-type loss function and has identical parts and requirements as \code{cofun}. The final standard deviation itself is calculated using the law of total variance. The default \code{vcofun} is:

\preformatted{
vcofun_default	<-	function(sqmudist,vcoparams)
{
  sqmudist=sqmudist
  vcob=vcoparams[1]
  vcos=vcoparams[2]
  vco		<-	vcob + vcos*sqmudist
  return(vco)
}
}

The defaults for the repair cost and cost variance are simple linear functions. For \code{crfun} it is

\deqn{c_{r} = c_{rb} + c_{rs} v,}{c_r = c_rb + c_rs v,}

where the notation are the same as before and "r" denotes repair. A custom function can be defined more freely here, but the first argument should be \eqn{v}{v} and the second a vector of further parameters.

The default function are:

\preformatted{
crfun_default	<-	function(mudist,crparams)
{
  mudist=mudist
  crb=crparams[1]
  crs=crparams[2]
  cr		<-	crb + crs*mudist
  return(cr)
}
}

\preformatted{
vcrfun_default	<-	function(mudist,vcrparams)
{
  mudist=mudist
  vcrb=vcrparams[1]
  vcrs=vcrparams[2]
  vcr		<-	vcrb + vcrs*mudist;
  return(vcr)
}
}

}
\value{
The value depends on the parameters:
If \code{h} and \code{k} are both of length 1, no costs are calculated (\code{COST}="no") and \code{shiftfun} is not "deg" then the value is the stationary distribution of the Markov chain which is a numeric vector of length \code{Vd}*2. The length is double of \code{Vd} because each state has an alarm and a non-alarm (out-of-control) version. If \code{shiftfun} is "deg" then the stationary distribution is always of length 4.
The probabilities in the stationary distribution are labeled. If \code{shiftfun} is not "deg" then there are multiple out-of-control and true alarm states. These are labeled with an index and the value the state represents. If \code{shiftfun} is "deg" then the out-of-control and true alarm states are at a distance \code{delta} from the target value. The in-control and the false alarm state are always at the target value.

If either \code{h} or \code{k} have length greater than 1, then the G-value (weighted average of average cost and cost standard deviation) is calculated for all given values without optimisation. The value of the function in this case is a data frame with \code{length(h)*length(k)} number of rows and three columns for \code{h}, \code{k} and the G-value.

If \code{h} and \code{k} are both of length 1 and costs are calculated (\code{COST}="yes" or "optim"), then the value of the function is a "Markov_chart" object, which is a list of length 4, detailing the properties of the control chart setup.
\item{Results }{Vector of G-value, expected cost, cost standard deviation and further process moments. Note that these further moments only take into account the process variation (i.e. the standard deviation), while the "Total cost std. dev." takes into account all sources of variance (e.g. the different costs that can occur due to being out-of-control).}
\item{Subcosts }{Vector of sub-costs that are parts of the total expected cost.}
\item{Parameters }{A vector that contains the time between samplings (\code{h}) and critical value (\code{k}) which was used in the control chart setup.}
\item{Stationary_distribution }{The stationary distribution of the Markov chain. See above.}

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

\examples{
#Defining parallel_opt parallel settings
#parallel_opt can also be left empty to be defined automatically by the function
require(parallel)
num_workers <-  min(c(detectCores(),2))
parall	    <-  list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)


#Fixed shift size (essentially Duncan's cycle model) - no optimisation
res <-  Markovchart(h=1, k=1, shiftfun="deg", sigma=1, s=0.2, delta=2.5, cs=1, crparams=20,
                    coparams=50, COST="yes")
res
\donttest{
#Fixed shift size (essentially Duncan's cycle model) - with optimisation
res <-  Markovchart(h=1, k=1, shiftfun="deg", sigma=1, s=0.2, delta=2.5, cs=1, crparams=20,
                    coparams=50, COST="optim", lower = c(0.01,0.01), upper = c(5,5),
                    parallel_opt=parall)
res

#Exponential shift - perfect repair - no optimisation - default cost functions - no cost calculation
res <-  Markovchart(h=1, k=1, shiftfun="exp", sigma=1,s=0.2, delta=2, RanRep=FALSE, COST="no",
		    cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
		    vcrparams=c(50000,-600000,1.5), p=0.9, Vd=30, V=18)
res
#Notice how the In-control and the False-alarm states have non-zero probabilities
#If the repair is random (RanRep=TRUE), then these states have zero probability
head(res)

#Exponential shift - no optimisation - default cost functions - no cost calculation
res <-  Markovchart(h=1, k=1, shiftfun="exp", sigma=1,s=0.2, delta=2, RanRep=TRUE, alpha=1, beta=3,
                    COST="no", cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
		    vcrparams=c(50000,-600000,1.5), p=0.9, Vd=30, V=18)
res

#Exponential shift - no optimisation - default cost functions
res <-  Markovchart(h=1, k=1, shiftfun="exp", sigma=1,s=0.2, delta=2, RanRep=TRUE, alpha=1, beta=3,
                    COST="yes", cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
		    vcrparams=c(50000,-600000,1.5), p=0.9, Vd=30, V=18)
res
}
#Exponential shift - with optimisation - default cost functions
parall	    <-  list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res <-  Markovchart(h=1, k=1, shiftfun="exp", sigma=1,s=0.2, delta=2, RanRep=TRUE, alpha=1, beta=3,
                    COST="optim", cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
		    vcrparams=c(50000,-600000,1.5), p=0.9, Vd=30, V=18, parallel_opt=parall)
res
\donttest{
#Exponential-geometric mixture shift - no optimisation -
#random sampling - custom repair variance function
vcrfun_new	<-	function(mudist,vcrparams)
{
  mudist=mudist
  vcrb=vcrparams[1]
  vcrs=vcrparams[2]
  vcrs2=vcrparams[3]

  vcr		<-	vcrb + vcrs/(mudist + vcrs2)
  return(vcr)
}

res <- Markovchart(h=1.5, k=2, shiftfun='exp-geo', sigma=1, s=0.2,
                   delta=1.2, probmix=0.7, probnbin=0.8, disj=2,
                   RanRep=TRUE, alpha=1, beta=3, RanSam=TRUE,
                   StateDep=TRUE, a=1, b=15, COST='yes', cs=1,
                   coparams=c(10,6), crparams=c(20,3),
                   vcoparams=c(10000,100), vcrfun=vcrfun_new,
                   vcrparams=c(50000,-600000,1.5), p=0.9,
                   Vd=100, V=8)
res

#Exponential shift - no optimisation  - vectorised
parall	    <-  list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
Gmtx	<-	Markovchart(h=seq(1,10,by=(10-1)/5),k=seq(0.1,5,by=(5-0.1)/5), shiftfun="exp",
                            sigma=1,s=0.2, delta=2, RanRep=TRUE, alpha=1, beta=3, COST="yes",
                            cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
                            vcrparams=c(50000,-600000,1.5), p=0.9, V=18, parallel_opt=parall)
Gmtx
}
}
\keyword{Markov chain}
\keyword{control chart}
