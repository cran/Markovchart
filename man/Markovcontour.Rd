\name{Markovcontour}
\alias{Markovcontour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Contour plot for control chart results.
}
\description{
Convenience function for plotting G-values in a contour plot as the function of the time between samplings and the critical value.
}
\usage{
Markovcontour(Gmtx, xlab = "Time between samplings", ylab = "Critical value",
              low = "white", mid = "#999999", high = "black", colour = "white",
              name = expression(atop(italic("G")*-value~per, unit~time)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gmtx}{
A data.frame with three columns (preferably created by the Markovchart function): time between samplings, critical value and the weighted mean of the expected cost and the cost standard deviation (G-values).
}
  \item{xlab}{
A title for the x axis.
}
  \item{ylab}{
A title for the x axis.
}
  \item{low}{
Colour for the low end of the gradient.
}
  \item{mid}{
Colour for the midpoint.
}
  \item{high}{
Colour for the high end of the gradient.
}
  \item{colour}{
Colour of the contour lines.
}
  \item{name}{
The name of the scale.
}
}
\value{
A plot object of class "gg" and "ggplot" produced using the \code{ggplot2} package.
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
\note{
The plot itself is made using the package \code{ggplot2} by Hadley Wickham et al. The text on the contour lines is added with the \code{geom_text_contour} function from the package \code{metR} by Elio Campitelli.
}

\seealso{
\code{\link{Markovchart}}
}

\examples{
\donttest{
#Defining parallel_opt parallel settings
#parallel_opt can also be left empty to be defined automatically by the function
require(parallel)
num_workers <- min(c(detectCores(),2))

parall	<-	list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
Gmtx	<-	Markovchart(h=seq(1,10,by=(10-1)/5),k=seq(0.1,5,by=(5-0.1)/5), shiftfun="exp",
                            sigma=1,s=0.2, delta=2, RanRep=TRUE, alpha=1, beta=3, COST="yes",
                            cs=1, coparams=c(10,3), crparams=c(1,2), vcoparams=c(8000,100),
                            vcrparams=c(50000,-600000,1.5), p=0.9, V=18, parallel_opt=parall)
Markovcontour(Gmtx)
}
}
\keyword{Markov chain}
\keyword{control chart}
\keyword{contour}
\keyword{plot}
