\name{LDL}
\alias{LDL}
\docType{data}
\title{
Aggregated low-density-lipoprotein patient data for control chart applications
}
\description{
  A data frame  containing aggregated low-density-lipoprotein (LDL) patient data gathered from various sources.
}
\usage{data("LDL")}
\format{
  A data frame with 1 observation on the following 12 variables.
  \describe{
    \item{\code{target_value}}{Target LDL value.  Set according to the European guideline for patients at risk. Garmendia et al. (2000)}
    \item{\code{standard_deviation}}{Process standard deviation. Estimated using real life data from Hungary, namely registry data through Healthware Consulting Ltd. (Budapest, Hungary).}
    \item{\code{expected_shift_size}}{Expected shift size, around 0.8 increase in LDL per year on average, due to the expected number of shifts per year. Estimated with the help of a health care professional from the Healthware Consulting Ltd.}
    \item{\code{num_exp_daily_shifts}}{We expect around 3 shifts per year on average. Estimated with the help of a health professional from the Healthware Consulting Ltd.}
    \item{\code{rep_size_first}}{First repair size distribution parameter. Estimated using an international study which included Hungary. Garmendia et al. (2000)}
    \item{\code{rep_size_second}}{Second repair size distribution parameter.}
    \item{\code{samp_prob_first}}{First sampling probability parameter. Patient non-compliance in LDL controlling medicine is quite high, and this is represented through the parametrisation of the logistic function. Lardizabal
and Deedwania (2010)}
    \item{\code{samp_prob_second}}{Second sampling probability parameter.}
    \item{\code{sampling_cost}}{Sampling cost in Euro. Estimated using the official LDL testing cost and visit cost in Hungary.}
    \item{\code{OOC_cost}}{Out-of-control operation (health care event) cost in Euro. Estimated using real world data of cardiovascular event costs from Hungary}
    \item{\code{base_rep_cost}}{Base repair (treatment) cost in Euro. Estimated using the simvastatin therapy costs in Hungary}
    \item{\code{prop_rep_cost}}{Shift-proportional (illness-proportional) repair cost in Euro. Estimated using the simvastatin therapy costs in Hungary.}
  }
}
\details{
It is very difficult to give a good estimate for the type and parameters of a distribution that properly models the non-compliance (sampling probability), thus the data here can at best be regarded as close approximations to a real-life situation. This is not a limiting factor, as patients themselves can have vast differences in their behaviour, so evaluation of different scenarios are often required. Since high LDL levels rarely produce noticeable symptoms, the sampling probability should only depend on the time between samplings (Institute for Quality and Efficiency in Health Care, 2017). Thus, the sampling probability parameters assume the use of a logistic function and not a beta distribution in the \code{Markovstat} function. It is important to note that the proportional costs usually assumed to increase according to a Taguchi-type loss function (squared loss), thus huge expenses can be generated if the patient’s health is highly out-of-control. For cost calculations, the 2021 March 21 EUR-HUF exchange rate was used (1 EUR = 369.05 HUF).
}
\source{
Dobi, B. and Zempléni, A. Markov chain-based cost-optimal control charts for health care data. Quality and Reliability Engineering International, 2019;35(5):1379–1395. https://doi.org/10.1002/qre.2518
}
\references{
Boekholdt SM, Arsenault BJ, Mora S, et al. Association of LDL cholesterol, non–HDL cholesterol, and apolipoprotein B levels with risk of cardiovascular events among patients treated with statins: a meta-analysis. J Am Med Assoc. 2012;307(12):1302-1309. https://doi.org/10. 1001/jama.2012.366

Garmendia F, Brown AS, Reiber I, Adams PC. Attaining United States and European guideline LDL-cholesterol levels with simvastatin in patients with coronary heart disease (the GOALLS study). Curr Med Res Opin. 2000;16(3):208-219. PMID: 11191012.

Lardizabal JA, Deedwania PC. Benefits of statin therapy and compliance in high risk cardiovascular patients. Vasc Health Risk Manag. 2010;6:843-853. https://doi/org/10.2147/VHRM.S9474

High cholesterol: Overview. Institute for Quality and Efficiency in Health Care. https://www.ncbi.nlm.nih.gov/books/NBK279318/ [10 September 2021] Bookshelf ID: NBK279318.
}
\examples{
#Print data
data("LDL")
LDL

###

# Run analysis from Dobi & Zempleni 2019.
# Note that the code below is generated with updated HUF-EUR rate and
# a more accurate R implementation than in the original paper.

stat_LDL_cost <- Markovstat(
  shiftfun = 'exp', h = 50, k = 3.15-LDL$target_value,
  sigma = LDL$standard_deviation,
  s = LDL$num_exp_daily_shifts,
  delta = LDL$expected_shift_size,
  RanRep = TRUE, alpha = LDL$rep_size_first, beta = LDL$rep_size_second,
  RanSam = TRUE, q = LDL$samp_prob_first, z = LDL$samp_prob_second,
  Vd = 100, V = 6-LDL$target_value)
\donttest{
#Defining parallel_opt parallel settings.
#parallel_opt can also be left empty to be defined automatically by the function.
require(parallel)
num_workers <- min(c(detectCores(),2))
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_cost <- Markovchart(
  statdist = stat_LDL_cost,
  OPTIM = TRUE, p = 1,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams = c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)

num_workers <-  min(c(detectCores(),2))
parall      <-  list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_cost_grid <- Markovchart(
  statdist = stat_LDL_cost,
  h=seq(50,75,2.5),
  k=seq(0.05,0.25,0.02),
  p = 1,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams = c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)

require(ggplot2)
plot(res_LDL_cost_grid,
     y = 'Expected \ndaily cost',
     mid = '#ff9494',
     high = '#800000',
     xlab = 'Days between samplings',
     ylab = 'Critical LDL increase') +
     geom_point(aes(x = res_LDL_cost$Parameters[[1]],
                    y = res_LDL_cost$Parameters[[2]]))
}
}
\keyword{datasets}
