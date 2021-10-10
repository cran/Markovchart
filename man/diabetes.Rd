\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
Pseudonymised and randomised time series dataset of diabetes patients for control chart applications
}
\description{
Pseudonymised and randomised time series data of 800 patients. The patients are divided into two main groups by therapy type: insulin analogues (artificial insulins) and glucagon-like peptides (GLP, promotes insulin secretion). The patients' well-being is indicated by the blood sugar (more accurately, the glicated haemoglobin - HbA1c) level.
}
\usage{data("diabetes")}
\format{
  A data frame with 87598 observations on the following 11 variables.
  \describe{
    \item{\code{ID}}{Patient ID}
    \item{\code{DATE}}{Date of the sampling/observation}
    \item{\code{AGE}}{Age of the patient}
    \item{\code{THERAPY}}{Therapy type}
    \item{\code{THERAPY_COST_EUR}}{Therapy cost}
    \item{\code{HEALTHCARE_BURDEN_EUR}}{Event (e.g. heart attack) cost for the health care provider}
    \item{\code{HBA1C_AVG}}{Blood sugar average for the 30-day sampling cycle}
    \item{\code{HBA1C_SD}}{Blood sugar standard deviation for the 30-day sampling cycle}
    \item{\code{SAMPLING_IN_MONTH}}{Number of sampling for the 30-day sampling cycle}
    \item{\code{ICD}}{Diabetes diagnosis type (International Classification of Diseases)}
    \item{\code{THERAPY_VECTOR}}{Therapy vector of the patient, i.e. taking into account the time the therapy lasts after initiation}
  }
}
\details{
The example data focuses on two therapy types: insulin analogues (artificial insulins) and glucagon-like peptides (GLP, promotes insulin secretion). Of course there are more treatment types, the database also lists oral antidiabetics (OAD) and human insulins, but we choose to make the data simpler by focusing on GLP and analogue therapies. For the sake of comparison the therapies are grouped in this way: the first group is insulin analogues with possible parallel OAD therapies but human insulin and GLP excluded. The second group is GLP therapies with possible parallel OAD and insulin analogue therapies but human insulin excluded. Essentially we are comparing the effect and cost of insulin analogues with the effect and cost of additional GLP therapies. For cost calculations, the 2021 March 21 EUR-HUF exchange rate was used (1 EUR = 369.05 HUF).

The example below contains a lengthy code which exemplifies the process of gathering useful data for control chart use. Detailed application of this data can be found in the package's vignette.
}
\source{
The original dataset is based on a month-aggregated time series data of diabetic patients from Hungary which was gathered from the period of 2007 September - 2017 September. The data came from two sources: the National Health Insurance Fund of Hungary and the South-Pest Central Hospital. The first source provided information about diagnoses, treatments, health care event and related costs while the latter provided laboratory data regarding blood sugar level. Patients with International Classification of Diseases (ICD) codes (diagnosis) of E10, E11 and E14, and at least one blood sugar measurement were included initially. Only the data of patients with at least one E11 (type II diabetes) diagnosis in the study period was kept. An additional homogenising filter was the requirement of age above 40 at the time of the first diagnosis. Disease progression and therapy effectiveness estimation required at least two blood sugar (HbA1c) measurements with simultaneous therapy data. A total of 4434 patients satisfied all conditions out of which 2151 had at least two HbA1c measurements.
}
\references{
https://ecmiindmath.org/2019/08/19/markov-chain-based-cost-optimal-control-charts/
}
\examples{
data("diabetes")
str(diabetes)

\dontrun{
##### Example of data assessment for control chart use #####
### Packages
require(zoo)
require(reshape2)

RANDOMISED_DATA <- diabetes

### Functions
weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
na.rm)
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

### Setting up data
# Way too high HbA1C levels are discarded as outliers
RANDOMISED_DATA$HBA1C_AVG[RANDOMISED_DATA$HBA1C_AVG>20 & !is.na(RANDOMISED_DATA$HBA1C_AVG)] <- NA

# Lowest HbA1c level taken into account
lowest <- 4

### Gathering data for several estimates
RANDOMISED_DATA <- RANDOMISED_DATA[RANDOMISED_DATA$ID %in%
RANDOMISED_DATA$ID[grepl("E11",RANDOMISED_DATA$ICD)],]

# Process standard deviation
sigma_param <- sigma <- sqrt(weighted.mean((RANDOMISED_DATA$HBA1C_SD[
RANDOMISED_DATA$SAMPLING_IN_MONTH>=2 & !is.na(RANDOMISED_DATA$SAMPLING_IN_MONTH)])^2,
RANDOMISED_DATA$SAMPLING_IN_MONTH[RANDOMISED_DATA$SAMPLING_IN_MONTH>=2 &
!is.na(RANDOMISED_DATA$SAMPLING_IN_MONTH)]))

IDLIST <- unique(RANDOMISED_DATA$ID[!is.na(RANDOMISED_DATA$HBA1C_AVG)][
duplicated(RANDOMISED_DATA$ID[!is.na(RANDOMISED_DATA$HBA1C_AVG)])])
IDLIST <- unique(RANDOMISED_DATA$ID[(RANDOMISED_DATA$ID %in% IDLIST) & RANDOMISED_DATA$AGE>39])

shiftdat <- NULL
stimedat <- NULL
repaidat <- NULL
deltats  <- NULL
deltaATC <- NULL
for(i in IDLIST)
{
	smalldat    <- RANDOMISED_DATA[RANDOMISED_DATA$ID==i,c("DATE","HBA1C_AVG","THERAPY_VECTOR")]
	smalldat    <- smalldat[!is.na(smalldat$DATE) & !is.na(smalldat$HBA1C_AVG),]
	patshiftdat <- as.data.frame(cbind(i,smalldat$DATE[2:dim(smalldat)[1]],diff(smalldat$DATE),
	diff(smalldat$HBA1C_AVG))[diff(smalldat$HBA1C_AVG)>2*sigma,,drop=FALSE])
	if(dim(patshiftdat)[1]>1) stimedat <- rbind(stimedat,cbind(i,diff(as.Date(patshiftdat$V2))))
	patrepaidat <- as.data.frame(cbind(i,diff(smalldat$DATE),(smalldat$HBA1C_AVG-lowest)[2:
    length(smalldat$HBA1C_AVG)]/(smalldat$HBA1C_AVG-lowest)[1:(length(smalldat$HBA1C_AVG)-1)],
		as.character(smalldat$THERAPY_VECTOR[1:(length(smalldat$THERAPY_VECTOR)-1)]))[
		(which(diff(smalldat$HBA1C_AVG)<(-2*sigma) &
		smalldat$HBA1C_AVG[1:(length(smalldat$HBA1C_AVG)-1)]>6 &
		smalldat$HBA1C_AVG[1:(length(smalldat$HBA1C_AVG)-1)]<=20)),,drop=FALSE])

	shiftdat <- rbind(shiftdat,patshiftdat)
	repaidat <- rbind(repaidat,patrepaidat)
	deltats  <- rbind(deltats,cbind(i,diff(as.Date(RANDOMISED_DATA$DATE[
	!is.na(RANDOMISED_DATA$HBA1C_AVG) & RANDOMISED_DATA$ID==i]))))
	try(deltaATC <- rbind(deltaATC,cbind(i,diff(as.Date(RANDOMISED_DATA$DATE[
	!is.na(RANDOMISED_DATA$THERAPY) & RANDOMISED_DATA$ID==i])))), silent=TRUE)
}
colnames(shiftdat) <- c("ID","TIME","TIMEDIFF","SHIFTSIZE")
colnames(deltats)  <- c("ID","DeltaT")
colnames(deltaATC) <- c("ID","deltaATC")

# delta: expected shift size
delta_param <- mean(shiftdat$SHIFTSIZE[shiftdat$TIMEDIFF>=90 & shiftdat$TIMEDIFF<184])

# Empirical distribution of elapsed times (between samplings)
summary(deltats[,2])
mean(deltats[,2])
median(deltats[,2])
sd(deltats[,2])

# s: expected number of shifts per unit time
stimedat           <- as.data.frame(stimedat)
colnames(stimedat) <- c("ID","TIMEDIFF")
s_param            <- 1/mean(stimedat$TIMEDIFF[stimedat$TIMEDIFF<367])

# alphas, betas: therapy effectiveness parameters
colnames(repaidat) <- c("ID","TIMEDIFF","REMAIN","THERAP")
repaidat$REMAIN    <- as.numeric(as.character(repaidat$REMAIN))
repaidat$TIMEDIFF  <- as.numeric(as.character(repaidat$TIMEDIFF))
repaidat$THERAP    <- as.character(repaidat$THERAP)
repaidat           <- repaidat[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184,]
repaidat$REMAIN[repaidat$REMAIN<0] <- 0

ther_eff <- as.data.frame(rbind(
cbind("ANALOGUE",repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
grepl("ANALOGUE",repaidat$THERAP) & !grepl("-H",repaidat$THERAP) & !grepl("GLP",repaidat$THERAP)]),
cbind("GLP",repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
grepl("GLP",repaidat$THERAP) & !grepl("-H",repaidat$THERAP)])))
ther_eff[,1]       <- factor(ther_eff[,1], levels = c("ANALOGUE", "GLP"))
ther_eff[,2]       <- as.numeric(as.character(ther_eff[,2]))
colnames(ther_eff) <- c("Type","Effectiveness")

ANALOGUE <- estBetaParams(mean(repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
grepl("ANALOGUE",repaidat$THERAP) & !grepl("-H",repaidat$THERAP) & !grepl("GLP",repaidat$THERAP)]),
			 var(repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
			 grepl("ANALOGUE",repaidat$THERAP) & !grepl("-H",repaidat$THERAP) &
			 !grepl("GLP",repaidat$THERAP)]))
GLP <- estBetaParams(mean(repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
                       grepl("GLP",repaidat$THERAP) & !grepl("-H",repaidat$THERAP)]),
                     var(repaidat$REMAIN[repaidat$TIMEDIFF>=90 & repaidat$TIMEDIFF<184 &
                       grepl("GLP",repaidat$THERAP) & !grepl("-H",repaidat$THERAP)]))

### Repair cost
HBA1C_fill <- NULL
for (i in unique(RANDOMISED_DATA$ID[!is.na(RANDOMISED_DATA$HBA1C_AVG)]))
{
  vec <- RANDOMISED_DATA$HBA1C_AVG[RANDOMISED_DATA$ID==i]
  vec[which(is.na(vec))[which(is.na(vec))<which(!is.na(vec))[1]]] <- vec[which(!is.na(vec))[1]]
  vec[which(is.na(vec))[which(is.na(vec))>which(!is.na(vec))[length(which(!is.na(vec)))]]] <-
    vec[which(!is.na(vec))[length(which(!is.na(vec)))]]
  vec <- na.approx(vec)
  HBA1C_fill <- rbind(HBA1C_fill,cbind(i,vec))

  smaldat <- RANDOMISED_DATA[RANDOMISED_DATA$ID==i,]
  smaldat$THERAPY_COST_EUR[smaldat$THERAPY_COST_EUR==0 & smaldat$THERAPY_VECTOR!=""] <- NA
  if(is.na(smaldat$THERAPY_COST_EUR[1])) smaldat$THERAPY_COST_EUR[1]                 <- 0
  new_burden <- na.locf(smaldat$THERAPY_COST_EUR)

  seged                     <- cbind(rle(is.na(smaldat$THERAPY_COST_EUR))[[2]],
                                     rle(is.na(smaldat$THERAPY_COST_EUR))[[1]])
  seged[,2][seged[,1]==0]   <- seged[,2][seged[,1]==0]-1
  seged[,2][seged[,1]==1]   <- seged[,2][seged[,1]==1]+1
  if(seged[length(seged[,1]),1]==0) seged[length(seged[,2]),2] <- seged[length(seged[,2]),2]+1
  seged2                    <- cbind(rep(seged[,1], seged[,2]),rep(seged[,2], seged[,2]))
  new_burden[seged2[,1]==1] <- new_burden[seged2[,1]==1]/seged2[,2][seged2[,1]==1]

  RANDOMISED_DATA$THERAPY_COST_EUR[RANDOMISED_DATA$ID==i]	<-	new_burden
}

RANDOMISED_DATA$HBA1C_fill <- NA
RANDOMISED_DATA$HBA1C_fill[RANDOMISED_DATA$ID%in%HBA1C_fill[,1]] <- HBA1C_fill[,2]
RANDOMISED_DATA$HBA1C_fill_filter <- RANDOMISED_DATA$HBA1C_fill
RANDOMISED_DATA$HBA1C_fill_filter[RANDOMISED_DATA$HBA1C_fill_filter>=10] <- NA

discparam    <-	150
cutHBA1C_AVG <-	cut(na.omit(RANDOMISED_DATA$HBA1C_fill_filter),breaks=discparam)
newlvls      <-	seq(min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)),
                    max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)),
                    (max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter))-
                      min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)))/discparam)[1:discparam] +
                    (max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter))-
                      min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)))/discparam/2
levels(cutHBA1C_AVG) <- newlvls
costs                <- cbind(as.numeric(as.character(cutHBA1C_AVG)),
                               RANDOMISED_DATA$THERAPY_COST_EUR[!is.na(
                               RANDOMISED_DATA$HBA1C_fill_filter)]/30,
                               as.character(RANDOMISED_DATA$THERAPY[
                               !is.na(RANDOMISED_DATA$HBA1C_fill_filter)]))
costs           <- as.data.frame(costs)
colnames(costs)	<- c("HBA1C","HC_BURDEN","THERAP")
costs$HBA1C     <- as.numeric(as.character(costs$HBA1C))
costs$HC_BURDEN	<- as.numeric(as.character(costs$HC_BURDEN))
costs$THERAP    <- as.character(costs$THERAP)

costs$THERAP[grepl("ANALOGUE", costs$THERAP) & !grepl("GLP", costs$THERAP)] <- "ANALOGUE"
costs$THERAP[grepl("GLP",costs$THERAP)]                                     <- "GLP"

cost.ANALOGUE <- as.data.frame(cbind(sort(unique(costs$HBA1C[costs$THERAP=="ANALOGUE"])),
                                as.numeric(tapply(costs$HC_BURDEN[costs$THERAP=="ANALOGUE"],
                                  costs$HBA1C[costs$THERAP=="ANALOGUE"],mean))))
colnames(cost.ANALOGUE) <- c("HBA1C","HC_BURDEN")

cost.GLP <-	as.data.frame(cbind(sort(unique(costs$HBA1C[costs$THERAP=="GLP"])),
                            as.numeric(tapply(costs$HC_BURDEN[costs$THERAP=="GLP"],
                              costs$HBA1C[costs$THERAP=="GLP"],mean))))
colnames(cost.GLP) <-	c("HBA1C","HC_BURDEN")

## ANALOGUE therapy
# Mean
cost.ANALOGUE           <- na.omit(as.data.frame(cbind(as.numeric(
                                     costs$HBA1C[costs$THERAP=="ANALOGUE"]),
                                     costs$HC_BURDEN[costs$THERAP=="ANALOGUE"])))
colnames(cost.ANALOGUE) <- c("HBA1C","HC_BURDEN")
cost.ANALOGUE           <- cost.ANALOGUE[order(cost.ANALOGUE$HBA1C),]
cost.ANALOGUE           <- cost.ANALOGUE[cost.ANALOGUE$HBA1C>lowest,]
cost.ANALOGUE$HBA1C     <- cost.ANALOGUE$HBA1C-min(lowest)

# Fit non-linear model
mod.ANALOGUE <- nls(HC_BURDEN ~  a + b/(HBA1C + c),
                    start = list(a = 5, b = -5, c = 1), cost.ANALOGUE,
                    control = list(maxiter = 50000, minFactor = 0.000000000000001))

cost_ANALOGUE_plotdat <- cbind("ANALOGUE",as.data.frame(cbind(seq(0,6,6/99),
                                predict(mod.ANALOGUE,
                                 newdata=data.frame(HBA1C = seq(0,6,6/99))))))

# Variance
cost_var.ANALOGUE  <- na.omit(as.data.frame(cbind(sort(unique(
                               costs$HBA1C[costs$THERAP=="ANALOGUE"])),
                               as.numeric(tapply(costs$HC_BURDEN[costs$THERAP=="ANALOGUE"],
                               costs$HBA1C[costs$THERAP=="ANALOGUE"],var)))))
colnames(cost_var.ANALOGUE)	<- c("HBA1C","HC_BURDEN")
cost_var.ANALOGUE           <- cost_var.ANALOGUE[cost_var.ANALOGUE$HBA1C>lowest,]
cost_var.ANALOGUE$HBA1C     <- cost_var.ANALOGUE$HBA1C-min(lowest)

# Fit non-linear model
mod_var.ANALOGUE <- nls(HC_BURDEN ~  a + b/(HBA1C + c),
                        start = list(a = 5, b = -3, c = 0.1),
                        cost_var.ANALOGUE[cost_var.ANALOGUE$HBA1C<10-lowest,],
                        control = list(maxiter = 50000, minFactor = 0.000000000000001))

cost_ANALOGUE_plotdat <- cbind(cost_ANALOGUE_plotdat,
                               cost_ANALOGUE_plotdat[,3] -
                                 sqrt(predict(mod_var.ANALOGUE,
                                   newdata=data.frame(HBA1C = seq(0,6,6/99)))),
                               cost_ANALOGUE_plotdat[,3] +
                                 sqrt(predict(mod_var.ANALOGUE,
                                   newdata=data.frame(HBA1C = seq(0,6,6/99)))))
colnames(cost_ANALOGUE_plotdat) <- c("Data","HbA1c","Therapy cost","low","high")

## GLP
# Mean
cost.GLP           <- na.omit(as.data.frame(cbind(as.numeric(
                                costs$HBA1C[costs$THERAP=="GLP"]),
                                costs$HC_BURDEN[costs$THERAP=="GLP"])))
colnames(cost.GLP) <- c("HBA1C","HC_BURDEN")
cost.GLP           <- cost.GLP[order(cost.GLP$HBA1C),]
cost.GLP           <- cost.GLP[cost.GLP$HBA1C>lowest,]
cost.GLP$HBA1C     <- cost.GLP$HBA1C-min(lowest)

# Simple linear
mod.GLP <- nls(HC_BURDEN ~ a + b * HBA1C,
               start = list(a = 1, b = 1), cost.GLP,
               control = list(maxiter = 50000, minFactor = 0.000000000000001))

cost_GLP_plotdat <-	cbind("GLP",as.data.frame(cbind(seq(0,6,6/99),
                     predict(mod.GLP, newdata=data.frame(HBA1C = seq(0,6,6/99))))))

# Variance
cost_var.GLP           <- na.omit(as.data.frame(cbind(sort(unique(
                                   costs$HBA1C[costs$THERAP=="GLP"])),
                                   as.numeric(tapply(costs$HC_BURDEN[costs$THERAP=="GLP"],
                                   costs$HBA1C[costs$THERAP=="GLP"],var)))))
colnames(cost_var.GLP) <- c("HBA1C","HC_BURDEN")
cost_var.GLP           <- cost_var.GLP[cost_var.GLP$HBA1C>lowest,]
cost_var.GLP$HBA1C     <- cost_var.GLP$HBA1C-min(lowest)

# Simple linear
mod_var.GLP <- nls(HC_BURDEN ~  a + b*(HBA1C),
                   start = list(a = 5, b = -3), cost_var.GLP,
                   control = list(maxiter = 50000, minFactor = 0.000000000000001))

cost_GLP_plotdat <- cbind(cost_GLP_plotdat,
                          cost_GLP_plotdat[,3] -
                           sqrt(predict(mod_var.GLP, newdata=data.frame(HBA1C = seq(0,6,6/99)))),
                          cost_GLP_plotdat[,3] +
                           sqrt(predict(mod_var.GLP, newdata=data.frame(HBA1C = seq(0,6,6/99)))))
colnames(cost_GLP_plotdat) <- c("Data","HbA1c","Therapy cost","low","high")

### Out-of-control cost
COST_CUMU<-NULL
for (i in unique(RANDOMISED_DATA$ID[!is.na(RANDOMISED_DATA$HEALTHCARE_BURDEN_EUR)]))
{
	vec       <- RANDOMISED_DATA$HEALTHCARE_BURDEN_EUR[RANDOMISED_DATA$ID==i]
	vec2      <- rollapply(vec,min(6,length(vec)),
	              sum,align="left",partial=TRUE)/
	               (pmin(length(vec)-(1:length(vec))+1,6)*30)
	COST_CUMU <- rbind(COST_CUMU,cbind(i,vec2))
}

RANDOMISED_DATA$COST_CUMU <- NA
RANDOMISED_DATA$COST_CUMU[RANDOMISED_DATA$ID%in%COST_CUMU[,1]] <- COST_CUMU[,2]

discparam    <- 150
cutHBA1C_AVG <- cut(na.omit(RANDOMISED_DATA$HBA1C_fill_filter),breaks=discparam)
newlvls      <- seq(min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)),
                    max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)),
                    (max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter))-
                      min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)))/discparam)[1:discparam] +
                      (max(na.omit(RANDOMISED_DATA$HBA1C_fill_filter))-
                        min(na.omit(RANDOMISED_DATA$HBA1C_fill_filter)))/discparam/2
levels(cutHBA1C_AVG) <- newlvls
ooc_costs <- cbind(round(as.numeric(as.character(cutHBA1C_AVG)),1),
              RANDOMISED_DATA$COST_CUMU[!is.na(RANDOMISED_DATA$HBA1C_fill_filter)])
ooc_costs <- as.data.frame(ooc_costs)

# Mean
disc_ooc_cost           <- as.data.frame(cbind(as.numeric(ooc_costs[,1]),ooc_costs[,2]))
colnames(disc_ooc_cost)	<- c("HBA1C","COST")
disc_ooc_cost           <- disc_ooc_cost[order(disc_ooc_cost$HBA1C),]
disc_ooc_cost           <- disc_ooc_cost[disc_ooc_cost$HBA1C >= lowest,]
disc_ooc_cost$HBA1C     <- disc_ooc_cost$HBA1C - lowest

mod.COST <- nls(COST ~ a + c*HBA1C^2, start = list(a = 1, c = 1), disc_ooc_cost)

cost_COST_plotdat <- cbind("Complications",as.data.frame(cbind(seq(0, 6, 6/99),
                           predict(mod.COST, newdata=data.frame(HBA1C = seq(0, 6, 6/99))))))

# Variance
disc_ooc_cost_var           <- as.data.frame(cbind(sort(unique(ooc_costs[,1])),
                                as.numeric(tapply(ooc_costs[,2],ooc_costs[,1],var))))
colnames(disc_ooc_cost_var) <- c("HBA1C","COST")

disc_ooc_cost_var       <- disc_ooc_cost_var[disc_ooc_cost_var$HBA1C>=lowest,]
disc_ooc_cost_var$HBA1C <- disc_ooc_cost_var$HBA1C-lowest

mod_var.COST <- nls(COST ~ a + c*HBA1C^2,
                    start = list(a = 0.5, c = 0.5), disc_ooc_cost_var,
                    control = list(maxiter = 50000, minFactor = 0.000000000000001))

cost_COST_plotdat <- cbind(cost_COST_plotdat,
                           cost_COST_plotdat[,3] -
                             sqrt(predict(mod_var.COST,
                               newdata=data.frame(HBA1C = seq(0,6,6/99)))),
                           cost_COST_plotdat[,3] +
                             sqrt(predict(mod_var.COST,
                               newdata=data.frame(HBA1C = seq(0,6,6/99)))))
colnames(cost_COST_plotdat) <- c("Data","HbA1c","Therapy cost","low","high")

cost_plots       <- rbind(cost_ANALOGUE_plotdat,cost_GLP_plotdat,cost_COST_plotdat)
cost_plots$HbA1c <- cost_plots$HbA1c + lowest
cost_plots[,"Therapy cost"]               <- cost_plots[,"Therapy cost"]/1
cost_plots[,"low"]                        <- cost_plots[,"low"]/1
cost_plots[,"low"][cost_plots[,"low"]<0]  <- 0
cost_plots[,"high"]                       <- cost_plots[,"high"]/1

cost_plots <- cost_plots

### Sampling cost: official, government-regulated prices related to sampling
### converted from HUF to EUR
sampling_cost=2875/369.05

### Empirical costs for comparison
# GLP
mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30 +
mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU) +
sampling_cost/mean(deltats[,2])

sd(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR/30 +
RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU +
sampling_cost/mean(deltats[,2]))

# ANALOGUE
mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30 +
mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU) +
sampling_cost/mean(deltats[,2])

sd(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR/30 +
RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU +
sampling_cost/mean(deltats[,2]))

### Empirical HbA1c distribution
# GLP
empi.GLP  <- RANDOMISED_DATA$HBA1C_fill[grepl("GLP", RANDOMISED_DATA$THERAPY) &
               RANDOMISED_DATA$HBA1C_fill>=4 & RANDOMISED_DATA$HBA1C_fill<=20]
cutHBA1C  <- cut(na.omit(empi.GLP),breaks=100)
newlvls   <- seq(min(na.omit(empi.GLP)),max(na.omit(empi.GLP)),
                     (max(na.omit(empi.GLP))-min(na.omit(empi.GLP)))/100)[1:100] +
                      (max(na.omit(empi.GLP))-min(na.omit(empi.GLP)))/100/2
levels(cutHBA1C)    <- newlvls
empi.GLP            <- as.data.frame(table(cutHBA1C)/sum(table(cutHBA1C)))
empi.GLP$cutHBA1C   <- as.numeric(as.character(empi.GLP$cutHBA1C))
empi.GLP            <- cbind("GLP", empi.GLP)
colnames(empi.GLP)  <- c("Therapy", "HbA1c", "Probability")

# ANALOGUE
empi.ANALOGUE   <- RANDOMISED_DATA$HBA1C_fill[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
                    !grepl("GLP", RANDOMISED_DATA$THERAPY) &
                      RANDOMISED_DATA$HBA1C_fill>=4 & RANDOMISED_DATA$HBA1C_fill<=20]
cutHBA1C        <- cut(na.omit(empi.ANALOGUE),breaks=100)
newlvls         <- seq(min(na.omit(empi.ANALOGUE)),
                       max(na.omit(empi.ANALOGUE)),
                       (max(na.omit(empi.ANALOGUE))-
                         min(na.omit(empi.ANALOGUE)))/100)[1:100] +
                          (max(na.omit(empi.ANALOGUE))-
                            min(na.omit(empi.ANALOGUE)))/100/2
levels(cutHBA1C)        <- newlvls
empi.ANALOGUE           <- as.data.frame(table(cutHBA1C)/sum(table(cutHBA1C)))
empi.ANALOGUE$cutHBA1C  <- as.numeric(as.character(empi.ANALOGUE$cutHBA1C))
empi.ANALOGUE           <- cbind("ANALOGUE", empi.ANALOGUE)
colnames(empi.ANALOGUE)	<- c("Therapy", "HbA1c", "Probability")

# Merging datasets
hba1c_empir <- rbind(empi.ANALOGUE, empi.GLP)

# The list of gathered data and statistics:
# sigma_param, s_param, delta_param, ANALOGUE
# GLP, mod.ANALOGUE, mod_var.ANALOGUE
# mod.GLP, mod_var.GLP, mod.COST, mod_var.COST
# cost_plots, sampling_cost, hba1c_empir

}
}
\keyword{datasets}
