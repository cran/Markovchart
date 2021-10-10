## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(Markovchart)
require(ggplot2)
require(grid)
require(gridExtra)
require(kableExtra)
require(reshape2)
require(zoo)
require(parallel)

set.seed(2020)


## ----lipidtable, echo = FALSE-------------------------------------------------

ldltable <- as.data.frame(cbind(
  c("3 mmol/l","0.1 mmol/l","0.8/3","1/120","0.027, 1.15","0.1, 30","5.01 EUR","4.60 EUR","9.97 EUR","7.48 EUR"),
  c("Target value"," Process standard deviation","Expected shift size, 0.8 increase in LDL per year on average","Expected number of shifts in a day, 3 shifts per year on average","Parameters of the repair size beta distribution","Parameters of the sampling probability logistic function","Sampling cost"," Shift-proportional daily out-of-control (OOC) cost","Base repair cost","Shift-proportional repair cost"),
  c("Set according to the European guideline for patients at risk","Estimated using real life data from Hungary, namely registry data through Healthware Consulting Ltd.","Estimated with the help of a health professional","Estimated with the help of a health professional","Estimated using an international study which included Hungary","Patient non-compliance in LDL controlling medicine is quite high, and this is represented through the parametrisation of the logistic function","Estimated using the LDL testing cost and visit cost in Hungary","Estimated using real world data of cardiovascular event costs from Hungary","Estimated using the simvastatin therapy costs in Hungary","Estimated using the simvastatin therapy costs in Hungary")))
colnames(ldltable) <- c("Parameter value","Meaning","Parameter source")

kable_styling(kbl(ldltable))


## ----app_LDL_cost,  echo=TRUE, fig.height=3, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
data("LDL")

stat_LDL_cost <- Markovstat(
  shiftfun = 'exp', h = 50, k = 0.15,
  sigma = LDL$standard_deviation, s = LDL$num_exp_daily_shifts,
  delta = LDL$expected_shift_size,
  RanRep = TRUE, alpha = LDL$rep_size_first, beta = LDL$rep_size_second,
  RanSam = TRUE, q = LDL$samp_prob_first, z = LDL$samp_prob_second,
  Vd = 50, V = 3)

#Defining parallel_opt parallel settings.
#parallel_opt can also be left empty to be defined automatically by the function.
num_workers <-  min(c(detectCores(),2))
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_cost <- Markovchart(
  statdist = stat_LDL_cost,
  OPTIM=TRUE, p = 1,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams = c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)
res_LDL_cost


## ----app_LDL_cost_plot,  echo=TRUE, fig.height=3, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_cost_grid <- Markovchart(
  statdist = stat_LDL_cost,
  h=seq(50,75,2.5),
  k=seq(0.05,0.25,0.02),
  p = 1,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams= c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)
plot(res_LDL_cost_grid,
     y = 'Expected \ndaily cost',
     mid = '#ff9494',
     high = '#800000',
     xlab = 'Days between samplings',
     ylab = 'Critical LDL increase') +
     geom_point(aes(x = res_LDL_cost$Parameters[[1]],
                    y = res_LDL_cost$Parameters[[2]]))

## ----app_LDL_G,  echo=TRUE, fig.height=3, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
stat_LDL_cost <- Markovstat(
  shiftfun = 'exp', h = 50, k = 0.15,
  sigma = LDL$standard_deviation, s = LDL$num_exp_daily_shifts,
  delta = LDL$expected_shift_size,
  RanRep = TRUE, alpha = LDL$rep_size_first, beta = LDL$rep_size_second,
  RanSam=TRUE, q = LDL$samp_prob_first, z = LDL$samp_prob_second,
  Vd = 50, V = 3)
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_G <- Markovchart(
  statdist = stat_LDL_cost,
  OPTIM=TRUE, p = 0.9,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams= c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)
res_LDL_G

## ----app_LDL_G_plot,  echo=TRUE, fig.height=3, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_LDL_G_grid <- Markovchart(
  statdist = stat_LDL_cost,
  h=seq(50,75,2.5),
  k=seq(0.05,0.25,0.02),
  p = 0.9,
  cs = LDL$sampling_cost,
  cf = LDL$base_rep_cost,
  coparams = c(0,LDL$OOC_cost),
  crparams= c(LDL$base_rep_cost,LDL$prop_rep_cost),
  parallel_opt = parall)
plot(res_LDL_G_grid,
		 y = 'Expected \ndaily cost',
		 mid = '#ff9494',
     high = '#800000',
     xlab = 'Days between samplings',
     ylab = 'Critical LDL increase',
		 nbreaks = 14) +
	   geom_point(aes(x = res_LDL_G$Parameters[[1]],
					          y = res_LDL_G$Parameters[[2]]))

## ----app_LDL_params,  echo=TRUE, fig.height=5, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
results <- NULL
statds  <- NULL
for	(i in 3:10)
{
	for	(j in c(0.9,1))
	{
	  parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
    res_LDL_optim <- Markovchart(
     statdist = stat_LDL_cost,
     OPTIM = TRUE,
     p = j,
     cs = LDL$sampling_cost,
     cf = LDL$base_rep_cost,
     coparams = c(0,i),
     crparams= c(LDL$base_rep_cost,LDL$prop_rep_cost),
     parallel_opt = parall)

   		results <- rbind(results,c(j,i,unlist(c(res_LDL_optim$Results[c(2,3)],res_LDL_optim$Parameters))))
   		statds  <- rbind(statds,res_LDL_optim$Stationary_distribution)
	}
}

results	<-	resultsbck		<-	as.data.frame(results)
colnames(results)	<-	c("p","co","EC","SDC","h","k")

results$p			<-	as.character(results$p)
results				<-	cbind(results[,1:2],melt(results[,3:6]))
results$variable	<-	as.character(results$variable)
results$variable[results$variable=="h"]		<-	"Days between samplings"
results$variable[results$variable=="k"]		<-	"Critical LDL increase"
results$variable[results$variable=="EC"]	<-	"Expected cost"
results$variable[results$variable=="SDC"]	<-	"Cost standard deviation"
results$variable = factor(results$variable, levels=c("Critical LDL increase","Days between samplings","Expected cost","Cost standard deviation"))
results$min	<-	NA

results$min[results$variable=="Days between samplings"]	<- min(results$value[results$variable=="Days between samplings"])
results$max[results$variable=="Days between samplings"]	<- max(results$value[results$variable=="Days between samplings"])
results$min[results$variable=="Critical LDL increase"]	<- min(results$value[results$variable=="Critical LDL increase"])
results$max[results$variable=="Critical LDL increase"]	<- max(results$value[results$variable=="Critical LDL increase"])
results$min[results$variable=="Expected cost"]	<- min(results$value[results$variable=="Expected cost"])
results$max[results$variable=="Expected cost"]	<- max(results$value[results$variable=="Expected cost"])
results$min[results$variable=="Cost standard deviation"]	<- min(results$value[results$variable=="Cost standard deviation"])
results$max[results$variable=="Cost standard deviation"]	<- max(results$value[results$variable=="Cost standard deviation"])

app_LDL_params_plot <- ggplot(results, aes(x=co, y=value, group=p)) +
		geom_line(aes(colour=p),size = 1.1) +
		facet_wrap(~variable, labeller = label_bquote(rows=.(variable)), scales="free_y", nrow=2) +
		scale_colour_manual(expression(italic(p)),values=c("black","#800000")) +
		geom_blank(aes(y = min)) +
		geom_blank(aes(y = max)) +
		xlab("Out-of-control cost") +
		ylab("Value") +
		theme_bw() +
		theme(strip.text.x = element_text(size = 11),
		strip.text.y = element_text(size = 11,angle = 0), legend.position="top")

app_LDL_params_plot

## ----diab_aggr,  include=FALSE, fig.height=3, fig.width=6, message=FALSE, warning=FALSE, eval=TRUE----
data("diabetes")

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
	patshiftdat	<- as.data.frame(cbind(i,as.data.frame(smalldat$DATE[2:dim(smalldat)[1]]),diff(smalldat$DATE),
	diff(smalldat$HBA1C_AVG))[diff(smalldat$HBA1C_AVG)>2*sigma,,drop=FALSE])
	if(dim(patshiftdat)[1]>1) stimedat <- rbind(stimedat,cbind(i,diff(as.Date(patshiftdat[,2]))))
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

## ----diabtable, echo = FALSE--------------------------------------------------

diabetestable <- as.data.frame(cbind(
                           c("Total (all have E11 and are over 40)","E11 ICD"),
                           c(800,492),
                           c(630,272),
                           c(170,99)))

colnames(diabetestable) <- c("Patient group","Total","Analogue", "GLP")

kable_styling(kbl(diabetestable))

## ----therap_plot_fig, echo = TRUE, fig.height = 3, fig.width = 6--------------
therap_plot	<-	ggplot(data.frame(x = c(0,1)), aes(x)) + 
						  stat_function(fun = dbeta, aes(colour = ' Analogue', linetype =' Analogue'), size = 1, args = list(shape1 = unlist(ANALOGUE)[1], shape2 = unlist(ANALOGUE)[2])) + 
						  stat_function(fun = dbeta, aes(colour = 'GLP', linetype = 'GLP'), size = 1, args = list(shape1 = unlist(GLP)[1], shape2 = unlist(GLP)[2]), alpha = 0.9) + 
						  scale_colour_manual('Therapy type', values = c('black', '#800000')) + 
            	scale_linetype_manual('Therapy type', values = c('solid', 'dashed')) + 
						  xlab('Therapy effectiveness \n(HbA1c level after therapy divided by the level before)') + 
						  ylab('Beta distribution density') + 
						  theme_bw() + 
						  theme(plot.title = element_text(hjust = 0.5, size = 11), legend.justification = c(1, 0),
												legend.position = c(0.20, 0.63), legend.title = element_blank(), legend.margin = margin(t = 0, r = 0.1, b = 0, l = 0, unit = 'cm'),
												legend.key = element_rect(fill = 'white', colour = 'white'))
therap_plot

## ----cost_plot_fig, echo = TRUE, fig.height = 3, fig.width = 6----------------
cost_plot <- ggplot(data = cost_plots, aes(x = HbA1c, y = `Therapy cost`)) + 
					geom_ribbon(aes(x = HbA1c, ymin = low, ymax = high, fill = '±SD'), alpha = 0.4) + 
					geom_line(aes(x = HbA1c, y = `Therapy cost`, col = 'Fitted line'), size = 1) + 
					facet_wrap(Data ~ . ,labeller = label_bquote(rows = .(Data)), ncol = 3) + 
					scale_color_manual(name = '', values = c('Fitted line' = '#800000')) + 
					scale_fill_manual(name = '', values = c('±SD' = 'grey70')) + 
          ylab('Daily cost in EUR') + 
					theme_bw() + 
					guides(color = guide_legend(order = 1),
                 fill = guide_legend(order = 2)) + 
					theme(legend.position = 'top')

cost_plot

## ----customfuns---------------------------------------------------------------
crfun_ANALOGUE <- function(mudist, crparams){
  mudist <- mudist
  crb    <- crparams[1]
  crs    <- crparams[2]
  crs2   <- crparams[3]
  cr     <- crb + crs / (mudist + crs2)
  return(cr)}
vcrfun_ANALOGUE <- function(mudist, vcrparams){
  mudist <- mudist
  vcrb   <- vcrparams[1]
  vcrs   <- vcrparams[2]
  vcrs2  <- vcrparams[3]
  vcr    <- vcrb + vcrs / (mudist + vcrs2)
  return(vcr)}

stat_ANALOGUE <- Markovstat(
  shiftfun = 'exp', h = 206.22, k = 3,
  sigma = sigma_param, s = s_param,
  delta = delta_param, RanRep = TRUE,
  alpha = as.numeric(ANALOGUE[1]),
  beta = as.numeric(ANALOGUE[2]),
  Vd = 100, V = 18)
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_ANALOGUE <- Markovchart(
  statdist = stat_ANALOGUE, p = 1,
  constantr = TRUE, ooc_rep = 1,
  cs = sampling_cost,
  coparams = summary(mod.COST)$coef[ , 1],
  crfun = crfun_ANALOGUE,
  crparams = summary(mod.ANALOGUE)$coef[ , 1],
  vcoparams = summary(mod_var.COST)$coef[ , 1],
  vcrfun = vcrfun_ANALOGUE,
  vcrparams = summary(mod_var.ANALOGUE)$coef[ , 1],
  parallel_opt = parall)

stat_GLP <- Markovstat(
  shiftfun = 'exp', h = 206.22, k = 3,
  sigma = sigma_param, s = s_param,
  delta = delta_param, RanRep = TRUE,
  alpha = as.numeric(GLP[1]),
  beta = as.numeric(GLP[2]),
  Vd = 100, V = 18)
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
res_GLP <- Markovchart(
  statdist = stat_GLP, p = 1,
  constantr = TRUE, ooc_rep = 1,
  cs = sampling_cost,
  coparams = summary(mod.COST)$coef[ , 1],
  crparams = summary(mod.GLP)$coef[ , 1],
  vcoparams = summary(mod_var.COST)$coef[ , 1],
  vcrparams = summary(mod_var.GLP)$coef[ , 1],
  parallel_opt = parall)

## ----statd_therapies_fig, echo = TRUE, fig.height = 3, fig.width = 6----------
int             <- seq(4, 22, by = (18 / (100 - 1))) - (18 / (100 - 1)) / 2
int[1]          <- 4

distance_dist_A <- res_ANALOGUE[[4]][c(1, (100 + 2):(100 * 2))] + 
                   res_ANALOGUE[[4]][2:(100 + 1)]
theo.ANALOGUE   <- cbind('Analogue', as.data.frame(cbind(int,distance_dist_A)))
colnames(theo.ANALOGUE) <- c('Therapy', 'HbA1c', 'Probability')

distance_dist_G    <- res_GLP[[4]][c(1,(100 + 2):(100 * 2))] + 
                   res_GLP[[4]][2:(100 + 1)]
theo.GLP			     <- cbind('GLP', as.data.frame(cbind(int,distance_dist_G)))
colnames(theo.GLP) <- c('Therapy', 'HbA1c', 'Probability')

hba1c_theo  <- rbind(theo.ANALOGUE, theo.GLP)
hba1c_empir$Therapy[hba1c_empir$Therapy=='ANALOGUE'] <- 'Analogue'

hba1c_empir$Density <- NA
hba1c_empir$Density[hba1c_empir$Therapy=='Analogue'] <-
  hba1c_empir$Probability[hba1c_empir$Therapy=='Analogue']/
  ((max(hba1c_empir$HbA1c[hba1c_empir$Therapy=='Analogue'])-min(hba1c_empir$HbA1c[hba1c_empir$Therapy=='Analogue']))/(100-1))
hba1c_empir$Density[hba1c_empir$Therapy=='GLP'] <-
  hba1c_empir$Probability[hba1c_empir$Therapy=='GLP']/
  ((max(hba1c_empir$HbA1c[hba1c_empir$Therapy=='GLP'])-min(hba1c_empir$HbA1c[hba1c_empir$Therapy=='GLP']))/(100-1))
hba1c_theo$Density  <- hba1c_theo$Probability/(18/(100-1))

statd_therapies	<-	ggplot() + 
						geom_bar(data = hba1c_empir, stat = 'identity', width = 0.01, aes(x = HbA1c, y = Density, fill = 'Empirical'),
							colour = 'black',
							alpha = .5) + 
            geom_line(data = hba1c_theo, aes(x = HbA1c, y = Density, col = 'Markovchart'),
							size = 1.2) + 
						facet_wrap(Therapy ~ . ,labeller = label_bquote(rows = .(Therapy)), ncol = 2) + 
						scale_x_continuous(breaks = seq(4, 22, by = 2)) + 
						xlab('HbA1c') + 
						ylab('Density') + 
						theme_bw() + 
						scale_colour_manual(values = c('#800000')) + 
						scale_fill_manual(values = c('black')) + 
						theme(plot.title = element_text(hjust = 0.5, size = 11), legend.justification = c(1,0),
							legend.position = c(0.99,0.60), legend.title = element_blank(), legend.margin = margin(t = 0, r = 0.1, b = 0, l = 0, unit = 'cm'),
							legend.key = element_rect(fill = 'white', colour = 'white'))

statd_therapies

## ----bothcosts, include = TRUE------------------------------------------------

### Empirical costs for comparison

# ANALOGUE
mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30 +
mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU) +
sampling_cost/mean(deltats[,2])

mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30
mean(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU)
sampling_cost/mean(deltats[,2])

sd(RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR/30 +
RANDOMISED_DATA[grepl("ANALOGUE", RANDOMISED_DATA$THERAPY) &
 !grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU +
sampling_cost/mean(deltats[,2]))

res_ANALOGUE

# GLP
mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30 +
mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU) +
sampling_cost/mean(deltats[,2])

mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR)/30
mean(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU)
sampling_cost/mean(deltats[,2])

sd(RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$THERAPY_COST_EUR/30 +
RANDOMISED_DATA[grepl("GLP", RANDOMISED_DATA$THERAPY),]$COST_CUMU +
sampling_cost/mean(deltats[,2]))

res_GLP


## ----Markov_matrices, echo = TRUE---------------------------------------------
parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
resmtx_ANALOGUE <- Markovchart(
	h = seq(90, 365, by = (365 - 90) / 19),
	k = seq(0.5, 6, by = (6 - 0.5) / 19),
	statdist = stat_ANALOGUE, p = 1,
	constantr = TRUE, ooc_rep = 1,
	cs = sampling_cost,
	coparams = summary(mod.COST)$coef[ , 1],
	crfun = crfun_ANALOGUE,
	crparams = summary(mod.ANALOGUE)$coef[ , 1],
	vcoparams = summary(mod_var.COST)$coef[ , 1],
	vcrfun = vcrfun_ANALOGUE,
	vcrparams = summary(mod_var.ANALOGUE)$coef[ , 1],
	parallel_opt = parall
)

parall <- list(cl=makeCluster(num_workers), forward=FALSE, loginfo=TRUE)
resmtx_GLP <- Markovchart(
	h = seq(90, 365, by = (365 - 90) / 19),
	k = seq(0.5, 6, by = (6 - 0.5) / 19),
	statdist = stat_GLP, p = 1,
	constantr = TRUE, ooc_rep = 1,
	cs = sampling_cost,
	coparams = summary(mod.COST)$coef[ , 1],
	crparams = summary(mod.GLP)$coef[ , 1],
	vcoparams = summary(mod_var.COST)$coef[ , 1],
	vcrparams = summary(mod_var.GLP)$coef[ , 1],
	parallel_opt = parall
)

## ----resmtx_ANALOGUE_fig, echo = TRUE, fig.height = 3, fig.width = 6, warning = FALSE----
resmtx_ANALOGUE[ , 2] <- resmtx_ANALOGUE[ , 2] + 4
plot(x = resmtx_ANALOGUE, y = 'Expected \ndaily cost',
     mid = '#ff9494', high = '#800000',
     xlab = 'Days between samplings',
     ylab = 'Critical HbA1c level')

## ----resmtx_GLP_fig, echo = TRUE, fig.height = 3, fig.width = 6, warning = FALSE----
resmtx_GLP[ , 2] <- resmtx_GLP[ , 2] + 4
plot(x = resmtx_GLP, y = 'Expected \ndaily cost',
		 mid = '#ff9494', high = '#800000',
     xlab = 'Days between samplings',
     ylab = 'Critical HbA1c level')

