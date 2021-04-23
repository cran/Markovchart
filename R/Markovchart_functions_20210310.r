
###GENERALISED_MARKOV_MODEL###

#ANCILLARY FUNCTIONS

qu0_exp.	<-	function(i,h,s,w,V,Vd,delta)
{
	if(i==0)	v	<-	dpois(0,s*h)
	if(i!=0)	v	<-	sum(dpois(w,s*h)*(pgamma((V/(Vd-1))*(i),w,rate=1/delta)-pgamma((V/(Vd-1))*(i-1),w,rate=1/delta)))
	return(v)
}
qu0_exp.	<-	Vectorize(qu0_exp.,"i")

qu0_mix.	<-	function(i,h,s,N,V,Vd,rate,probmix,probnbin,disj)
{
	if(i==0)	v	<-	dpois(0,s*h)
	if(i!=0)	v	<-	sum(dpois(N,s*h)*(mixeddist((V/(Vd-1))*(i),rate,probmix,probnbin,N,disj)-mixeddist((V/(Vd-1))*(i-1),rate,probmix,probnbin,N,disj)))
	return(v)
}
qu0_mix.	<-	Vectorize(qu0_mix.,"i")

qu_exp.	<-	function(i,h,s,w,delta,int2)
{
	if(i==1)	v	<-	dpois(0,s*h) + sum(dpois(w,s*h)*pgamma(int2[i],w,rate=1/delta))
	if(i!=1)	v	<-	sum(dpois(w,s*h)*(pgamma(int2[i],w,rate=1/delta)-pgamma(int2[i-1],w,rate=1/delta)))
	return(v)
}
qu_exp.	<-	Vectorize(qu_exp.,"i")

qu_mix.	<-	function(i,h,s,N,rate,probmix,probnbin,disj,int2)
{
	if(i==1)	v	<-	dpois(0,s*h) + sum(dpois(N,s*h)*mixeddist(int2[i],rate,probmix,probnbin,N,disj))
	if(i!=1)	v	<-	sum(dpois(N,s*h)*(mixeddist(int2[i],rate,probmix,probnbin,N,disj)-mixeddist(int2[i-1],rate,probmix,probnbin,N,disj)))
	return(v)
}
qu_mix.	<-	Vectorize(qu_mix.,"i")

jloop.	<-	function(j,nov,quv,mtxvec,index)
{
  res=as.vector(c(NA,NA))
  nov=nov
  quv=quv
  mtxvec=mtxvec
  index=index
  kvec          <-  quv[j-((1:min(index,j))-1)] * mtxvec[1:min(index,j)]
  res[1]		<-	sum((1-nov[j])	* kvec)
  res[2]		<-	sum(nov[j]		* kvec)
  return(res)
}
jloop.	<-	Vectorize(jloop.,"j")

roundany	<-	function(x, accuracy, fun=round){fun(x/accuracy)*accuracy}

#P(gamma + negbin < z) when the number of shifts is given
pgamnbin	<-	function(z,shape,rate,prob,n,disj)
{
	n=n
	z=z
	shape=shape
	rate=rate
	prob=prob
	if(roundany(z,disj,floor)<((n-shape)*disj))	res	<-	0
	else{
		if(z<disj | (n-shape)==0)						res	<-	pgamma(z,shape=n,rate=rate)
		else{ 
														res	<-	sum(pgamma(z-seq((n-shape)*disj,roundany(z,disj,floor),disj),shape=shape,rate=rate)*dnbinom(0:((roundany(z,disj,floor)-(n-shape)*disj)/disj),size=n-shape,prob=prob))
		}
	}
	return(res)
}
pgamnbin	<-	Vectorize(pgamnbin,"shape")

#This is the mixture distribution of the distance from mu_0 (the convolution of the mixture distribution of single shifts)
mixeddist	<-	function(z,rate,probmix,probnbin,N,disj)
{
	N=N
	z=z
	rate=rate
	probmix=probmix
	probnbin=probnbin
	disj=disj
	sum(dbinom(N-(0:N),N,prob=probmix)*pgamnbin(z=z,shape=0:N,rate=rate,prob=probnbin,n=N,disj))
}
mixeddist	<-	Vectorize(mixeddist,"N")

#DEFAULT COST FUNCTIONS

crfun_default	<-	function(mudist,crparams)
{
  mudist=mudist
  crb=crparams[1]
  crs=crparams[2]
  
  cr		<-	crb + crs*mudist
  return(cr)
}
crfun_default	<-	Vectorize(crfun_default,"mudist")

vcrfun_default	<-	function(mudist,vcrparams)
{
  mudist=mudist
  vcrb=vcrparams[1]
  vcrs=vcrparams[2]
  
  vcr		<-	vcrb + vcrs*mudist
  return(vcr)
}
vcrfun_default	<-	Vectorize(vcrfun_default,"mudist")

cofun_default	<-	function(sqmudist,coparams)
{
  sqmudist=sqmudist
  cob=coparams[1]
  cos=coparams[2]
  
  co		<-	cob + cos*sqmudist
  return(co)
}
cofun_default	<-	Vectorize(cofun_default,"sqmudist")

vcofun_default	<-	function(sqmudist,vcoparams)
{
  sqmudist=sqmudist
  vcob=vcoparams[1]
  vcos=vcoparams[2]
  
  vco		<-	vcob + vcos*sqmudist
  return(vco)
}
vcofun_default	<-	Vectorize(vcofun_default,"sqmudist")

#S3 METHODS

print.Markov_stationary <- function(x,...) {
	print(x[1],...)
}

print.Markov_chart <- function(x,...) {
	print(x[1:3],...)
}

print.Markov_sim <- function(x,...) {
	print(str(x),...)
}

plot.Markov_grid <- function(x,y=expression(atop(italic("G")*-value~per, unit~time)),xlab="Time between samplings",ylab="Critical value",low="white",mid="#999999",high="black",colour="white",...)
{
	if(!inherits(x, "Markov_grid") | dim(x)[2]!=3)	stop("x should be a Markov_grid data.frame with three columns (preferably created by the Markovchart function): time between samplings, critical value and the weighted mean of the expected cost and the cost standard deviation.")
	colnames(x)	<-	c("h","k","value")
	
	ggplot(x, aes(h, k, z=value)) + 
		geom_raster(aes(fill = value), interpolate=TRUE) + 
		scale_fill_gradient2(low=low, mid=mid, high=high, midpoint=median(x$value), name=y) + 
		geom_contour(colour = colour) + 
		geom_text_contour(size = 3.5, stroke = 0.1, breaks = pretty(x[,"value"],20)) + 
		xlab(xlab) +
		ylab(ylab) +
		theme_minimal()
}


#STATIONARY DISTRIBUTION

Markovstat	<-	function(shiftfun=c("exp","exp-geo","deg"),h,k,sigma,s,delta,probmix=1,probnbin=0.5,disj=1,
						 RanRep=FALSE,alpha=NULL,beta=NULL,RanSam=FALSE,StateDep=FALSE,a=NULL,b=NULL,q=NULL,z=NULL,
						 Vd=100,V=NULL,Qparam=30)
{
	#PARAMETERS
	if(h<=0| k<0)							stop("Non-applicable h or k parameters were given. h is the time between samplings and k is the critical value. Both should be positive numbers (k is allowed to be 0, but not negative, since positive shifts are assumed).")
	if(sigma<=0 | s<=0 | delta<=0)			stop("Non-applicable sigma, s or delta parameters were given. sigma is the process variance, s is the expected number of shifts in a sampling interval and delta is the exponential shift size parameter. These should all be positive numbers.")
	if(shiftfun=="deg" & (RanRep | RanSam))	warning("Random repair and random sampling is not implemented for degenerate ('deg') shift size distribution. Values of RanRep and RanSam are ignored.")


	if(Vd!=round(Vd))
	{
		Vd		<-	round(Vd)
		warning("The Vd discretisation parameter was rounded, as not an integer was given.")
	}

	if(Qparam!=round(Qparam))
	{
		Qparam	<-	round(Qparam)
		warning("The Qparam discretisation parameter was rounded, as not an integer was given.")
	}
	
	############STATIONARY DISTRIBUTION CALCULATION############
	
	paramlist <- list(h=h,k=k,shiftfun=shiftfun,sigma=sigma,s=s,delta=delta,probmix=probmix,probnbin=probnbin,disj=disj,
					  RanRep=RanRep,alpha=alpha,beta=beta,RanSam=RanSam,StateDep=StateDep,a=a,b=b,q=q,z=z,
					  Vd=Vd,V=V,Qparam=Qparam)
	
	if(shiftfun=="exp" | shiftfun=="exp-geo")
	{
		if(V<=0)											stop("Non-applicable V parameter was given. V is the maximum distance from the target value taken into account. This should be a positive number.")
		if(RanRep & (is.null(alpha) | is.null(beta))) 		stop("Parameters alpha or beta are missing. If RanRep=TRUE, then the beta distribution parameters (alpha and beta) must be given.")
		if(RanSam & StateDep & (is.null(a) | is.null(b)))	stop("Parameters a or b are missing. If RanRep=TRUE and StateDep=TRUE, then the beta distribution parameters (a and b) must be given.")
		if(RanSam & !StateDep & (is.null(q) | is.null(z)))	stop("Parameters q or z are missing. If RanRep=TRUE and StateDep=FALSE, then the logistic function parameters (q and z) must be given.")
		if(Vd<3)											stop("Non-applicable Vd parameter was given. This is a discretisation parameter (number of states after disretisation), which should be an integer value greater than 2.")
		if(probmix<0 | probmix>1)							stop("Invalid probmix parameter. probmix is the weight of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1.")
		if(probnbin<0 | probnbin>1)							stop("Invalid probnbin parameter. probnbin is the probability parameter of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1.")
		if(Qparam<1)										stop("Non-applicable Qparam parameter was given. This is the number of shifts considered within a sampling interval, which should be an integer value greater than 0.")

	
		if(is.null(V))	stop("Value of V is NULL. V must be given as a positive number in case of random shift size (i.e. when shiftfun='exp-geo' or shiftfun='exp').")

		int2	<-	seq(0,V,by=(V/(Vd-1))) + (V/(Vd-1))*0.5
	
		if(shiftfun=="exp")
		{
			w	<-	1:Qparam
			#Evaluation of function
			quv0				<-	qu0_exp.(0:(Vd-1),h,s,w,V,Vd,delta)
			#Evaluation of function
			quv					<-	qu_exp.(1:Vd,h,s,w,delta,int2)
		}
		
		if(shiftfun=="exp-geo")
		{
			if(disj<=0)	stop("Non-applicable disj parameter was given. disj is the size of a discrete jump in case of exponential-geometric mixture shift distribution and should be a positive number.")

			rate    <- 1/delta
			N	    <- 1:Qparam
			#Evaluation of function
			quv0				<-	qu0_mix.(0:(Vd-1),h,s,N,V,Vd,rate,probmix,probnbin,disj)
			#Evaluation of function
			quv					<-	qu_mix.(1:Vd,h,s,N,rate,probmix,probnbin,disj,int2)
		}
		
		quv0[length(quv0)]	<-	quv0[length(quv0)] + (1-sum(quv0))
		quv[length(quv)]	<-	quv[length(quv)] + (1-sum(quv))
	
		nov		<-	pnorm(k-((V/(Vd-1))*(0:(Vd-1))-((V/(Vd-1))*0.5)),0,sigma)
		nov[1]	<-	pnorm(k,0,sigma)
		
		if(RanRep)
		{
			mtx	<-	NULL
			for(j in seq(1,Vd,1))
			{
				v	  <-	NULL
				v	  <-	c(pbeta(((0:j)+1)/(j-0.5),alpha,beta)-pbeta((0:j)/(j-0.5),alpha,beta),rep(0,Vd-j))
				mtx	<-	rbind(mtx,v)
			}
			mtx	<-	mtx[,1:Vd]
			mtx	<-	rbind(c(1,rep(0,Vd)), cbind(0,mtx)[1:(Vd-1),])[,1:Vd]
		} else
		{
			mtx	<-	cbind(1,matrix(0,Vd,Vd-1))
		}
		
		###
		
		mx_alarm	<-	matrix(numeric(0),Vd,Vd)
		mx_ooc		<-	matrix(numeric(0),Vd,Vd)
		mx_alarm_o	<-	matrix(numeric(0),Vd,Vd)
		mx_ooc_o	<-	matrix(numeric(0),Vd,Vd)
	
		ione        <-  sum(mtx[1,]!=0)+1
		kvec_disc1	<-	numeric(min(ione,Vd))
		kvec_disc2	<-	numeric(min(ione,Vd))
		for(m in 1:ione)
		{
		Vdindex       <-  Vd-(m-1)
			quv_w		<-	c(quv0[1:(Vdindex-1)], quv0[Vdindex]+(1-sum(quv0[1:Vdindex])))
			kvec_disc1[m]	<-	(1-nov[Vd])*	quv_w[Vdindex]
			kvec_disc2[m]	<-	nov[Vd]*	  	quv_w[Vdindex]
		}
		row_alarm	<-	numeric(Vd)
		row_ooc		<-	numeric(Vd)
		jloopres  	<-  jloop.(1:(Vd-1),nov,quv,mtxvec=mtx[1,],index=ione)
		row_alarm[1:(Vd-1)]	<-	jloopres[1,]
		row_ooc[1:(Vd-1)]	<-	jloopres[2,]
		row_alarm[Vd]	  <-	sum(kvec_disc1*mtx[1,1:min(ione,Vd)])
		row_ooc[Vd]		  <-	sum(kvec_disc2*mtx[1,1:min(ione,Vd)])
		mx_alarm[1,]	  <-	row_alarm
		mx_ooc[1,]		  <-	row_ooc
		mx_alarm_o[1,]	  <-	c(quv0[0:(Vd-1)],quv0[Vd]+(1-sum(quv0[0:(Vd)])))*(1-nov)
		mx_ooc_o[1,]	  <-	c(quv0[0:(Vd-1)],quv0[Vd]+(1-sum(quv0[0:(Vd)])))*nov
	
	
		for(i in 2:Vd)
		{
		inum        <-  sum(mtx[i,]!=0)+1
		mtxvec      <-  mtx[i,]
			kvec_disc1	<-	numeric(min(inum,Vd))
			kvec_disc2	<-	numeric(min(inum,Vd))
			for(m in 1:inum)
			{
			Vdindex       	<-  Vd-(m-1)
				quv_w			<-	c(quv[1:(Vdindex-1)], quv[Vdindex]+(1-sum(quv[1:Vdindex])))
				kvec_disc1[m]	<-	(1-nov[Vd])*	quv_w[Vdindex]
				kvec_disc2[m]	<-	nov[Vd]*	  	quv_w[Vdindex]
			}
			row_alarm	<-	numeric(Vd)
			row_ooc		<-	numeric(Vd)
			jloopres  	<-  jloop.(1:(Vd-1),nov,quv,mtxvec=mtxvec,index=inum)
			row_alarm[1:(Vd-1)]	<-	jloopres[1,]
			row_ooc[1:(Vd-1)]	<-	jloopres[2,]
			row_alarm[Vd]	  <-	sum(kvec_disc1*mtx[i,1:min(inum,Vd)])
			row_ooc[Vd]		  <-	sum(kvec_disc2*mtx[i,1:min(inum,Vd)])
			mx_alarm[i,]	  <-	row_alarm
			mx_ooc[i,]		  <-	row_ooc
			mx_alarm_o[i,]	<-	c(rep(0,(i-1)),quv[0:(Vd-(i-1)-1)],quv[Vd-(i-1)]+(1-sum(quv[0:(Vd-(i-1))])))*(1-nov)
			mx_ooc_o[i,]  	<-	c(rep(0,(i-1)),quv[0:(Vd-(i-1)-1)],quv[Vd-(i-1)]+(1-sum(quv[0:(Vd-(i-1))])))*nov
		}
	
		mx_f	<-	cbind(rbind(mx_alarm,mx_alarm_o),rbind(mx_ooc,mx_ooc_o))
		mx_fin	<-	cbind(mx_f[,1], mx_f[,(Vd+1)], mx_f[,2:Vd], mx_f[,(Vd+2):(Vd*2)])
		mx_fin	<-	rbind(mx_fin[1,], mx_fin[(Vd+1),], mx_fin[2:Vd,], mx_fin[(Vd+2):(Vd*2),])
		mx_fin	<-	mx_fin*(1/apply(mx_fin,1,sum))
		
		###
		
		#RSI
		if(RanSam)
		{
			if(StateDep)
			{
				viv	<-	pbeta(((1:Vd)-1/2)/Vd,a/h,b)
			} else
			{
				viv	<-	rep(1/(1+exp(-q*(h-z))), Vd)
			}
		} else
		{
			viv	<-	rep(1,Vd)
		}
		
		mx_fin_n	<-	sweep(mx_fin[,2:(Vd+1)],MARGIN=2,viv,`*`)
		mx_fin_n	<-	cbind((mx_fin[,c(1,(Vd+2):(Vd*2))] + sweep(mx_fin[,2:(Vd+1)],MARGIN=2,1-viv,`*`))[,1],mx_fin_n)
		mx_fin_n	<-	cbind(mx_fin_n,(mx_fin[,c(1,(Vd+2):(Vd*2))] + sweep(mx_fin[,2:(Vd+1)],MARGIN=2,1-viv,`*`))[,2:Vd])
		
		mx_fin	<-	mx_fin_n
		mx_fin	<-	mx_fin*(1/apply(mx_fin,1,sum))
		
		###
		mx_fin_full	<-	mx_fin
		if(RanRep==TRUE)	mx_fin	<-	mx_fin[3:nrow(mx_fin),3:ncol(mx_fin)]
		
		eigenvec        <-  as.numeric(eigen(t(mx_fin))$vectors[,1])
		statd			<-	eigenvec/sum(eigenvec)
		statd[statd<0]	<-	0
		if(RanRep==TRUE)	statd	<-	c(0,0,statd)
		names(statd)	<-	c("In-control","False-alarm",paste("Out-of-control",1:(Vd-1),c(round(int2[1:(Vd-2)],3),paste(round(V/(Vd-1)*(Vd-2),3),"+",sep="")),sep="_"),paste("True-alarm",1:(Vd-1),c(round(int2[1:(Vd-2)],3),paste(round(V/(Vd-1)*(Vd-2),3),"+",sep="")),sep="_"))	
		
		paramlist <- c(paramlist, list(mtx=mtx), list(viv=viv))
	}
	
	#####
	if(shiftfun=="deg")
	{
		mx_fin	<-	t(matrix(
		c((1-pexp(h,s))*pnorm(k,0,sigma),	(1-pexp(h,s))*(1-pnorm(k,0,sigma)),	pexp(h,s)*pnorm(k-delta,0,sigma),	pexp(h,s)*(1-pnorm(k-delta,0,sigma)),
		(1-pexp(h,s))*pnorm(k,0,sigma),	(1-pexp(h,s))*(1-pnorm(k,0,sigma)),	pexp(h,s)*pnorm(k-delta,0,sigma),	pexp(h,s)*(1-pnorm(k-delta,0,sigma)),
									0,						   			 0,			  pnorm(k-delta,0,sigma),				1-pnorm(k-delta,0,sigma),
		(1-pexp(h,s))*pnorm(k,0,sigma),	(1-pexp(h,s))*(1-pnorm(k,0,sigma)),	pexp(h,s)*pnorm(k-delta,0,sigma),	pexp(h,s)*(1-pnorm(k-delta,0,sigma))),
		4,4))
	
		mx_fin_full		<-	mx_fin
		eigenvec        <-  as.numeric(eigen(t(mx_fin))$vectors[,1])
		statd			<-	eigenvec/sum(eigenvec)
		statd[statd<0]	<-	0
		names(statd)	<-	c("In-control","False-alarm","Out-of-control","True-alarm")
	}
	
	mx_fin_full				<-	as.data.frame(mx_fin_full)
	colnames(mx_fin_full)	<-	names(statd)
	rownames(mx_fin_full)	<-	names(statd)

	############STATIONARY DISTRIBUTION CALCULATION END############
	
	results			<-	list(Stationary_distribution=statd, Transition_matrix=mx_fin_full, Param_list=paramlist)
	class(results)	<-	c("Markov_stationary", class(results))
		
	return(results)
}


#COST CALCULATION
costfunexp	<-	function(fp=NULL,statdist,p=1,constantr=FALSE,ooc_rep=0,cs=NULL,cofun=cofun_default,
						 coparams=NULL,crfun=crfun_default,crparams=NULL,vcofun=vcofun_default,vcoparams=c(0,0),
						 vcrfun=vcrfun_default,vcrparams=c(0,0),ALL=TRUE,Markovstat_fun=Markovstat)
{
	#PARAMETERS
	inherited                           <- c(statdist[[3]],statd=statdist[1])
	names(inherited)[length(inherited)] <- "statd"
	list2env(inherited,environment())
	
	if(p<0 | p>1)		stop("Invalid p parameter. p is the weight of the cost expectation in the calculation of the G-value and should be between 0 and 1.")
	if(shiftfun!="exp")	stop("shiftfun is not 'exp', i.e., the stationary distribution is not compatible with the cost calculation method.")

	if(!is.null(fp)) {
		h=fp[1]					#(Prescribed) time between samplings
		k=fp[2]					#Critical value
		
		actualstatd	<-	Markovstat_fun(h=h,k=k,shiftfun="exp",sigma=sigma,s=s,delta=delta,probmix=probmix,probnbin=probnbin,disj=disj,
									   RanRep=RanRep,alpha=alpha,beta=beta,RanSam=RanSam,StateDep=StateDep,a=a,b=b,q=q,z=z,
									   Vd=Vd,V=V,Qparam=Qparam)
		inherited                           <- c(actualstatd[[3]],statd=actualstatd[1])
		names(inherited)[length(inherited)] <- "statd"
		list2env(inherited,environment())
	}

	INC	<-	c(1,rep(0,(Vd*2-1)))
	FA	<-	c(1,rep(0,(Vd*2-1)))
	TA	<-	cbind(mtx[2:Vd,1],matrix(0,Vd-1,Vd),mtx[2:Vd,2:Vd])
	OOC	<-	cbind(matrix(0,Vd-1,Vd+1),diag(1,Vd-1,Vd-1))
	
	MintaUtanMtx	<-	rbind(INC,FA,TA,OOC)
	MintaUtanMtx	<-	cbind(MintaUtanMtx[,1] + MintaUtanMtx[,2], MintaUtanMtx[,3:(Vd+1)] + MintaUtanMtx[,(Vd+2):(Vd*2)])
	
	#
	
	int			<-	seq(0,V,by=(V/(Vd-1))) - (V/(Vd-1))/2
	int[1]		<-	0		
	allsqmeans	<-	s*h*delta*(((s*h*delta)/3)+delta+int)+int^2
	Averages	<-	as.vector(apply(sweep(MintaUtanMtx,MARGIN=2,allsqmeans,`*`),1,sum))
	
	distance_dist	<-	statd[c(1,(Vd+2):(Vd*2))] + statd[2:(Vd+1)]
	
	cr			<-	crfun(int,crparams)
	cr[cr<0]	<-	0
	cr			<-	c(cr[1]*ooc_rep,cr[1],cr[2:Vd],cr[2:Vd]*ooc_rep)
	vcr			<-	vcrfun(int,vcrparams)
	vcr[vcr<0]	<-	0
	vcr			<-	c(vcr[1]*ooc_rep,vcr[1],vcr[2:Vd],vcr[2:Vd]*ooc_rep)
	if(!constantr)
	{
		cr					<-	cr/h
		vcr					<-	vcr/h^2
	}
	co			<-	cofun(Averages,coparams)
	vco			<-	vcofun(Averages,vcoparams)
	viv			<-	c(viv[1],viv[1],viv[2:length(viv)],viv[2:length(viv)])
	
	values		<-	(cs*(1/h)*viv) + co + cr
	if(sum(values<0)>0)	warning("Negative cost(s) were calculated from the given cost function and parameters. This implies income instead of expenses.")
	
	varvec		<-	vco + vcr
	if(sum(varvec==Inf)>0){
		varvec[varvec==Inf]	<-	max(varvec[varvec!=Inf])
		warning("Infinite repair/out-of-control cost variance was calculated. Substituting with the second largest variance in the vector.")
	}
	
	#PROCESS VARIANCE
	VARC				<-	sum(values^2*statd)-sum(values*statd)^2
	if(VARC<0) {
		VARC	<-	0
		warning("Negative process variance was estimated, check the model parameters.")
	}
	#TOTAL VARIANCE
	ALLSDC	<-	sqrt(sum(varvec * statd) + VARC)
	#AVERAGE COST
	EC		<-	sum(values*statd)
	
	mom2				<-	sum(values^2*statd)
	mom3				<-	sum(values^3*statd)
	mom4				<-	sum(values^4*statd)
	
	G	<-	p*EC + (1-p)*ALLSDC
	
	######
	
	results	<-	list(Results=as.data.frame(t(c(G,EC,ALLSDC,sqrt(VARC),mom2,mom3,mom4))),Subcosts=as.data.frame(t(c(cs/h,sum(cr*statd),sum(co*statd)))),Parameters=as.data.frame(t(c(h,k))),Stationary_distribution=statd)
	colnames(results[[1]])	<-	c("G-value","Expected cost (C)","Total cost std. dev.","Cost std. dev. due to process var.","Second process moment","Third process moment","Fourth process moment")
	colnames(results[[2]])	<-	c("Sampling cost","Repair cost","OOC cost")
	colnames(results[[3]])	<-	c("Time between samplings (h)","Critical value (k)")
	class(results)			<-	c("Markov_chart", class(results))
	
	if(!ALL)	return(G)
	if(ALL)		return(results)

}

costfunexpgeo	<-	function(fp=NULL,statdist,p=1,constantr=FALSE,ooc_rep=0,cs=NULL,cofun=cofun_default,
							 coparams=NULL,crfun=crfun_default,crparams=NULL,vcofun=vcofun_default,vcoparams=c(0,0),
							 vcrfun=vcrfun_default,vcrparams=c(0,0),ALL=TRUE,Markovstat_fun=Markovstat,qu_mix=qu_mix.)
{
	#PARAMETERS
	inherited                           <- c(statdist[[3]],statd=statdist[1])
	names(inherited)[length(inherited)] <- "statd"
	list2env(inherited,environment())
	
	if(p<0 | p>1)			stop("Invalid p parameter. p is the weight of the cost expectation in the calculation of the G-value and should be between 0 and 1.")
	if(shiftfun!="exp-geo")	stop("shiftfun is not 'exp-geo', i.e., the stationary distribution is not compatible with the cost calculation method.")
	
	if(!is.null(fp)) {
		h=fp[1]					#(Prescribed) time between samplings
		k=fp[2]					#Critical value
		
		actualstatd	<-	Markovstat_fun(h=h,k=k,shiftfun="exp-geo",sigma=sigma,s=s,delta=delta,probmix=probmix,probnbin=probnbin,disj=disj,
									   RanRep=RanRep,alpha=alpha,beta=beta,RanSam=RanSam,StateDep=StateDep,a=a,b=b,q=q,z=z,
									   Vd=Vd,V=V,Qparam=Qparam)
		inherited                           <- c(actualstatd[[3]],statd=actualstatd[1])
		names(inherited)[length(inherited)] <- "statd"
		list2env(inherited,environment())
	}
	
	int2	<-	seq(0,V,by=(V/(Vd-1))) + (V/(Vd-1))*0.5
	rate	<-	1/delta
	N		<-	1:Qparam
	
	INC	<-	c(1,rep(0,(Vd*2-1)))
	FA	<-	c(1,rep(0,(Vd*2-1)))
	TA	<-	cbind(mtx[2:Vd,1],matrix(0,Vd-1,Vd),mtx[2:Vd,2:Vd])
	OOC	<-	cbind(matrix(0,Vd-1,Vd+1),diag(1,Vd-1,Vd-1))
	
	MintaUtanMtx	<-	rbind(INC,FA,TA,OOC)
	MintaUtanMtx	<-	cbind(MintaUtanMtx[,1] + MintaUtanMtx[,2], MintaUtanMtx[,3:(Vd+1)] + MintaUtanMtx[,(Vd+2):(Vd*2)])
	
	#
	
	int			<-	seq(0,V,by=(V/(Vd-1))) - (V/(Vd-1))/2
	int[1]		<-	0
	
	qvec<-NULL
	for(delt in seq(0,h,by=h/9))
	{
		quv					<-	qu_mix(1:Vd,h=delt,s,N,rate,probmix,probnbin,disj,int2)
		quv[length(quv)]	<-	quv[length(quv)] + (1-sum(quv))
		qvec				<-	rbind(qvec,quv)
	}
	
	sqdistances	<-	(t(matrix(int , length(int) , length(int))) + int)^2
	
	meansatt<-NULL
	for (i in 1:nrow(qvec))
	{
		meansatt	<-	rbind(meansatt,apply(sqdistances*qvec[i,],2,sum))
	}
	allsqmeans	<-	apply(meansatt,2,sum)/10
	
	Averages	<-	as.vector(apply(sweep(MintaUtanMtx,MARGIN=2,allsqmeans,`*`),1,sum))
	
	distance_dist	<-	statd[c(1,(Vd+2):(Vd*2))] + statd[2:(Vd+1)]
	
	cr			<-	crfun(int,crparams)
	cr[cr<0]	<-	0
	cr			<-	c(cr[1],cr[1],cr[2:Vd],cr[2:Vd])
	vcr			<-	vcrfun(int,vcrparams)
	vcr[vcr<0]	<-	0
	vcr			<-	c(vcr[1],vcr[1],vcr[2:Vd],vcr[2:Vd])
	if(!constantr)
	{
		cr[-(3:(Vd+1))]		<-	0
		cr					<-	cr/h
		vcr[-(3:(Vd+1))]	<-	0
		vcr					<-	vcr/h^2
	}
	co			<-	cofun(Averages,coparams)
	vco			<-	vcofun(Averages,vcoparams)
	viv			<-	c(viv[1],viv[1],viv[2:length(viv)],viv[2:length(viv)])
	
	values		<-	(cs*(1/h)*viv) + co + cr
	if(sum(values<0)>0)	warning("Negative cost(s) were calculated from the given cost function and parameters. This implies income instead of expenses.")
	
	varvec		<-	vco + vcr
	if(sum(varvec==Inf)>0){
		varvec[varvec==Inf]	<-	max(varvec[varvec!=Inf])
		warning("Infinite repair/out-of-control cost variance was calculated. Substituting with the second largest variance in the vector.")
	}
	
	#PROCESS VARIANCE
	VARC				<-	sum(values^2*statd)-sum(values*statd)^2
	if(VARC<0) {
		VARC	<-	0
		warning("Negative process variance was estimated, check the model parameters.")
	}
	#TOTAL VARIANCE
	ALLSDC	<-	sqrt(sum(varvec * statd) + VARC)
	#AVERAGE COST
	EC		<-	sum(values*statd)
	
	mom2				<-	sum(values^2*statd)
	mom3				<-	sum(values^3*statd)
	mom4				<-	sum(values^4*statd)
	
	G	<-	p*EC + (1-p)*ALLSDC
	
	######
	
	results	<-	list(Results=as.data.frame(t(c(G,EC,ALLSDC,sqrt(VARC),mom2,mom3,mom4))),Subcosts=as.data.frame(t(c(cs/h,sum(cr*statd),sum(co*statd)))),Parameters=as.data.frame(t(c(h,k))),Stationary_distribution=statd)
	colnames(results[[1]])	<-	c("G-value","Expected cost (C)","Total cost std. dev.","Cost std. dev. due to process var.","Second process moment","Third process moment","Fourth process moment")
	colnames(results[[2]])	<-	c("Sampling cost","Repair cost","OOC cost")
	colnames(results[[3]])	<-	c("Time between samplings (h)","Critical value (k)")
	class(results)			<-	c("Markov_chart", class(results))
	
	if(!ALL)	return(G)
	if(ALL)		return(results)
}


costfundeg	<-	function(fp=NULL,statdist,cs=NULL,crparams=NULL,cf=crparams,coparams=NULL,p=1,ALL=TRUE,Markovstat_fun=Markovstat)
{
	#PARAMETERS
	inherited                           <- c(statdist[[3]],statd=statdist[1])
	names(inherited)[length(inherited)] <- "statd"
	list2env(inherited,environment())
	
	if(p<0 | p>1)		stop("Invalid p parameter. p is the weight of the cost expectation in the calculation of the G-value and should be between 0 and 1.")
	if(shiftfun!="deg")	stop("shiftfun is not 'deg', i.e., the stationary distribution is not compatible with the cost calculation method.")
	
	if(!is.null(fp)) {
		h=fp[1]					#(Prescribed) time between samplings
		k=fp[2]					#Critical value
		
		actualstatd	<-	Markovstat_fun(h=h,k=k,shiftfun="deg",sigma=sigma,s=s,delta=delta,probmix=probmix,probnbin=probnbin,disj=disj,
									  RanRep=RanRep,alpha=alpha,beta=beta,RanSam=RanSam,StateDep=StateDep,a=a,b=b,q=q,z=z,
									  Vd=Vd,V=V,Qparam=Qparam)
		inherited                           <- c(actualstatd[[3]],statd=actualstatd[1])
		names(inherited)[length(inherited)] <- "statd"
		list2env(inherited,environment())
	}
	
	Tau=(1-(1+h*s)*exp(sigma)^(-h*s))/(s*(1-exp(sigma)^(-h*s)))
	B=(h-Tau)/h
	
	values	<-	c(cs/h,(cs/h)+(cf/h), (cs/h)+coparams, (cs/h)+(coparams*B+crparams/h))
	if(sum(values<0)>0)	warning("Negative cost(s) were calculated from the given cost function and parameters. This implies income instead of expenses.")
	
	#PROCESS VARIANCE
	VARC				<-	sum(values^2*statd)-sum(values*statd)^2
	if(VARC<0) {
		VARC	<-	0
		warning("Negative process variance was estimated, check the model parameters.")
	}
	
	#AVERAGE COST
	EC		<-	sum(values*statd)
	
	mom2				<-	sum(values^2*statd)
	mom3				<-	sum(values^3*statd)
	mom4				<-	sum(values^4*statd)
	
	G	<-	p*EC + (1-p)*sqrt(VARC)
	
	results					<-	list(Results=as.data.frame(t(c(G,EC,sqrt(VARC),mom2,mom3,mom4))),Subcosts=as.data.frame(t(c(cs/h,statd[2]*(cf/h),statd[3]*coparams,statd[4]*(coparams*B+crparams/h)))),Parameters=as.data.frame(t(c(h,k))),Stationary_distribution=statd)
	colnames(results[[1]])	<-	c("G-value","Expected cost (C)","Cost std. dev. due to process var.","Second process moment","Third process moment","Fourth process moment")
	colnames(results[[2]])	<-	c("In-control cost","False-alarm cost","Out-of-control cost","True-alarm cost")
	colnames(results[[3]])	<-	c("Time between samplings (h)","Critical value (k)")
	class(results)			<-	c("Markov_chart", class(results))
	
	if(!ALL)	return(G)
	if(ALL)		return(results)
}


Markovchart	<-	function(statdist,h=NULL,k=NULL,OPTIM=FALSE,p=1,
						 constantr=FALSE,ooc_rep=0,cs=NULL,cofun=cofun_default,coparams=NULL,crfun=crfun_default,crparams=NULL,
						 cf=crparams, vcofun=vcofun_default,vcoparams=c(0,0),vcrfun=vcrfun_default,vcrparams=c(0,0),
						 method=c("L-BFGS-B","Nelder-Mead","BFGS","CG","SANN","Brent"),parallel_opt=NULL,silent=TRUE,...)
{
	list2env(list(...),environment())
	if(is.null(h)) h <- statdist[[3]]$h
	if(is.null(k)) k <- statdist[[3]]$k
	shiftfun <- statdist[[3]]$shiftfun

	if(sum(h<=0)>0 | sum(k<0)>0)							stop("Non-applicable h or k parameters were given. h is the time between samplings and k is the critical value. Both should be positive numbers (k is allowed to be 0, but not negative, since positive shifts are assumed).")
	if((length(h)>1 | length(k)>1) & OPTIM)					warning("h or k is of length greater than one. Output will be a data.frame of G-values for all given h and k parameter values, OPTIM value as TRUE is ignored.")
	if(is.null(cs) | is.null(coparams) | is.null(crparams))	stop("Some necessary cost parameters are missing. Check the documentation for the list of necessary parameters.")

	if(OPTIM){
		if(!exists("gr",where=environment()))		gr 		<-	NULL
		if(!exists("lower",where=environment()))	lower	<-	c(0.00000000001,0)
		if(!exists("upper",where=environment()))	upper	<-	Inf
		if(!exists("control",where=environment()))	control	<-	list()
		if(!exists("hessian",where=environment()))	hessian	<-	FALSE
		
		if(length(lower)==1) {
			if(lower<=0) {
			lower		<-	c(0.00000000001,0)
				warning("The value of parameter lower is adjusted to c(0.00000000001,0). Parameter h must be positive and k must be non-negative.")
			}
		} else {
			if(lower[1]<=0) {
				lower[1]	<-	0.00000000001
				warning("The value of parameter lower[1] is adjusted to 0.00000000001. Parameter h must be positive.")
			}
			
			if(lower[2]<0) {
				lower		<-	0
				warning("The value of parameter lower[2] is adjusted to 0. Parameter k must be non-negative.")
			}
		}
	}

	if(shiftfun=="exp"){
		if(length(h)>1 | length(k)>1) {
			if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
			registerDoParallel(parallel_opt$cl)
			Gmtx	<-	foreach(i=h, .combine='rbind') %:%
							foreach(j=k, .combine='rbind') %dopar%{
								c(i,j,costfunexp(fp=c(i,j),
											statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											vcrfun=vcrfun,vcrparams=vcrparams,ALL=FALSE)
									)
			}
			stopCluster(parallel_opt$cl)
			
			Gmtx			<-	as.data.frame(Gmtx)
			rownames(Gmtx)	<-	NULL
			colnames(Gmtx)	<-	c("h","k","value")
			res				<-	Gmtx
			class(res)		<-	c("Markov_grid", class(res))
		
		}	else {
			if(OPTIM){
				if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
				res_optim	<-	optimParallel(par=c(h,k), fn=costfunexp, gr=gr, lower=lower, upper=upper, method=method, control=control,
										      hessian=hessian, parallel=parallel_opt, statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
										      coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
										      vcrfun=vcrfun,vcrparams=vcrparams,ALL=FALSE)
				stopCluster(parallel_opt$cl)
				
				res			<-	costfunexp(fp=res_optim$par,
										   statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
										   coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
										   vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
			}	else {
				if(sum(h!=statdist[[3]]$h)>0 | sum(k!=statdist[[3]]$k)>0){
					res			<-	costfunexp(fp=c(h,k),statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											   coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											   vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
				}	else {
					res			<-	costfunexp(statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											   coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											   vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
				}
			}
		}
	}
	if(shiftfun=="exp-geo"){
		if(length(h)>1 | length(k)>1) {
			if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
			registerDoParallel(parallel_opt$cl)
			Gmtx	<-	foreach(i=h, .combine='rbind') %:%
							foreach(j=k, .combine='rbind') %dopar%{
								c(i,j,costfunexpgeo(fp=c(i,j),
											statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											vcrfun=vcrfun,vcrparams=vcrparams,ALL=FALSE)
									)
			}
			stopCluster(parallel_opt$cl)
			
			Gmtx			<-	as.data.frame(Gmtx)
			rownames(Gmtx)	<-	NULL
			colnames(Gmtx)	<-	c("h","k","value")
			res				<-	Gmtx
			class(res)		<-	c("Markov_grid", class(res))
		
		}	else {
			if(OPTIM){
				if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
				res_optim	<-	optimParallel(par=c(h,k), fn=costfunexpgeo, gr=gr, lower=lower, upper=upper, method=method, control=control,
											hessian=hessian, parallel=parallel_opt,statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											vcrfun=vcrfun,vcrparams=vcrparams,ALL=FALSE)
				stopCluster(parallel_opt$cl)
				
				res			<-	costfunexpgeo(fp=res_optim$par,
											statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
											coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
											vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
			}	else {
				if(sum(h!=statdist[[3]]$h)>0 | sum(k!=statdist[[3]]$k)>0){
					res			<-	costfunexpgeo(fp=c(h,k),statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
												  coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
												  vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
				}	else {
					res			<-	costfunexpgeo(statdist=statdist,p=p,constantr=constantr,ooc_rep=ooc_rep,cs=cs,cofun=cofun,
												  coparams=coparams,crfun=crfun,crparams=crparams,vcofun=vcofun,vcoparams=vcoparams,
												  vcrfun=vcrfun,vcrparams=vcrparams,ALL=TRUE)
				}
			}
		}
	}
	if(shiftfun=="deg"){
		if(length(h)>1 | length(k)>1) {
			if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
			registerDoParallel(parallel_opt$cl)
			Gmtx	<-	foreach(i=h, .combine='rbind') %:%
							foreach(j=k, .combine='rbind') %dopar%{
								c(i,j,costfundeg(fp=c(i,j),
										statdist=statdist,cs=cs,crparams=crparams,cf=crparams,coparams=coparams,p=p,ALL=FALSE)
									)
			}
			stopCluster(parallel_opt$cl)
			
			Gmtx			<-	as.data.frame(Gmtx)
			rownames(Gmtx)	<-	NULL
			colnames(Gmtx)	<-	c("h","k","value")
			res				<-	Gmtx
			class(res)		<-	c("Markov_grid", class(res))
		
		}	else {
			if(OPTIM){
				if(is.null(parallel_opt))	parallel_opt	<-	list(cl=makeCluster(max(c(detectCores()-1,1))),forward=FALSE,loginfo=TRUE)
				res_optim	<-	optimParallel(par=c(h,k), fn=costfundeg, gr=gr, lower=lower, upper=upper, method=method, control=control,
											hessian=hessian, parallel=parallel_opt, statdist=statdist,cs=cs,crparams=crparams,cf=crparams,coparams=coparams,p=p,ALL=FALSE)
				stopCluster(parallel_opt$cl)
				
				res			<-	costfundeg(fp=res_optim$par,statdist=statdist,cs=cs,crparams=crparams,cf=crparams,coparams=coparams,p=p,ALL=TRUE)
			}	else {
				if(sum(h!=statdist[[3]]$h)>0 | sum(k!=statdist[[3]]$k)>0){
					res			<-	costfundeg(fp=c(h,k),statdist=statdist,cs=cs,crparams=crparams,cf=crparams,coparams=coparams,p=p,ALL=TRUE)
				}	else {
					res			<-	costfundeg(statdist=statdist,cs=cs,crparams=crparams,cf=crparams,coparams=coparams,p=p,ALL=TRUE)
				}
			}
		}
	}
	called <- match.call()   
    if(!silent)	print(called)
	return(res)
}


#SIMULATION

Markovsim	<-	function(shiftfun=c("exp","exp-geo"),num=100,h,k,sigma,s,delta,probmix=1,probnbin=0.5,disj=1,RanRep=FALSE,
						 alpha=NULL,beta=NULL,RanSam=FALSE,StateDep=FALSE,a=NULL,b=NULL,q=NULL,z=NULL,detail=100,Vd=50,V,
						 burnin=1)
{
	if(h<=0 | k<0)										stop("Non-applicable h or k parameters were given. h is the time between samplings and k is the critical value. Both should be positive numbers (k is allowed to be 0, but not negative, since positive shifts are assumed).")
	if(sigma<=0 | s<=0 | delta<=0 | V<=0)				stop("Non-applicable sigma, s, delta or V parameters were given. sigma is the process variance, s is the expected number of shifts in a sampling interval, delta is the exponential shift size parameter and V is the maximum distance from the target value taken into account. These should all be positive numbers.")
	if(RanRep & (is.null(alpha) | is.null(beta))) 		stop("Parameters alpha or beta are missing. If RanRep=TRUE, then the beta distribution parameters (alpha and beta) must be given.")
	if(RanSam & StateDep & (is.null(a) | is.null(b)))	stop("Parameters a or b are missing. If RanRep=TRUE and StateDep=TRUE, then the beta distribution parameters (a and b) must be given.")
	if(RanSam & !StateDep & (is.null(q) | is.null(z)))	stop("Parameters q or z are missing. If RanRep=TRUE and StateDep=FALSE, then the logistic function parameters (q and z) must be given.")
	if(Vd<3)											stop("Non-applicable Vd parameter was given. This is a discretisation parameter (number of states after disretisation), which should be an integer value greater than 2.")
	if(disj<=0)											stop("Non-applicable disj parameter was given. disj is the size of a discrete jump in case of exponential-geometric mixture shift distribution and should be a positive number.")
	if(probmix<0 | probmix>1)							stop("Invalid probmix parameter. probmix is the weight of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1..")
	if(probnbin<0 | probnbin>1)							stop("Invalid probnbin parameter. probnbin is the probability parameter of the geometric distribution in case of exponential-geometric mixture shift distribution and should be between 0 and 1.")
	if(num<=0)											stop("Non-applicable num parameter was given. num is the length of the simulation measured in the number of samplings, and should be a positive integer.")
	if(!(shiftfun=="exp" | shiftfun=="exp-geo"))		stop("shiftfun must be either exp or exp-geo.")
	if(detail<2)										stop("Invalid detail parameter was given. detail is the number of data points simulated within a sampling interval and should be an integer greater than one.")
	if(burnin<0 | burnin>=num)							stop("Invalid burnin parameter was given. burnin is the number of samplings deemed as a burn-in period and should be an integer greater than one but less than the number of the total simulated sampling intervals.")

	num		<-	round(num)
	detail	<-	round(detail)
	burnin	<-	round(burnin)

	if(shiftfun=="exp")
	{
		#EXP simulation is creating detail number data points per sampling interval
		eventvec	<-	"start"
		xbase		<-	0
		x			<-	NULL
		for (i in 1:num)
		{
			if(i>1)												xbase		<-	x[length(x)]
			if(i==2)											eventvec	<-	eventvec[2]
			if(eventvec[length(eventvec)]=="alarm" & !RanRep)	xbase		<-	x[length(x)]*0 else	{
				if(eventvec[length(eventvec)]=="alarm")			xbase		<-	x[length(x)]*rbeta(1,alpha,beta) 
			}
			
			stvec		<-	rexp(h*s*detail,s)
			shifttimes	<-	round(cumsum(stvec)[cumsum(stvec)<h],digits=2)
			shiftsizes 	<-	cumsum(rexp(length(shifttimes),1/delta))
			if(length(shifttimes)!=0)
			{
				onetr	<-	NULL
				onetr	<-	c(onetr, rep(xbase,round(shifttimes[1]*detail)))
				for (j in 1:length(shifttimes))
				{
					if(j!=length(shifttimes))	onetr	<-	c(onetr, rep(xbase+shiftsizes[j],round((shifttimes[j+1]-shifttimes[j])*detail)))
					if(j==length(shifttimes))	onetr	<-	c(onetr, rep(xbase+shiftsizes[j],round((h-shifttimes[j])*detail)))
				}
				x	<-	c(x,onetr)
			}	else	{x	<-	c(x,rep(xbase,round(h*detail)))}
			
			if(RanSam) 
			{
				if(StateDep)
				{
					probab	<-	pbeta(x[length(x)]/V,a/h,b)
				}else{
					probab	<-	1/(1+exp(-q*(x[length(x)]-z)))
				}
			}else{
					probab	<-	1
			}
			
			if(rbinom(n=1, size=1, prob=probab)==1)
			{
				eventvec	<-	c(eventvec,"success")
				if(rbinom(n=1, size=1, prob=(1-pnorm(k-x[length(x)],0,sigma)))==1)	eventvec[length(eventvec)]	<-	"alarm"
			} else {
				eventvec	<-	c(eventvec,"failure")
			}
		}
	}
	
	if(shiftfun=="exp-geo")
	{
		#EXP-GEO simulation is creating detail number data points per sampling interval
		eventvec	<-	"start"
		xbase		<-	0
		x			<-	NULL
		for (i in 1:num)
		{
			if(i>1)												xbase		<-	x[length(x)]
			if(i==2)											eventvec	<-	eventvec[2]
			if(eventvec[length(eventvec)]=="alarm" & !RanRep)	xbase		<-	x[length(x)]*0 else	{
				if(eventvec[length(eventvec)]=="alarm")			xbase		<-	x[length(x)]*rbeta(1,alpha,beta) 
			}
			
			stvec		<-	rexp(h*s*detail,s)
			shifttimes	<-	round(cumsum(stvec)[cumsum(stvec)<h],digits=2)
			if(length(shifttimes)!=0)
			{
			pattern		<-	rbinom(length(shifttimes),1,probmix)
			exps		<-	rexp(sum(pattern==0),1/delta)
			geoms		<-	disj*(rgeom(sum(pattern==1),probnbin)+1)
			shifts		<-	rep(0,length(shifttimes))
			shifts[pattern==0]	<-	exps
			shifts[pattern==1]	<-	geoms
			shiftsizes 	<-	cumsum(shifts)
		
				onetr	<-	NULL
				onetr	<-	c(onetr, rep(xbase,round(shifttimes[1]*detail)))
				for (j in 1:length(shifttimes))
				{
					if(j!=length(shifttimes))	onetr	<-	c(onetr, rep(xbase+shiftsizes[j],round((shifttimes[j+1]-shifttimes[j])*detail)))
					if(j==length(shifttimes))	onetr	<-	c(onetr, rep(xbase+shiftsizes[j],round((h-shifttimes[j])*detail)))
				}
				x	<-	c(x,onetr)
			}	else	{x	<-	c(x,rep(xbase,round(h*detail)))}
			
			if(RanSam) 
			{
				if(StateDep)
				{
					probab	<-	pbeta(x[length(x)]/V,a/h,b)
				}else{
					probab	<-	1/(1+exp(-q*(x[length(x)]-z)))
				}
			}else{
					probab	<-	1
			}
		
			if(rbinom(n=1, size=1, prob=probab)==1)
			{
				eventvec	<-	c(eventvec,"success")
				if(rbinom(n=1, size=1, prob=(1-pnorm(k-x[length(x)],0,sigma)))==1)	eventvec[length(eventvec)]	<-	"alarm"
			} else {
				eventvec	<-	c(eventvec,"failure")
			}
		}
	}
	
	
	burntin			<-	x[seq(detail,num*detail,detail)][(burnin+1):length(x[seq(detail,num*detail,detail)])]
	burntevent		<-	eventvec[(burnin+1):length(eventvec)]
	int				<-	seq(0,V,by=(V/(Vd-1)))
	discr_sim_alarm	<-	NULL
	discr_sim_ooc	<-	NULL
	for (i in 1:(length(int)-1))
	{
		discr_sim_alarm	<-	c(discr_sim_alarm, length(burntin[burntin>int[i] & burntin<=int[i+1] & burntevent=="alarm"])/length(burntin))
		discr_sim_ooc	<-	c(discr_sim_ooc,   length(burntin[burntin>int[i] & burntin<=int[i+1] & burntevent!="alarm"])/length(burntin))
	}
	discr_sim_alarm[length(discr_sim_alarm)]	<-	length(burntin[burntin>int[length(int)-1] & burntevent=="alarm"])/length(burntin)
	discr_sim_ooc[length(discr_sim_ooc)]		<-	length(burntin[burntin>int[length(int)-1] & burntevent=="alarm"])/length(burntin)
	discr_sim									<-	c(sum(burntin==0 & burntevent!="alarm")/length(burntin),sum(burntin==0 & burntevent=="alarm")/length(burntin),discr_sim_alarm,discr_sim_ooc)
	
	int2				<-	seq(0,V,by=(V/(Vd-1))) + (V/(Vd-1))*0.5
	names(discr_sim)	<-	c("In-control","False-alarm",paste("Out-of-control",1:(Vd-1),c(round(int2[1:(Vd-2)],3),paste(round(V/(Vd-1)*(Vd-2),3),"+",sep="")),sep="_"),paste("True-alarm",1:(Vd-1),c(round(int2[1:(Vd-2)],3),paste(round(V/(Vd-1)*(Vd-2),3),"+",sep="")),sep="_"))	

	
	res	        <-	list(Value_at_samplings=x[seq(detail,num*detail,detail)], Sampling_event=eventvec, Simulation_data=x, Stationary_distribution=discr_sim)
	class(res)	<-	c("Markov_sim", class(res))
	return(res)
}


utils::globalVariables(c("i","j","h","k","value","shiftfun","s","delta","probmix","probnbin","disj","RanRep","RanSam","StateDep","a","b","z","Vd","V","Qparam","statd","shiftfun","s","delta","probmix","probnbin","disj","RanRep","RanSam","StateDep","a","b","z","Vd","V","Qparam","mtx","statd"))