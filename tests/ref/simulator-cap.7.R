### initialisation

# initialisation of age-specific vectors
L<-Linf*(1-exp(-k*(1:ages-a0)))		# length-at-age (cm)
W<-theta*L^g/1000					# weight-at-age (kg)

# initialisation of other data structures
N<-matrix(0,nrow=years+1,ncol=ages)	# numbers at age
N[1,1]<-N.ini
# for(i in 2:ages) N[1,i]<-N.ini*exp(-sum(12*M[1:(i-1)])) 	# numbers at age, year 1
Res<-data.frame(year=1:years,N=NA,TSB=NA,SSB=NA,Y=NA,R=NA,E=NA,Ereal=NA,Fratio=NA,F.hypo=NA,F.hypo.SSB=NA,status=NA,
	D.sea=NA,D.shr=NA,D.tot=NA,P.sea=NA,P.shr=NA,P.tot=NA) # results

# find out carrying capacity (Z=M)
# independent of fishing parameters, so this needs to be done only once
N.<-matrix(0,nrow=10051,ncol=ages)	# this calculation should not be influenced by the length of the main simulation
N.[1,]<-N[1,]
for(year in 1:(nrow(N.)-1)){				
	SSB<-sum(N.[year,]*W*mat*exp(-M))
	N.[year+1,1]<-r0*SSB/(1+SSB/Bhalf)*exp(rnorm(n=1,mean=-sf^2/2,sd=sf))
	N.[year+1,2:ages]<-N.[year,1:(ages-1)]*exp(-12*Mi[1:(ages-1)])*(1-mat[1:(ages-1)])
}
K<-apply(N.[51:10050,],2,mean) # skip first 50 years as the transient
rm(N.)

### functions
# image that uses gray scale as the default
Image<-function(...){image(...,col=gray((32:0)/32)); box()}

# function needed to find fishing mortality corresponding to a certain effort
simul<-function(F.max,parms){ # parms=N.rel,E.real
	E<-parms[1]^(1-b)*F.max*(exp(-(F.max+M[ages])*(1-b))-1)/(q*(F.max+M[ages])*(b-1))
	(E-parms[2])^2
}

# Baranov catch equation - solve F when M, Y and B0 are known
Bce<-function(F,parms){
	Y<-parms[1,1]
	N0<-parms[,2]
	Y.<-sum(F/(F+M)*(1-exp(-F-M))*N0*W)
	sum((Y-Y.)^2)
}

# main simulation function
simulate<-function(h.prop,LF50,display=FALSE){
	F.max<- -log(1-h.prop)
	age.fish<-min((1:ages)[1/(1+exp(-SF*(L-LF50)))>0.5])
	F2E.factor<-ifelse(b==1,F.max/q,F.max*(exp(-(F.max+M[ages])*(1-b))-1)/(q*(F.max+M[ages])*(b-1))) # Liu & Heino 2014
	sel.fact<-1/(1+exp(-SF*(L-LF50)))
	for(year in 1:years){ # 
		if(h.prop>0){ # & sum(N[year,]*W*mat*exp(-M))>Blim){	# fishing is possible
			# determine the effective fishing mortality
			N.rel<-mean(N[year,age.fish:ages]/K[age.fish:ages])
			if(sum(N[year,]*W*mat*exp(-M-F.max*sel.fact))>Blim){
				F.max.<-F.max	# normal F
				E<-N.rel^(1-b)*F2E.factor			# effort corresponding to F (F2E.factor contains only terms that are year-independent)
				status<-2 # nominal F attempted
			}else{				# reduced F when full F would bring SSB<200kt
				F.max.<- max(0,-mean(M[3:6]) -log(Blim/sum(N[year,3:6]*W[3:6]*mat[3:6])))
				E<-ifelse(F.max.==0,0,N.rel^(1-b)*ifelse(b==1,F.max./q,F.max.*(exp(-(F.max.+M[ages])*(1-b))-1)/(q*(F.max.+M[ages])*(b-1))))
				status<-1 # reduced F attempted
			}
			D.sea.demand<-d.sea*E				# corresponding employment demand
			D.eff<-d.max*D.sea.demand/(d.max+D.sea.demand)	# effective man hours
			E.real<-D.eff/d.sea					# realized effort (~vessel days)
			F.real<-ifelse(F.max.==0,0,ifelse(b==1,q*E.real,optimize(f=simul,interval=c(0,F.max.),parms=c(N.rel,E.real))$minimum))
      F.real=F.max # Jaideep
      
			# population dynamics
			F<-F.real*sel.fact					# fishing mortality at age
			Z<-F+M								# total mortality for mature
			Y<-sum(F/Z*(1-exp(-Z))*N[year,]*W*mat)	# mature catch
			# calculate the hypothetical F corresponding to unselective harvesting of all age groups
			F.hypo<-optimize(f=Bce,interval=c(0,F.max),parms=data.frame(Y,N[year,]))$minimum
			F.hypo.SSB<-optimize(f=Bce,interval=c(0,F.max),parms=data.frame(Y,mat*N[year,]))$minimum
			TSB<-sum(N[year,]*W*exp(-Z))
			SSB<-sum(N[year,]*W*mat*exp(-Z))
		}else{
			Y<-0; D.sea.demand<-0; E<-0; E.real<-0; F.real<-0; F.hypo<-0; F.hypo.SSB<-0; status<-0
			TSB<-sum(N[year,]*W*exp(-M))
			SSB<-sum(N[year,]*W*mat*exp(-M))
		}
		# N[year+1,1]<-min(r0*SSB/(1+SSB/Bhalf)*exp(rnorm(n=1,mean=-sf^2/2,sd=sf)),max.R)
		N[year+1,1]<-min(r0*SSB/(1+SSB/Bhalf),max.R)
		N[year+1,2:ages]<-N[year,1:(ages-1)]*exp(-12*Mi[1:(ages-1)])*(1-mat[1:(ages-1)])

		# economic calculations
		profit.sea<-Y*price.sea - scale.catch.cap*(D.sea.demand*s.sea + E*varC.sea) - scale.catch.all*fixedC.sea
		profit.shr<-Y*(price.shr-price.sea) - Y*d.shr*s.shr - scale.catch.all*fixedC.shr
		Res[year,]<-c(year,sum(N[year,]),TSB,SSB,Y,N[year+1,1],E,E.real,F.real/F.max,F.hypo,F.hypo.SSB,status,
			D.sea.demand,Y*d.shr,D.sea.demand+Y*d.shr,profit.sea,profit.shr,profit.sea+profit.shr)
	}
	if(display==TRUE){
		op<-par(mfrow=c(2,2),pty="s",mar=c(4,4,.5,.5))
		plot(Res$R/1e9,type="l",xlab="",ylab="Recruitment (billions)")
		plot(Res$SSB/1e6,type="l",xlab="",ylab="SSB (1000t)")
		plot(Res$Y/1e6,type="l",xlab="",ylab="Total catch weight (1000t)")
		Image(x=1:years,y=1:ages,z=log(N[1:years,]),xlab="",ylab=""); legend("top","Age structure",bty="n") # non-normalized
		par(op)
	}
	list(N=N[1:years,],summaries=Res)
}

hvec = seq(0,0.8,0.1)
ovec = numeric(length(hvec))
for(i in 1:length(hvec)){
  res = simulate(hvec[i], LF50, F)
  ovec[i] = mean(res$summaries$SSB[80:100])
# plot(y=res$summaries$SSB/1e9, x=res$summaries$year)
}
plot(ovec~hvec)
