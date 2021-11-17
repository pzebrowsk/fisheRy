#--------------------------------------------------------------------------------------------------------------------------
m <- matrix(data=NA, nrow=5, ncol=4)  ## 5 SHs rows, 4 UCs columns
m[1,] <- c(0.3, 0.0, 0.7, 0.0)  ## ind              ## original stakeholder preference (SP) vectors
m[2,] <- c(0.5, 0.1, 0.1, 0.3)  ## art
m[3,] <- c(0.2, 0.5, 0.0, 0.3)  ## employ
m[4,] <- c(0.2, 0.0, 0.6, 0.2)  ## prof
m[5,] <- c(0.1, 0.2, 0.2, 0.5)  ## con

#--------------------------------------------------------------------------------------------------------------------------
gm<-function(data,alpha=1,na.rm=FALSE){ # generalized mean, with some violence done to input data unless alpha=1
	if(na.rm==TRUE) data<-data[!is.na(data)]
	#if(alpha!=1) data<-pmax(0,data) # force values to be non-negative for alpha!=1 - prevents errors, but the results may be useless
	if(alpha==0){prod(data)^(1/length(data))}
	else{(sum(data^alpha)/length(data))^(1/alpha)}
}

# backward-looking exponential avarage for a finite time window
EAB<-function(data,T50,period=NA,slide=FALSE){
# T50 = half time, period = length of the observation window (if NA, then all)
# if slide=TRUE, then returns vector
	if(T50<0.1) stop("Too short halving time\n")
	if(!is.na(period)) data<-data[(length(data)-period):length(data)]
	n<-length(data)
	i<-1:n
	if(slide){
		x<-sum(data*2^((i-n)/T50))/sum(2^((i-n)/T50))
		while(length(data)>1){
			data<-data[1:(n-1)]
			n<-length(data)
			i<-1:n
			x<-c(sum(data*2^((i-n)/T50))/sum(2^((i-n)/T50)),x)
		}
		return(x)
	}
	else sum(data*2^((i-n)/T50))/sum(2^((i-n)/T50))
}

# forward-looking exponential avarage for a finite time window
EAF<-function(data,T50,period=NA,slide=FALSE){
# T50 = half time, period = length of the observation window (if NA, then all)
# if slide=TRUE, then returns vector
	if(T50<0.1) stop("Too short halving time\n")
	if(!is.na(period)) data<-data[1:period]
	n<-length(data)
	i<-1:n
	if(slide){
		x<-sum(data*2^((n-i)/T50))/sum(2^((n-i)/T50))
		while(length(data)>1){
			data<-data[2:n]
			n<-length(data)
			i<-1:n
			x<-c(x,sum(data*2^((n-i)/T50))/sum(2^((n-i)/T50)))
		}
		return(x)
	}
	else sum(data*2^((n-i)/T50))/sum(2^((n-i)/T50))
}

#--------------------------------------------------------------------------------------------------------------------------
normalize_max<-function(a,debug=FALSE,SSB=NA){
# this normalizes utility components so that the maximal case-specific mean is 1 (negative values stay negative)
# if debug is TRUE, intermediate results are kept (raw utility components and their maxima)
# if SSB is given, then the normalizing takes place with respect to that value (necessary when non-profitable options have been removed, and data do not contain non-fished populations)

# 1) calculate first replicate- and management option-specific means
	means<-aggregate(a[,c("y","e","p","s")], list(rep=a$it,ms=a$ms,hl=a$hl),mean,na.rm=T)
# 2) then find maximal means for each replicate
	maxs<-aggregate(means[,c("y","e","p","s")],list(rep=means$rep),max,na.rm=T)
	if(min(maxs)<=0) stop('Non-positive maximum - normalization not possible!\n')
# 3) normalize
	tmp<-merge(a,maxs,by.x="it",by.y="rep")
	# because of R's way of dealing with non-unique names, now original values have suffix .x, maxs suffix .y
	tmp$y<-tmp$y.x/tmp$y.y
	tmp$e<-tmp$e.x/tmp$e.y
	tmp$p<-tmp$p.x/tmp$p.y
	if(is.na(SSB)) tmp$s<-tmp$s.x/tmp$s.y else tmp$s<-tmp$s.x/SSB
	if(!debug) tmp<-tmp[,c("it","yr","ms","hl","y","e","p","s")]
	if(debug) print(apply(tapply(tmp$y,list(tmp$it,tmp$hl,tmp$ms),mean,na.rm=T),1,max,na.rm=T))
	tmp[order(tmp$it,tmp$ms,tmp$hl,tmp$yr),] #### NOTE NB!!!! need to take out tmp$yr and remember to aggregate over years if you want maxes to =1
}

#--------------------------------------------------------------------------------------------------------------------------
SHutility<-function(d,m,first.year=NA,T50=NA){ 
# calculate stakeholder utilities

# drop transient years, if desired
	if(!is.na(first.year)) d<-d[d$yr>=first.year,]

	# calculate SH utilities
	d$ind <- (m[1,1] * d$y) +  (m[1,2] * d$e) +  (m[1,3] * d$p) +  (m[1,4] * d$s)
	d$art <- (m[2,1] * d$y) +  (m[2,2] * d$e) +  (m[2,3] * d$p) +  (m[2,4] * d$s)
	d$emp <- (m[3,1] * d$y) +  (m[3,2] * d$e) +  (m[3,3] * d$p) +  (m[3,4] * d$s)
	d$pro <- (m[4,1] * d$y) +  (m[4,2] * d$e) +  (m[4,3] * d$p) +  (m[4,4] * d$s)
	d$con <- (m[5,1] * d$y) +  (m[5,2] * d$e) +  (m[5,3] * d$p) +  (m[5,4] * d$s)
	d
}

#--------------------------------------------------------------------------------------------------------------------------
SHsatisfaction.EA.exp<-function(d,T50=Inf,T50.discounting=Inf,q=0){
	# "shifting baselines" calculation if T50 is specified,  Equation 3d
	# for time scales of "Inf", these collapse to the simple case without any discounting
	if(T50!=Inf){
		d<-d[order(d$it,d$ms,d$hl,d$yr),]
		years<-(max(d$yr)-min(d$yr)+1)
		tmp<-as.vector(as.matrix(d[,c("ind","art","emp","pro","con")]))
		tmp<-matrix(tmp,byrow=FALSE,nrow=years)
		eab.inf<-apply(tmp,2,EAB,T50=1e99,slide=TRUE)
		eab.T50<-apply(tmp,2,EAB,T50=T50,slide=TRUE)
		tmp<-tmp*2^(q*(eab.inf-eab.T50))
		tmp<-matrix(tmp,byrow=FALSE,nrow=nrow(d))
		d[,c("ind","art","emp","pro","con")]<-tmp
	}
	# find the control-option specific time-averaged maxima
	d$id<-paste(d$it,d$hl,d$ms) # unique "case" identity
	# case-specific means - direct mean or discounted from the present
	if(T50.discounting==Inf) means<-aggregate(d[,c("ind","art","emp","pro","con")],list(id=d$id),mean,na.rm=TRUE) # direct mean
	else means<-aggregate(d[,c("ind","art","emp","pro","con")],list(id=d$id),EAF,T50=T50.discounting) # discounted from the present
	means$it<-aggregate(d$it,list(id=d$id),mean)$x # replicate numbers
	maxs<-aggregate(means[,c("ind","art","emp","pro","con")],list(it=means$it),max, na.rm=TRUE) # identify max within a replicate
	d<-merge(d,maxs,by="it")
	
	# normalization with respect to the maximum - these are year-specific SH satisfactionss
	d$ind <- d$ind.x/d$ind.y
	d$art <- d$art.x/d$art.y
	d$emp <- d$emp.x/d$emp.y
	d$pro <- d$pro.x/d$pro.y
	d$con <- d$con.x/d$con.y
	# case-specific time-averaged satisfaction - direct mean or discounted from the present
	if(T50.discounting==Inf) aggregate(d[,c("ind","art","emp","pro","con")],list(it=d$it,ms=d$ms,hl=d$hl),mean,na.rm=T) # direct mean
	else aggregate(d[,c("ind","art","emp","pro","con")],list(it=d$it,ms=d$ms,hl=d$hl),EAF,T50=T50.discounting) # discounted from the present
}

#--------------------------------------------------------------------------------------------------------------------------
JSS<-function(d,all=FALSE){
# calculate JSS
	d$SHUmin <- apply(d[,c("ind","art","emp","pro","con")],1,min)
	if(all==TRUE){
		d$SHUn4 <- apply(d[,c("ind","art","emp","pro","con")],1,gm,alpha=-4)
		d$SHUhm <- apply(d[,c("ind","art","emp","pro","con")],1,gm,alpha=-1)
		d$SHUmean <- apply(d[,c("ind","art","emp","pro","con")],1,gm,alpha=1)
		}
	d
}

#----------- F U N C T I O N S   F O R   F I G  4 C,D ----------------------

# function that calculates realized JSS for replicates and shorter time windows
# realized = satisfaction relative to what is possible within the whole time period

split.JSS<-function(data,period=10){ # period = time window for satisfaction assessment - integer fraction of 50

	u.std<-normalize_max(data)

	# remove cases with negative longterm profit - on average, over replicates
	u.std$id<-paste(u.std$hl,u.std$ms) # unique "case" identity
	mean.p<-aggregate(u.std$p, list(id=u.std$id),mean,na.rm=T)
	u.std<-merge(u.std,mean.p,by="id")
	u.std<-u.std[u.std$x>0,!names(u.std) %in% c("id","x")]

	su<-SHutility(u.std,m) #equation 2b
	su$key<-paste(su$it,su$ms,su$hl)

	ss<-SHsatisfaction.EA.exp(su) #equation 3d
	jss<-JSS(ss) #equation 4

	# identify maximum jss for each replicate
	jss<-aggregate(jss$SHUmin,list(it=jss$it,ms=jss$ms,hl=jss$hl),mean) # JSS for each ctrl option and replicate
	winners<-aggregate(jss$x,list(it=jss$it),max)
	jss<-merge(jss,winners,by="it")
	jss<-jss[jss$x.x==jss$x.y,]
	jss$key<-paste(jss$it,jss$ms,jss$hl)

	# mean utilities for the cases that correspond to max JSS
	ss.bis<-su[su$key %in% jss$key,]
	ss.bis<-aggregate(ss.bis[,c("ind","art","emp","pro","con")],list(it=ss.bis$it,ms=ss.bis$ms,hl=ss.bis$hl,period=(ss.bis$yr-1)%/%period),mean,na.rm=TRUE)

	# scale relative to max JSS for the whole 50-yr period
	ss.bis<-merge(ss.bis,jss[,c("x.x","it")],by="it")
	ss.bis[,c("ind","art","emp","pro","con")]<-sweep(ss.bis[,c("ind","art","emp","pro","con")],1,ss.bis$x.x,"/")

	# calculate JSS
	jss<-JSS(ss.bis)
	jss<-jss[order(jss$it,jss$period),]
	return(jss$SHUmin)
}


#--------------------------------------------------------------------------------------------------------------------------
legend.col <- function(col, lev){
## continous legend function from: http://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/

	opar <- par
	 
	n <- length(col)
	 
	bx <- par("usr")
	 
	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n
	 
	xx <- rep(box.cx, each = 2)
	 
	par(xpd = TRUE)
	for(i in 1:n){
	 
	yy <- c(box.cy[1] + (box.sy * (i - 1)),
	box.cy[1] + (box.sy * (i)),
	box.cy[1] + (box.sy * (i)),
	box.cy[1] + (box.sy * (i - 1)))
	polygon(xx, yy, col = col[i], border = col[i])
	 
	}
	par(new = TRUE)
	plot(0, 0, type = "n",
	ylim = c(min(lev), max(lev)),
	yaxt = "n", ylab = "",
	xaxt = "n", xlab = "",
	frame.plot = FALSE)
	axis(side = 4, las = 2, tick = FALSE, line = .25)
	par <- opar
}  