rm(list=ls())
library(rfish)
library(tidyverse)

fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
#pop = new(Population, fish)

# sim = new(Simulator)
# res_ibm = sim$simulate(pop, 45, 0, 200, T)

pop = new(Population, fish)
# h = 0
# pop$mort_fishing_mature = -log(1-h)
# pop$mort_fishing_immature = 0
pop$set_harvestProp(0)
pop$set_superFishSize(2e3)
pop$init(1000)
pop$get_state()

nsteps = 200
nfish = numeric(nsteps)
nfish[1]=1000
ssb = numeric(nsteps)
dat = data.frame(ssb=0, yield=0, emp_sea=0, emp_shore=0, profit_sea=0, profit_shr=0)
for (i in 2:nsteps){
  v = pop$update()
  dat = rbind(dat, v)
  nfish[i] = pop$nfish()
}
dat = dat[-1,]

d = pop$get_state()
dist = table(d$age, d$length)

#K_ibm = pop$calcK()

par(mfrow = c(2,2), mar=c(4,4,1,1))
plot(nfish~seq(1,nsteps,1), ylab="Number of superfish", xlab="Year")
image(x=as.numeric(rownames(dist)), y = as.numeric(colnames(dist)), z=log(1+3*log(dist)), col=scales::viridis_pal()(100), xlab="", ylab="")
plot(dat$ssb/1e9, ylab="SSB (MT)", xlab="Year")
#plot(K_ibm[-1], log="xy", xlab="Age", ylab="Frequency")

K = matrix(nrow=10, ncol=10)
hvec = seq(0,1,length.out=10)
l50vec = seq(20,160,length.out=10)
for (ih in 1:length(hvec)){
  for (il50 in 1:length(l50vec)){
    pop$set_harvestProp(hvec[ih])
    pop$set_minSizeLimit(l50vec[il50])
    K[ih,il50] = pop$fishableBiomass()
  }
}

Bfishable_nonfished = K[1,1:10]/1e9
pop_noFish = pop

pop = new(Population, fish)
# h = 0
# pop$mort_fishing_mature = -log(1-h)
# pop$mort_fishing_immature = 0
pop$set_harvestProp(0.3)
pop$set_minSizeLimit(45)
pop_noFish$set_harvestProp(0.3)
pop_noFish$set_minSizeLimit(45)
pop$K = pop_noFish$fishableBiomass()
pop$set_superFishSize(5e2)
pop$init(1000)
pop$get_state()

nsteps = 200
nfish = numeric(nsteps)
nfish[1]=1000
ssb = numeric(nsteps)
dat = data.frame(ssb=0, yield=0, emp_sea=0, emp_shore=0, profit_sea=0, profit_shr=0)
for (i in 2:nsteps){
  v = pop$update()
  dat = rbind(dat, v)
  nfish[i] = pop$nfish()
}
dat = dat[-1,]

d = pop$get_state()
dist = table(d$age, d$length)

#K_ibm = pop$calcK()

par(mfrow = c(2,3), mar=c(4,4,1,1))
plot(nfish~seq(1,nsteps,1), ylab="Number of superfish", xlab="Year")
image(x=as.numeric(rownames(dist)), y = as.numeric(colnames(dist)), z=log(1+3*log(dist)), col=scales::viridis_pal()(100), xlab="", ylab="")
plot(dat$ssb/1e9, ylab="SSB (MT)", xlab="Year")
plot(dat$yield/1e9, ylab="Yield (MT)", xlab="Year", col="cyan3", type="l")
matplot(cbind(dat$emp_sea, dat$emp_shore), ylab="Employment at sea/shore (Person years)", xlab="Year", col=c("seagreen", "bisque3"), type="l", lty=1, lwd=2)
matplot(cbind(dat$profit_sea, dat$profit_shr)/1e9, ylab="Profit at sea/shore (Billion NOK)", xlab="Year", col=c("seagreen", "bisque3"), type="l", lty=1, lwd=2)

#plot(K_ibm[-1], log="xy", xlab="Age", ylab="Frequency")

K = matrix(nrow=10, ncol=10)
hvec = seq(0,1,length.out=10)
l50vec = seq(20,160,length.out=10)
for (ih in 1:length(hvec)){
  for (il50 in 1:length(l50vec)){
    pop$set_harvestProp(hvec[ih])
    pop$set_minSizeLimit(l50vec[il50])
    K[ih,il50] = pop$fishableBiomass()
  }
}

Bfishable_h0.1 = K[1,1:10]/1e9

# matplot(x = l50vec, y= cbind(Bfishable_nonfished, Bfishable_h0.1), col=c("black", "blue"), type="l", lty=1, ylab="Fishable Biomass (MT)", xlab="L50")
