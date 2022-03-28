rm(list=ls())
library(rfish)
library(tidyverse)

setwd("~/codes/rfish")

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")

library(rfish)

fish = new(Fish)
fish$par$s0 = 0.09637
fish$par$Bhalf_growth = 100e11

pop_K = new(Population, fish)
pop_K$set_superFishSize(1e6)
K_ibm = pop_K$noFishingEquilibriate()


sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop = new(Population, fish)
pop$par$Bhalf = 365426284/4.7
pop$set_superFishSize(2e5)

nsteps = 200
h = 0.3
lf = 45

res_ibm = sim$simulate(pop, lf, h, nsteps, T)
res = simulate(h, lf, F)

d = pop$get_state()
table(d$age)
par(mfrow = c(2,3), mar=c(4,4,1,1))
# plot(nfish~seq(1,nsteps,1), ylab="Number of superfish", xlab="Year")  
# hist(d$age, breaks = seq(0.5,32.5,by=1), main="", xlab="Age")
# plot(y=ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year")
  
ssb.max = max(c(res_ibm$ssb/1e9, res$summaries$SSB/1e9))
plot(y=res_ibm$ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year", col="cyan3", type="l", ylim=c(0,ssb.max))
points(y=res$summaries$SSB/1e9, x=res$summaries$year, type="l")

yield.max = max(c(res_ibm$yield/1e9, res$summaries$Y/1e9))
plot(y=res_ibm$yield/1e9, x=seq(1,nsteps,1), ylab="Yield (MT)", xlab="Year", col="cyan3", type="l", ylim=c(0,yield.max))
points(y=res$summaries$Y/1e9, x=res$summaries$year, type="l")

emp.max = max(c(res_ibm$employment.sea, res$summaries$D.sea))
plot(y=res_ibm$employment.sea, x=seq(1,nsteps,1), ylab="Employment at sea (person-years)", xlab="Year", col="cyan3", type="l", ylim=c(0,emp.max))
points(y=res$summaries$D.sea, x=res$summaries$year, type="l")

emp.max = max(c(res_ibm$employment.shore, res$summaries$D.shr))
plot(y=res_ibm$employment.shore, x=seq(1,nsteps,1), ylab="Employment on shore (person-years)", xlab="Year", col="cyan3", type="l", ylim=c(0,emp.max))
points(y=res$summaries$D.shr, x=res$summaries$year, type="l")

p.sea.max = max(c(res_ibm$profit.sea/1e9, res$summaries$P.sea/1e9))
p.sea.min = min(c(res_ibm$profit.sea/1e9, res$summaries$P.sea/1e9))
plot(y=res_ibm$profit.sea/1e9, x=seq(1,nsteps,1), ylab="Profit at sea (Billions NOK)", xlab="Year", col="cyan3", type="l", ylim=c(p.sea.min,p.sea.max))
points(y=res$summaries$P.sea/1e9, x=res$summaries$year, type="l")

p.shr.max = max(c(res_ibm$profit.shore/1e9, res$summaries$P.shr/1e9))
p.shr.min = min(c(res_ibm$profit.shore/1e9, res$summaries$P.shr/1e9))
plot(y=res_ibm$profit.shore/1e9, x=seq(1,nsteps,1), ylab="Profit at shore (Billions NOK)", xlab="Year", col="cyan3", type="l", ylim=c(p.shr.min,p.shr.max))
points(y=res$summaries$P.shr/1e9, x=res$summaries$year, type="l")

par(mfrow=c(1,1))
d = pop$get_state()
d1 = d %>% group_by(age) %>% summarize(mat = length(which(isMature))/length(isMature))
plot(mat, type="l")
points(d1$mat~d1$age, type="o", col="cyan3")
# d = get_stocks()
# pnow = maturation_probability(1:30, 10, maturation_steepness(d$pmrn_width, d$pmrn_envelope), d$pmrn_slope, d$pmrn_intercept)
# pnext = maturation_probability(2:31, 10, maturation_steepness(d$pmrn_width, d$pmrn_envelope), d$pmrn_slope, d$pmrn_intercept)
# plot(x=1:30, y=(pnext-pnow)/(1-pnow))
# plot(x=1:30, y=pnow)
# 
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# a = -3.5
# b = 0.363
# c = -0.97
# d = 0.0059
# e = -0.44
# curve(5/(1+exp(-(a+b*log(1)+c*log(500)+d*500+e*log(x)))), from=0.1, to=10, xlab="Growth")
# curve(5/(1+exp(-(a+b*log(20)+c*log(x)+d*x+e*log(.1)))), from=10, to=600, xlab="Diameter")
# curve(5/(1+exp(-(a+b*log(x)+c*log(500)+d*500+e*log(5)))), from=0.01, to=1, xlab="Light")
# 

#### Simulate vs hvec #####
library(rfish)

fish.par = new(FishParams)
fish.par$s0 = 0.09637
fish.par$Bhalf_growth = 100e11

fish = new(Fish)
fish.par$flag = 1
fish.par$use_old_model_mat = T
fish$par = fish.par
pop_K = new(Population, fish)
pop_K$set_superFishSize(1e6)
K_ibm = pop_K$noFishingEquilibriate()

h = 0.3
lf = 45

nsteps = 200
sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop.par = new(PopulationParams)
pop.par$Bhalf = 365426284/4.7

pop = new(Population, fish)
pop$par = pop.par
pop$set_superFishSize(2e5)
hvec = seq(0, 0.99, length.out = 20)
dat_ibm = data.frame(matrix(ncol=7, nrow=0))
colnames(dat_ibm) =c("ssb", "yield", "employment.sea", "employment.shore", "profit.sea", "profit.shore", "tsb")
dat = data.frame(matrix(ncol=7, nrow=0))
colnames(dat) =c("ssb", "yield", "employment.sea", "employment.shore", "profit.sea", "profit.shore", "tsb")
for (ih in 1:length(hvec)){
  res_ibm = sim$simulate(pop, lf, hvec[ih], nsteps, T)
  res = simulate(hvec[ih], lf, F)
  v_ibm = colMeans(res_ibm[151:200,])/c(1e9, 1e9, 1, 1, 1e9, 1e9)
  v = colMeans((res$summaries %>% select(SSB, Y, D.sea, D.shr, P.sea, P.shr, TSB))[151:200,])/c(1e9, 1e9, 1, 1, 1e9, 1e9)
  dat_ibm[nrow(dat_ibm)+1,] = v_ibm
  dat[nrow(dat)+1,] = v
}

par(mfrow = c(2,3), mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.5)
matplot(ylab="ssb (MT)",y=cbind(dat_ibm$ssb, dat$ssb, dat_ibm$tsb, dat$tsb), x=hvec, type="l", col=c("cyan2", "black", "orange", "magenta"), lwd=2, lty=1, xlab="Harvest prop")
matplot(ylab="yield (MT)",y=cbind(dat_ibm$yield, dat$yield), x=hvec, type="l", col=c("cyan2", "black"), lwd=2, lty=1, xlab="Harvest prop")
matplot(ylab="emp at sea (PY)",y=cbind(dat_ibm$employment.sea, dat$employment.sea), x=hvec, type="l", col=c("cyan2", "black"), lwd=2, lty=1, xlab="Harvest prop")
matplot(ylab="emp on shore (PY)",y=cbind(dat_ibm$employment.shore, dat$employment.shore), x=hvec, type="l", col=c("cyan2", "black"), lwd=2, lty=1, xlab="Harvest prop")
matplot(ylab="profit at sea (Bn NOK)",y=cbind(dat_ibm$profit.sea, dat$profit.sea), x=hvec, type="l", col=c("cyan2", "black"), lwd=2, lty=1, xlab="Harvest prop")
matplot(ylab="profit on shore (Bn NOK)",y=cbind(dat_ibm$profit.shore, dat$profit.shore), x=hvec, type="l", col=c("cyan2", "black"), lwd=2, lty=1, xlab="Harvest prop")

