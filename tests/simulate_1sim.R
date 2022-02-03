rm(list=ls())
library(rfish)
library(tidyverse)

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")

library(rfish)
fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
pop_K = new(Population, fish)
pop_K$set_superFishSize(500)
K_ibm = pop_K$noFishingEquilibriate()

h = 0.3
lf = 45

nsteps = 200
sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop = new(Population, fish)
pop$set_superFishSize(500)
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
