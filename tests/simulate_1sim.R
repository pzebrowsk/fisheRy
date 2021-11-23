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
pop = new(Population, fish)
K_ibm = pop$calcK()
h = 0.5

nsteps = 200
sim = new(Simulator)
res_ibm = sim$simulate(pop, h, nsteps, T)
res = simulate(h, LF50, F)

d = pop$get_state()
table(d$age)
par(mfrow = c(2,3), mar=c(4,4,1,1))
# plot(nfish~seq(1,nsteps,1), ylab="Number of superfish", xlab="Year")  
hist(d$age, breaks = seq(0.5,32.5,by=1), main="", xlab="Age")
# plot(y=ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year")
  
ssb.max = max(c(res_ibm$ssb/1e9, res$summaries$SSB/1e9))
plot(y=res_ibm$ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year", col="red", ylim=c(0,ssb.max))
points(y=res$summaries$SSB/1e9, x=res$summaries$year, type="l")

yield.max = max(c(res_ibm$yield/1e9, res$summaries$Y/1e9))
plot(y=res_ibm$yield/1e9, x=seq(1,nsteps,1), ylab="Yield (MT)", xlab="Year", col="red", ylim=c(0,yield.max))
points(y=res$summaries$Y/1e9, x=res$summaries$year, type="l")

emp.max = max(c(res_ibm$employment, res$summaries$Ereal))
plot(y=res_ibm$employment, x=seq(1,nsteps,1), ylab="Employment (person-years)", xlab="Year", col="red", ylim=c(0,emp.max))
points(y=res$summaries$Ereal, x=res$summaries$year, type="l")

p.sea.max = max(c(res_ibm$profit.sea/1e9, res$summaries$P.sea/1e9))
p.sea.min = min(c(res_ibm$profit.sea/1e9, res$summaries$P.sea/1e9))
plot(y=res_ibm$profit.sea/1e9, x=seq(1,nsteps,1), ylab="Profit at sea (Billions NOK)", xlab="Year", col="red", ylim=c(p.sea.min,p.sea.max))
points(y=res$summaries$P.sea/1e9, x=res$summaries$year, type="l")

p.shr.max = max(c(res_ibm$profit.shore/1e9, res$summaries$P.shr/1e9))
p.shr.min = min(c(res_ibm$profit.shore/1e9, res$summaries$P.shr/1e9))
plot(y=res_ibm$profit.shore/1e9, x=seq(1,nsteps,1), ylab="Profit at shore (Billions NOK)", xlab="Year", col="red", ylim=c(p.shr.min,p.shr.max))
points(y=res$summaries$P.shr/1e9, x=res$summaries$year, type="l")

