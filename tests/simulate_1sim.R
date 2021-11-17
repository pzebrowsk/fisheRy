rm(list=ls())
library(rfish)
library(tidyverse)
a = 1:30
l = numeric(30)
f = new(Fish)

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")


library(rfish)
fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
h = 0.5
pop = new(Population, fish)

pop$set_harvestProp(h)

Kibm = pop$calcK()

pop$init(1000)
pop$get_state()
pop$print_summary()

nsteps = 200
nfish = numeric(nsteps)
nfish[1]=1000
ssb = numeric(nsteps)
for (i in 2:nsteps){
  ssb[i] = pop$update()
  nfish[i] = pop$nfish()
}

ssb
d = pop$get_state()
table(d$age)
par(mfrow = c(2,2), mar=c(4,4,1,1))
plot(nfish~seq(1,nsteps,1), ylab="Number of superfish", xlab="Year")
hist(d$age, breaks = seq(0.5,32.5,by=1), main="", xlab="Age")
# plot(y=ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year")

res = simulate(h, LF50, F)
ssb.max = max(c(ssb/1e9, res$summaries$SSB/1e9))
plot(y=ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year", col="red", ylim=c(0,ssb.max))
points(y=res$summaries$SSB/1e9, x=res$summaries$year, type="l")
