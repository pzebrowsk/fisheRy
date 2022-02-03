rm(list=ls())
library(rfish)
library(tidyverse)

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")

#### MAIN ####

fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
pop_K = new(Population, fish)
pop_K$set_superFishSize(500)
K_ibm = pop_K$noFishingEquilibriate()

hist(pop_K$get_state()$age, ylab="Frequency", xlab="Age", breaks=21)
hist(pop_K$get_state()$length, ylab="Frequency", xlab="length", breaks=21)

# plot(K/1e6, log="y", type="l")
# points(K_ibm[-1])

sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop = new(Population, fish)
pop$set_superFishSize(500)

res_ibm_sq = sim$simulate(pop, 45, 0.41, 200, T)
res_sq = simulate(0.41, 45, F)

# ### Continuation from Status quo (test)
# res_ib_csq = sim$simulate(pop, 45, 0.2, 200, F)
# with(rbind(res_ibm_sq, res_ib_csq), plot(ssb/1e9, type="l", ylab="SSB (MT)", xlab="Yr"))
# ###

pop_ref = pop

# hvec = seq(0.01,0.45,0.02)
hvec = seq(0.01,0.8,length.out = 20)
lfvec = seq(20,160,length.out = 20)


niter = 10
jss_arr = array(dim=c(3, length(hvec), length(lfvec), niter))
full_arr = array(dim=c(c(50, length(hvec), length(lfvec), 4, niter)))
for (iter in 1:niter){
  
  pop = pop_ref
  
  res_ibm_full = sim$simulate_multi_2d(pop, lfvec, hvec, 50, F)
  arr = array(data=res_ibm_full, dim=c(50, length(hvec), length(lfvec), 4))
  full_arr[,,,,iter] = arr
  
  #sim = new(Simulator)
  d = sim$max_avg_utils_2d(c(4,length(lfvec),length(hvec),50), res_ibm_full)
  utils = array(data=d, dim=c(length(hvec),length(lfvec), 4))
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  for (i in 1:4){
    image(x=hvec, y=lfvec, z=utils[,,i], col=scales::viridis_pal()(100), zlim=c(0,1), main=c("ssb","yield","employment","profit")[i])
  }
  
  #sim = new(Simulator)
  ss = sim$stakeholder_satisfaction_2d(c(4,length(lfvec),length(hvec),50), res_ibm_full)
  scs = array(data=ss, dim=c(length(hvec),length(lfvec), 5))
  par(mfrow=c(2,3), mar=c(4,4,4,1))
  for (i in 1:5){
    image(x=hvec, y=lfvec, z=scs[,,i], col=scales::viridis_pal()(100), zlim=c(0,1), main=c("Industrial", "Artisanal", "Emp. policymakers", "Profit. policymakers", "Conservationists")[i])
  }
  
  JSS_min = apply(scs, c(1,2), min)
  JSS_mean = apply(scs, c(1,2), mean)
  JSS_hmean = 1/apply(1/scs, c(1,2), mean)
  
  jss_arr[1,,,iter] = JSS_mean
  jss_arr[2,,,iter] = JSS_hmean
  jss_arr[3,,,iter] = JSS_min
  
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  for (i in 1:3){
    image(x=hvec, y=lfvec, z=jss_arr[i,,,iter], col=rainbow(start = 0, end = 0.75, n=100), zlim=c(0,1), main=c("JSS mean","JSS hmean","JSS min")[i])
  }
    
}

JSS_ensavg = apply(jss_arr, MARGIN = c(1,2,3), FUN = mean)
par(mfrow=c(2,2), mar=c(4,4,4,1))
for (i in 1:3){
  image(x=hvec, y=lfvec, z=JSS_ensavg[i,,], col=scales::viridis_pal()(100), zlim=c(0,1), main=c("JSS mean","JSS hmean","JSS min")[i])
}

