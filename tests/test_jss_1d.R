rm(list=ls())
library(rfish)
library(tidyverse)

#### MAIN ####

fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
fish$par$s0 = 0.09637
fish$par$Bhalf_growth = 100e11

pop_K = new(Population, fish)
pop_K$set_superFishSize(2e6)
K_ibm = pop_K$noFishingEquilibriate()

hist(pop_K$get_state()$age, ylab="Frequency", xlab="Age", breaks=21)
hist(pop_K$get_state()$length, ylab="Frequency", xlab="length", breaks=21)

# plot(K/1e6, log="y", type="l")
# points(K_ibm[-1])

sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop = new(Population, fish)
pop$set_superFishSize(2e6)

res_ibm_sq = sim$simulate(pop, 45, 0.41, 200, T)
#res_sq = simulate(0.41, 45, F)

# ### Continuation from Status quo (test)
# res_ib_csq = sim$simulate(pop, 45, 0.2, 200, F)
# with(rbind(res_ibm_sq, res_ib_csq), plot(ssb/1e9, type="l", ylab="SSB (MT)", xlab="Yr"))
# ###

pop_ref = pop

# hvec = seq(0.01,0.45,0.02)
hvec = seq(0.0,0.99,length.out = 20)
lfvec = seq(20,160,length.out = 20)

# res_ibm_full = sim$simulate_multi(pop, hvec, 50, F)

niter = 1
jss_arr = array(dim=c(3, length(hvec), niter))
full_arr = array(dim=c(c(50, length(hvec), 4, niter)))
for (iter in 1:niter){

  pop = pop_ref

  res_ibm_full = sim$simulate_multi(pop, hvec, 50, F)
  arr = res_ibm_full #array(data=res_ibm_full, dim=c(50, length(hvec), 4))
  full_arr[,,,iter] = arr

  #sim = new(Simulator) 
  d = sim$max_avg_utils(c(4,length(hvec),50), res_ibm_full)

  par(mfrow=c(1,2), mar=c(4,4,1,1))
  matplot(y=matrix(data=d, ncol=4), x=hvec, lty=1, type="l", col=c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"), lwd=3, ylim=c(0,1))
  abline(h=0, col="grey")
  plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
  legend(x = 0.7, y = 1, legend = c("ssb", "yield", "employment", "profit"), fill = c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"))

  #sim = new(Simulator)
  ss = sim$stakeholder_satisfaction(c(4,length(hvec),50), res_ibm_full)
  par(mfrow=c(1,2), mar=c(4,4,1,1))
  matplot(y=matrix(data=ss, ncol=5), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("grey", "darkgoldenrod1", "dodgerblue3", "coral1","darkgreen"))
  abline(h=0, col="grey")
  plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
  legend(x = 0.7, y = 1, legend = c("Industrial", "Artisanal", "Emp. policymakers", "Profit. policymakers", "Conservationists"), fill=c("grey", "darkgoldenrod1", "dodgerblue3", "coral1","darkgreen"))

  mat = matrix(data=ss, ncol=5)
  JSS_min = apply(mat, 1, min)
  JSS_mean = apply(mat, 1, mean)
  JSS_hmean = 1/apply(1/mat, 1, mean)

  jss_arr[1,,iter] = JSS_mean
  jss_arr[2,,iter] = JSS_hmean
  jss_arr[3,,iter] = JSS_min

  matplot(y=cbind(JSS_mean, JSS_hmean, JSS_min), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("darkgoldenrod1", "dodgerblue3", "darkgreen"))
  plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
  legend(x = 0.7, y = 1, legend = c("Mean", "Harmonic mean", "Min"), fill=c("darkgoldenrod1", "dodgerblue3","darkgreen"))
}

JSS_ensavg = apply(jss_arr, MARGIN = c(1,2), FUN = mean)
matplot(y=t(JSS_ensavg), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("darkgoldenrod1", "dodgerblue3", "darkgreen"))
plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
legend(x = 0.7, y = 1, legend = c("Mean", "Harmonic mean", "Min"), fill=c("darkgoldenrod1", "dodgerblue3","darkgreen"))



