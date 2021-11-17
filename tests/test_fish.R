library(rfish)
library(tidyverse)
a = 1:30
l = numeric(30)
f = new(Fish)

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")

simulate_ibm = function(pop, h){
  # h = 0.5
  pop$set_harvestProp(h)
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
}

fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
pop = new(Population, fish)
K_ibm = pop$calcK()
hvec = seq(0.1,0.8,0.1)
ovec_ibm = numeric(length(hvec))
ovec = ovec_ibm
for(i in 1:length(hvec)){
  res_ibm = simulate_ibm(pop, hvec[i])
  res = simulate(hvec[i], LF50, F)
  ovec_ibm[i] = mean(res_ibm[61:100], na.rm=T)/1e9
  ovec[i] = mean(res$summaries$SSB[61:100])/1e9
  # plot(y=res$summaries$SSB/1e9, x=res$summaries$year)
}
points(ovec~hvec, ylim=c(0,8), xlab="Harvest fraction", ylab="SSB", type="l", col="green4")
points(ovec_ibm~hvec, col="green3", pch=1, cex=1.5)
# 
# ssb_mean = numeric(10)
# ssb_var = numeric(10)
# nvec = seq(4,0.1,length.out=10)*1e8
# nfish_vec = numeric(10)
# for (iter in 1:10){
#   pop = new(Population)
#   pop$n = nvec[iter]
#   h = 0
#   pop$mort_fishing_mature = -log(1-h)
#   pop$mort_fishing_immature = 0
#   pop$init(100)
#   pop$get_state()
#   
#   nsteps = 100
#   nfish = numeric(nsteps)
#   nfish[1]=100
#   ssb = numeric(nsteps)
#   for (i in 2:nsteps){
#     ssb[i] = pop$update()
#     nfish[i] = pop$nfish()
#   }
#   
#   ssb
#   ssb_mean[iter] = mean(ssb[61:100]/1e9)
#   ssb_var[iter] = var(ssb[61:100]/1e9)
#   nfish_vec[iter] = nfish[100]
# }
# par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
# plot(ssb_mean~nfish_vec, log="xy", xlab = "Number of superfish", ylab = "Mean Equil SSB")
# plot(ssb_var~nfish_vec, log="xy", xlab = "Number of superfish", ylab = "Variance in Equil SSB")
# mod = lm(log(ssb_var)~log(nfish_vec))
# mtext(text = sprintf("m = %1.2f", mod$coefficients[2]), side = 1, line = -1.5, adj = 0.05)
# 
#   
# 
# 
# 
# 
# pop = new(Population)
# pop$init(100)
# f = new(Fish)
# x = 1:30
# y = numeric(30)
# for (i in 1:30){
#   f$set_age(i)
#   y[i] = f$length
# }
# plot(y~x)
# plot(y=sapply(y, FUN = pop$selectivity), x=y)
# plot(y=sapply(y, FUN = pop$selectivity), x=x)
# 
# # 
# # sapply(X = L, FUN = pop$selectivity)




