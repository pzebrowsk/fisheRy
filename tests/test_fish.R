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

plot(table(pop_K$get_state()$age))
# plot(K/1e6, log="y", type="l")
# points(K_ibm[-1])

# hvec = seq(0.01,0.45,0.02)
hvec = seq(0.01,0.8,length.out = 20)
lfvec = seq(20,160,length.out = 20)

sim = new(Simulator, fish)
sim$setNaturalPopulation(pop_K)

pop = new(Population, fish)
pop$set_superFishSize(500)

res_ibm_sq = sim$simulate(pop, 45, 0.41, 200, T)
res_sq = simulate(0.41, 45, F)

pop_ref = pop

niter = 10
jss_arr = array(dim=c(3, length(hvec), length(lfvec), niter))
full_arr = array(dim=c(c(50, length(hvec), length(lfvec), 4, niter)))
for (iter in 1:niter){
  
  pop = pop_ref
  
  res_ibm_full = sim$simulate_multi_2d(pop, lfvec, hvec, 50, F)
  arr = array(data=res_ibm_full, dim=c(50, length(hvec), length(lfvec), 4))
  full_arr[,,,,iter] = arr
  
  sim = new(Simulator)
  d = sim$max_avg_utils_2d(c(4,length(lfvec),length(hvec),50), res_ibm_full)
  utils = array(data=d, dim=c(length(hvec),length(lfvec), 4))
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  for (i in 1:4){
    image(x=hvec, y=lfvec, z=utils[,,i], col=rainbow(start = 0, end = 0.75, n=100), zlim=c(0,1), main=c("ssb","yield","employment","profit")[i])
  }
  
  sim = new(Simulator)
  ss = sim$stakeholder_satisfaction_2d(c(4,length(lfvec),length(hvec),50), res_ibm_full)
  scs = array(data=ss, dim=c(length(hvec),length(lfvec), 5))
  par(mfrow=c(2,3), mar=c(4,4,4,1))
  for (i in 1:5){
    image(x=hvec, y=lfvec, z=scs[,,i], col=rainbow(start = 0, end = 0.75, n=100), zlim=c(0,1), main=c("ssb","yield","employment","profit")[i])
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


# Working 1D code
# niter = 10
# jss_arr = array(dim=c(3, length(hvec), niter))
# full_arr = array(dim=c(c(50, length(hvec), 4, niter)))
# for (iter in 1:niter){
# 
#   pop = pop_ref
#   
#   res_ibm_full = sim$simulate_multi(pop, hvec, 50, F)
#   arr = array(data=res_ibm_full, dim=c(50, length(hvec), 4))
#   full_arr[,,,iter] = arr
#     
#   sim = new(Simulator)
#   d = sim$max_avg_utils(c(4,length(hvec),50), res_ibm_full)
#   
#   par(mfrow=c(1,2), mar=c(4,4,1,1))
#   matplot(y=matrix(data=d, ncol=4), x=hvec, lty=1, type="l", col=c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"), lwd=3, ylim=c(0,1))
#   abline(h=0, col="grey")
#   plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
#   legend(x = 0.7, y = 1, legend = c("ssb", "yield", "employment", "profit"), fill = c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"))
#   
#   sim = new(Simulator)
#   ss = sim$stakeholder_satisfaction(c(4,length(hvec),50), res_ibm_full)
#   par(mfrow=c(1,2), mar=c(4,4,1,1))
#   matplot(y=matrix(data=ss, ncol=5), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("grey", "darkgoldenrod1", "dodgerblue3", "coral1","darkgreen"))
#   abline(h=0, col="grey")
#   plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
#   legend(x = 0.7, y = 1, legend = c("Industrial", "Artisanal", "Emp. policymakers", "Profit. policymakers", "Conservationists"), fill=c("grey", "darkgoldenrod1", "dodgerblue3", "coral1","darkgreen"))
#     
#   mat = matrix(data=ss, ncol=5)
#   JSS_min = apply(mat, 1, min)
#   JSS_mean = apply(mat, 1, mean)
#   JSS_hmean = 1/apply(1/mat, 1, mean)
#   
#   jss_arr[1,,iter] = JSS_mean
#   jss_arr[2,,iter] = JSS_hmean
#   jss_arr[3,,iter] = JSS_min
#   
#   matplot(y=cbind(JSS_mean, JSS_hmean, JSS_min), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("darkgoldenrod1", "dodgerblue3", "darkgreen"))
#   plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
#   legend(x = 0.7, y = 1, legend = c("Mean", "Harmonic mean", "Min"), fill=c("darkgoldenrod1", "dodgerblue3","darkgreen"))
# }
# 
# JSS_ensavg = apply(jss_arr, MARGIN = c(1,2), FUN = mean)
# matplot(y=t(JSS_ensavg), x=hvec, lty=1, type="l", lwd=3, ylim=c(0,1), col=c("darkgoldenrod1", "dodgerblue3", "darkgreen"))
# plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
# legend(x = 0.7, y = 1, legend = c("Mean", "Harmonic mean", "Min"), fill=c("darkgoldenrod1", "dodgerblue3","darkgreen"))


#### DEBUG ####

format_ibm = function(mat){
  d = data.frame(mat)
  colnames(d) = c("ssb", "yield", "employment", "profit")
  d
}

# ssb1 = sim$simulate(pop, 0, 200, T)
# ssb2 = sim$simulate(pop, 0.2, 100, F)
# plot(c(ssb1, ssb2), type="l")

plot_compare = function(res_ibm, res, nsteps){
  ssb.max = max(c(res_ibm$ssb/1e9, res$summaries$SSB/1e9))
  plot(y=res_ibm$ssb/1e9, x=seq(1,nsteps,1), ylab="SSB (MT)", xlab="Year", col="red", ylim=c(0,ssb.max))
  points(y=res$summaries$SSB/1e9, x=res$summaries$year, type="l")
  
  yield.max = max(c(res_ibm$yield/1e9, res$summaries$Y/1e9))
  plot(y=res_ibm$yield/1e9, x=seq(1,nsteps,1), ylab="Yield (MT)", xlab="Year", col="red", ylim=c(0,yield.max))
  points(y=res$summaries$Y/1e9, x=res$summaries$year, type="l")
  
  emp.max = max(c(res_ibm$employment.sea+res_ibm$employment.shore, res$summaries$D.tot))
  plot(y=res_ibm$employment.sea+res_ibm$employment.shore, x=seq(1,nsteps,1), ylab="Employment (person-years)", xlab="Year", col="red", ylim=c(0,emp.max))
  points(y=res$summaries$D.tot, x=res$summaries$year, type="l")

  profit.max = max(c((res_ibm$profit.sea+res_ibm$profit.shore)/1e9, res$summaries$P.tot/1e9))
  profit.min = min(c((res_ibm$profit.sea+res_ibm$profit.shore)/1e9, res$summaries$P.tot/1e9))
  plot(y=(res_ibm$profit.sea+res_ibm$profit.shore)/1e9, x=seq(1,nsteps,1), ylab="Profit (Billions NOK)", xlab="Year", col="red", ylim=c(profit.min,profit.max))
  points(y=res$summaries$P.tot/1e9, x=res$summaries$year, type="l")
  
}

par(mfrow = c(2,2), mar=c(4,4,1,1))
for(i in 1:length(hvec)){
#  res_ibm = sim$simulate(pop, hvec[i], 200, T)
  res = simulate(hvec[i], LF50, F)
  res_ibm = format_ibm(arr[,,i])
  plot_compare(res_ibm, res, 50)
}


dir = "~/Documents/fish_project/"
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to=dir)
plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, dir, row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0(dir, 1:length(sorted.png.names), ".png")
# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
file.rename(from=sorted.png.names, to=numbered.png.names)

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




