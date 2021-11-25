rm(list=ls())
library(rfish)
library(tidyverse)

source("tests/ref/parameters.cod.R")
source("tests/ref/simulator.7.R")

fish.par = new(FishParams)
fish = new(Fish)
fish.par$flag = 1
fish$par = fish.par
pop = new(Population, fish)
K_ibm = pop$calcK()
# hvec = seq(0.01,0.45,0.02)
hvec = seq(0.01,0.5,0.02)

sim = new(Simulator)
sim$simulate(pop, 0.45, 200, T)

res_ibm_full = sim$simulate_multi(pop, hvec, 50, F)
arr = array(data=res_ibm_full, dim=c(50, 4, length(hvec)))

sim = new(Simulator)
d = sim$max_avg_utils(c(length(hvec),4,50), res_ibm_full)

par(mfrow=c(1,2), mar=c(4,4,1,1))
matplot(y=matrix(data=d, byrow=T, ncol=4), x=hvec, lty=1, type="l", col=c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"), lwd=3, ylim=c(0,1))
abline(h=0, col="grey")
plot(1,1, cex=0.01, xlab = "", ylab = "", axes = F)
legend(x = 0.7, y = 1, legend = c("ssb", "yield", "employment", "profit"), fill = c("darkgreen", "darkgoldenrod1", "dodgerblue3", "coral1"))

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
  
  emp.max = max(c(res_ibm$employment, res$summaries$Ereal))
  plot(y=res_ibm$employment, x=seq(1,nsteps,1), ylab="Employment (person-years)", xlab="Year", col="red", ylim=c(0,emp.max))
  points(y=res$summaries$Ereal, x=res$summaries$year, type="l")

  profit.max = max(c(res_ibm$profit/1e9, res$summaries$P.tot/1e9))
  profit.min = min(c(res_ibm$profit/1e9, res$summaries$P.tot/1e9))
  plot(y=res_ibm$profit/1e9, x=seq(1,nsteps,1), ylab="Profit (Billions NOK)", xlab="Year", col="red", ylim=c(profit.min,profit.max))
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




