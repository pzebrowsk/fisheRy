rm(list=ls())
library(fisheRy)
library(tidyverse)

    
plot(x=1, y=NA, xlim=c(0,31), ylim=c(0,200), xlab="age", ylab = "length")

# NEW MODEL
mat_ages = numeric(0)
mort = data.frame(a=0, m=0, l=0)
N = 100
for (f in 1:N){
  #cat("f = ", f, "\n")
  fish.par = new(FishParams)
  fish.par$s0 = 0.09637
  # fish.par$Bhalf_growth = 1e11
  
  fish = new(Fish)
  fish.par$flag = 1
  fish.par$growth_model = 1
  fish$par = fish.par
  fish$init(1.93e3, 5.61)
    
  dat = data.frame(t_birth=0,	age=0,	isMature=0,	isAlive=0,	length=0,	weight=0)
  dat = rbind(dat, fish$get_state())
  dat = dat[-1,]
  for (i in 1:30){
    #cat("rmat prob: age = ", fish$age, ", p = ", 1 / (1 + exp( -(fish$length - (1.604122 * fish$age + 18.399575)) / 1.93031)), "\n")
    #cat("f = ", f, "\n")
    dat = rbind(dat, fish$get_state())
    mort = rbind(mort, c(a=fish$age,  fish$naturalMortalityRate(), fish$length))
    fish$updateMaturity()
    fish$grow(1.93e3, 5.61)
    fish$set_age(fish$age+1)
  }
  mat_ages = c(mat_ages, which(dat$isMature>0)[1])
  points(dat$length~dat$age, type="l", col=rainbow(N)[f])
}
table(mat_ages)

# OLD MODEL
N=1
for (f in 1:N){
  #cat("f = ", f, "\n")
  fish.par = new(FishParams)
  fish.par$s0 = 0.09637
  # fish.par$Bhalf_growth = 1e11

  fish = new(Fish)
  fish.par$flag = 1
  fish.par$growth_model = 0
  fish$par = fish.par
  fish$init(1.93e3, 5.61)
  
  dat = data.frame(t_birth=0,	age=0,	isMature=0,	isAlive=0,	length=0,	weight=0)
  dat = rbind(dat, fish$get_state())
  dat = dat[-1,]
  for (i in 1:30){
    #cat("rmat prob: age = ", fish$age, ", p = ", 1 / (1 + exp( -(fish$length - (1.604122 * fish$age + 18.399575)) / 1.93031)), "\n")
    #cat("f = ", f, "\n")
    dat = rbind(dat, fish$get_state())
    mort = rbind(mort, c(a=fish$age,  fish$naturalMortalityRate(), fish$length))
    fish$updateMaturity()
    fish$grow(1.93e3, 5.61)
    fish$set_age(fish$age+1)
  }
  mat_ages = c(mat_ages, which(dat$isMature>0)[1])
  points(dat$length~dat$age, type="l", lwd=2, col="black")
}



source("tests/ref/parameters.cod.R")
M<-c(M,rep(M[length(M)],times=ages-length(M))); names(M)<-1:ages
with(mort %>% filter(a<=20 & a>1), plot(m~a, ylim=c(0,1)))
points(M, type="l", col="blue")

length_age = mort %>% group_by(a) %>% summarize(l = mean(l))

# plot(x=1, y=NA, xlim=c(0,35), ylim=c(0,1), xlab="length", ylab = "mat prob")
# for (age in 1:15){
#   curve(1 / (1 + exp( -(x - (1.604122 * age + 18.399575)) / 1.93031)), 0, 35, col=rainbow(15)[age], add=T)
# }
# 
# plot(x=1, y=NA, xlim=c(0,15), ylim=c(0,1), xlab="age", ylab = "mat prob")
# for (length in 1:35){
#   curve(1 / (1 + exp( -(length - (1.604122 * x + 18.399575)) / 1.93031)), 0, 35, col=rainbow(35)[length], add=T)
# }

# mat = matrix(nrow=20, ncol=160)
# for (age in 1:20){
#   for (length in 1:160){
#     mat[age, length] = 1 / (1 + exp( -(length - (1.604122 * age + 18.399575)) / 1.93031))
#   }   
# }
# image(x=1:15, y=1:35, z=(mat), col=scales::viridis_pal()(100), xlab="age", ylab="length")
# 

