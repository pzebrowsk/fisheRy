rm(list=ls())
library(rfish)
library(tidyverse)


plot(x=1, y=NA, xlim=c(0,21), ylim=c(0,50), xlab="age", ylab = "length")

mat_ages = numeric(0)
for (f in 1:500){
  #cat("f = ", f, "\n")
  fish.par = new(FishParams)
  fish = new(Fish)
  fish.par$flag = 1
  fish$par = fish.par
  
  dat = data.frame(t_birth=0,	age=0,	isMature=0,	isAlive=0,	length=0,	weight=0)
  dat = rbind(dat, fish$get_state())
  dat = dat[-1,]
  for (i in 1:20){
    #cat("rmat prob: age = ", fish$age, ", p = ", 1 / (1 + exp( -(fish$length - (1.604122 * fish$age + 18.399575)) / 1.93031)), "\n")
    #cat("f = ", f, "\n")
    fish$updateMaturity()
    fish$grow()
    fish$set_age(fish$age+1)
    dat = rbind(dat, fish$get_state())
  }
  mat_ages = c(mat_ages, which(dat$isMature>0)[1])
  points(dat$length~dat$age, type="l", col=rainbow(500)[f])
}
table(mat_ages)

plot(x=1, y=NA, xlim=c(0,35), ylim=c(0,1), xlab="length", ylab = "mat prob")
for (age in 1:15){
  curve(1 / (1 + exp( -(x - (1.604122 * age + 18.399575)) / 1.93031)), 0, 35, col=rainbow(15)[age], add=T)
}

plot(x=1, y=NA, xlim=c(0,15), ylim=c(0,1), xlab="age", ylab = "mat prob")
for (length in 1:35){
  curve(1 / (1 + exp( -(length - (1.604122 * x + 18.399575)) / 1.93031)), 0, 35, col=rainbow(35)[length], add=T)
}

mat = matrix(nrow=15, ncol=35)
for (age in 1:15){
  for (length in 1:35){
    mat[age, length] = 1 / (1 + exp( -(length - (1.604122 * age + 18.399575)) / 1.93031))
  }   
}
image(x=1:15, y=1:35, z=(mat), col=scales::viridis_pal()(100), xlab="age", ylab="length")


