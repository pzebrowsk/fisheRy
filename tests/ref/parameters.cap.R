### parameters
# implementation parameters and age ranges
ages<-6			# maximum age - AFWG2014 Table 9.1
years<-100		# simulation years
skip<-50		# transient years to be skipped
N.ini<-1e8		# initial numbers at age 1

# biological parameters
r0<-2503.285		# maximum age 1 recruits per kg SSB (=1e3*762/(71.4+2.33*100), from Mikkelsen & Pedersen (2004) Table 1, model 1-1-1, for relatively low herring abundance)
max.R<-10*700*1e9 # cap for R - 10 times the highest R on record
Bhalf<-662*1e6	# half-saturation constant in recruitment ***tuned***
Linf<-24.7		# VB: asymptotic length - estimated from IMR data (autumn)
a0<- -0.6387	# accounts for autumn timing of the survey
k<-0.248		# 
theta<-0.001178	# weight-length: condition - estimated from IMR data
g<-3.515		# weight-length: allometric exponent - estimated from IMR data
# age-specific biological parameters
Mi<-c(0.0572,0.0572,0.0830,0.0830,0.0830,0.0830) # natural mortality: immature, per month, latest benchmark (AFWG2008)
M<-c(0.1713,0.1713,0.2489,0.2489,0.2489,0.2489) # natural mortality: mature, per season (3 month prior to spawning), latest benchmark (AFWG2008)
mat<-c(0,0.0099,0.1918,0.6840,0.6629,1) # maturity-at-age, from autumn survey data, adjusted to spring
sf<-1.362229		# recruitment noise, from residuals of the S-R model

# fishing parameters
Blim<-200e6		# escapement target in the end of the fishing season
F.max<-0.2998	# de facto fishing mortality when fishing is allowed
LF50<-14.5		# length at 50% selectivity - from Bjarte Bogstad pers. comm. L50=14.5cm, L99=16cm 
SF<-log(99)/(16-14.5)	# selection factor  - from Bjarte Bogstad pers. comm. L50=14.5cm, L99=16cm
b<-0.25			# concentration profile parameter
q<-4.65e-05		# scaling parameter relating to catchability and density ***tuned***

# employment parameters for CAPELIN
d.sea<-0.045	# FTE/vessel day
d.max<-30000	# maximum available FTE at sea
d.shr<-8.23e-7	# FTE/kg

# economic parameters for CAPELIN
price.sea<-1.81 # landing price NOK/kg
price.shr<-9.52 # selling price NOK/kg

s.sea<-1353000			# employment cost sea NOK/FTE
s.shr<-367000			# employment cost shore NOK/FTE, , data avg(2001-2011)
fixedC.sea<-294443000	# fixed costs sea NOK (=average per unit * #units)
fixedC.shr<-12690757000	# fixed costs shore NOK? (=average per unit * #units), data avg(2001-2011)
varC.sea<-61000			# variable costs NOK/vessel day
scale.catch.all<-0.128	# percentage of total pelagic catch that is capelin - all years
scale.catch.cap<-0.192	# percentage of total pelagic catch that is capelin - capelin years