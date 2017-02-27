source("./R/model_functions.R")

library(lhs)

# Explore model output

# Parameter values
betaA1.hcw2p <- 0.195 #patient S to A transmission rate (from hcw to patient) in setting 1 etc
betaA1.p2hcw <- 1.95 #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
betaB1.hcw2p <- 0.1  #patient S to A transmission rate (from hcw to patient) in setting 1 etc
betaB1.p2hcw <- 1.0  #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
betaA2 <- 0.00186
betaB2 <- 0.0032
betaA3 <- 0.00186
betaB3 <- 0.0032
theta <- 0  # Assumed lowered transmission when already colonised with one of the strains
betaAB1.hcw2p<- betaB1.hcw2p*theta  #A to B replace rate in setting 1 etc
betaBA1.hcw2p <- betaA1.hcw2p*theta 
betaAB2 <-betaB2*theta  
betaBA2 <-betaA2*theta
betaAB3 <-betaB3*theta
betaBA3 <-betaA3*theta
muA1<- 0.0025 # A clearance rate in setting 1 etc
muB1<- 0.025
muA2<- 0.0025
muB2<- 0.0025
muA3<- 0.0025
muB3<- 0.0025
mS1<-0.1 # hospital discharge rate for those in S1
mA1<-0.1 # etc
mB1<-0.1
rho<-20
n.hcw<-100
hhfreq<- 0.4 # hand hygiene frequency (from 0 to 0.9999) [and prob of hand hygiene after a contact]
c<-10 # contact rate: number of contacts a patient on average gets per day 
f<-1
f23<-10000/100000 # N2 over N3
f32<-1

# initial values
S1.t0<-999 # initial value for S1 etc
A1.t0<-1
B1.t0<-0
S2.t0<-5000 # initial value for S2 etc
A2.t0<-0
B2.t0<-0
A3.t0<-0
B3.t0<-1
S3.t0<-100000 - S2.t0 - S1.t0 - A1.t0 - B3.t0 # initial value for S3 etc
H.A1.t0<-0 # initial levels of HCW carriage
H.B1.t0 <-0
total.t0=sum(c(S1.t0,A1.t0,B1.t0,S2.t0,A2.t0,B2.t0,S3.t0,A3.t0,B3.t0))
time.int = 365000
newhhfreq = 0.5

R0etc<-calcR0.fqr(S1.t0,A1.t0, B1.t0,S2.t0,A2.t0, B2.t0,S3.t0,A3.t0, B3.t0,mS1,mA1,mB1,rho,muA1,muA2,muA3,betaA1.hcw2p, betaA1.p2hcw,betaA2,betaA3,n.hcw,hhfreq,c,f23,f32)
R0etc
R0Betc<-calcR0B.fqr(S1.t0,A1.t0, B1.t0,S2.t0,A2.t0, B2.t0,S3.t0,A3.t0, B3.t0,mS1,mA1,mB1,rho,muB1,muB2,muB3,betaB1.hcw2p, betaB1.p2hcw,betaB2,betaB3,n.hcw,hhfreq,c,f23,f32)
R0Betc

par(mfrow=c(1,1))
manipulate(
  plot.model(S.betaA1.hcw2p,S.betaA1.p2hcw ,S.betaA2,S.betaA3,S.betaB1.hcw2p, S.betaB1.p2hcw,S.betaB2,S.betaB3, S.hhfreq, S.time.int,S.newhhfreq,time=40000),
  S.betaB2=slider(0,0.02,0.0032,step=0.00001),
  S.betaB3=slider(0,0.02,0.0032,step=0.00001),
  S.newhhfreq=slider(0,0.99,0.4,step=0.01),
  S.betaA1.hcw2p=slider(0,.5,0.195,step=.001),
  S.betaA1.p2hcw=slider(0,5,1.95,step=.01), 
  S.betaA2=slider(0,0.1,.00186,step=.00001),
  S.betaA3=slider(0,0.1,.00186,step=.00001),
  S.betaB1.hcw2p=slider(0,.4,.1,step=.001),
  S.betaB1.p2hcw=slider(0,4,1.0,step=0.01),
  S.hhfreq=slider(0,0.99,.4,step=0.01),
  S.time.int=slider(100,30000,25000,step=100)
) 

par(mfrow=c(1,1))
manipulate(
  plot.model.inc(S.betaA1.hcw2p,S.betaA1.p2hcw,S.betaA2,S.betaA3,S.betaB1.hcw2p, S.betaB1.p2hcw,S.betaB2,S.betaB3, S.hhfreq,S.time.int, S.newhhfreq,time=60000),
  S.betaB2=slider(0,0.02,0.0032,step=0.00001),
  S.betaB3=slider(0,0.02,0.0032,step=0.00001),
  S.newhhfreq=slider(0,0.99,0.4,step=0.01),
  S.betaA1.hcw2p=slider(0,.5,0.195,step=.001),
  S.betaA1.p2hcw=slider(0,5,1.95,step=.01), 
  S.betaA2=slider(0,0.1,.00186,step=.00001),
  S.betaA3=slider(0,0.1,.00186,step=.00001),
  S.betaB1.hcw2p=slider(0,.4,.1,step=.001),
  S.betaB1.p2hcw=slider(0,4,1.0,step=0.01),
  S.hhfreq=slider(0,0.99,.4,step=0.01),
  S.time.int=slider(100,30000,25000,step=100)
) 

