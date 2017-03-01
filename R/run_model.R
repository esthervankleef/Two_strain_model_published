source("./R/model_functions.R")

library(lhs)
script.name.mrsa = "hosp_com_mrsa"

sens.values = read.csv("./Sensitivity_values/Values_sens_an.csv",sep=";")

names(sens.values) = c("fh","betaA_hcw2p","betaA2","betaB_hcw2p","betaB2","R0A","R0B","baseA","baseB")

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

# Fraction of patients in community
par(mfrow=c(1,1))
manipulate(
  f.trans(S.betaA1.hcw2p,S.betaA1.p2hcw,S.betaA2,S.betaA3,S.betaB1.hcw2p,betaB1.p2hcw,S.betaB2,S.betaB3,time=40000),
  S.betaA1.hcw2p=slider(0,.5,0.195,step=.001),
  S.betaA1.p2hcw=slider(0,5,1.95,step=.01), 
  S.betaA2=slider(0,0.1,.00145,step=.0001),
  S.betaA3=slider(0,0.1,.00145,step=.0001),
  S.betaB1.hcw2p=slider(0,.1,.05,step=.001),
  S.betaB1.p2hcw=slider(0,1,.5,step=0.01),
  S.betaB2=slider(0,0.02,0.0026,step=0.0001),
  S.betaB3=slider(0,0.02,0.0026,step=0.0001)
)


newhhfreq.vec<-c(.2,.25,.3,.35,.4,.45, 0.5, 0.55,.6) # 

l<-length(newhhfreq.vec)

matA1.mrsa<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matB1.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matA2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matB2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matA.all.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matB.all.mrsa<-matrix(rep(NA,73001*l),ncol=l)

matIA1.mrsa<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matIB1.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matIA2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matIB2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)

matdIA1.mrsa<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matdIB1.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matdIA2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)
matdIB2and3.mrsa<-matrix(rep(NA,73001*l),ncol=l)

# STORE PREVALENCE & INCIDENCE - SA
for(i in 1:length(newhhfreq.vec)){
  print(newhhfreq.vec[i])
  newhhfreq<-newhhfreq.vec[i]
  time.int<-365000
  modpars <- c(betaA1.hcw2p, betaA1.p2hcw,betaB1.hcw2p,betaB1.p2hcw,betaA2,betaB2,betaA3,betaB3,betaAB1.hcw2p,betaBA1.hcw2p,betaAB2,betaBA2,betaAB3,betaBA3,muA1,muB1,muA2,muB2,muA3,muB3,mS1,mA1,mB1,rho, n.hcw,hhfreq,c, time.int,newhhfreq,f23,f32)   
  mod.init <-c(S1.t0,A1.t0,B1.t0,S2.t0,A2.t0,B2.t0,S3.t0,A3.t0,B3.t0,H.A1.t0,H.B1.t0,total.t0)
  mod.t <- seq(0,730000,by=10) 
  mod.sol <- lsoda(mod.init,mod.t,modCDI.dyn,modpars)
  N1<-S1.t0+A1.t0+B1.t0
  N2<-S2.t0+A2.t0+B2.t0
  N3<-S3.t0+A3.t0+B3.t0
  TIME <- mod.sol[,1] 
  A1 <- mod.sol[,3]/N1 
  B1 <- mod.sol[,4]/N1 
  A2and3<- (mod.sol[,6]+mod.sol[,9])/(N2+N3)
  B2and3<- (mod.sol[,7]+mod.sol[,10])/(N2+N3)
  A.all.mrsa<-(mod.sol[,3]+mod.sol[,6]+mod.sol[,9])/(N1+N2+N3)
  B.all.mrsa<-(mod.sol[,4]+mod.sol[,7]+mod.sol[,10])/(N1+N2+N3)
  
  matA1.mrsa[,i]<-A1
  matB1.mrsa[,i]<-B1
  matA2and3.mrsa[,i]<-A2and3
  matB2and3.mrsa[,i]<-B2and3
  matA.all.mrsa[,i] <- A.all.mrsa
  matB.all.mrsa[,i] <- B.all.mrsa
  
  # INCIDENCE
  IA1.mrsa <- mod.sol[,14] 
  IB1.mrsa <- mod.sol[,17] 
  IA2and3.mrsa<- (mod.sol[,15]+mod.sol[,16])
  IB2and3.mrsa<- (mod.sol[,18]+mod.sol[,19])
  matIA1.mrsa[,i]<-IA1.mrsa
  matIB1.mrsa[,i]<-IB1.mrsa
  matIA2and3.mrsa[,i]<-IA2and3.mrsa
  matIB2and3.mrsa[,i]<-IB2and3.mrsa
  
  # INCIDENCE PER DAY
  IdA1.mrsa <- mod.sol[,20] 
  IdB1.mrsa <- mod.sol[,23] 
  IdA2and3.mrsa<- (mod.sol[,21]+mod.sol[,22])
  IdB2and3.mrsa<- (mod.sol[,24]+mod.sol[,25])
  matdIA1.mrsa[,i]<-IdA1.mrsa
  matdIB1.mrsa[,i]<-IdB1.mrsa
  matdIA2and3.mrsa[,i]<-IdA2and3.mrsa
  matdIB2and3.mrsa[,i]<-IdB2and3.mrsa
  
  
} 

# Run model to find correct values for sensitivity analysis
time = 30000
time.int = time-365
lhs = randomLHS(1,2)
lhs[,1] = c(0.215)
lhs[,2] = c(0.00234)
out = data.frame(matrix(nrow=1, ncol=6))
R0.find = find.R0(lhs,out,time,time.int,"B",hhfreq = 0.4) 

# Run sensitivity analysis
t=730000
time.int<-365000
sets = length(sens.values[,1])
sets2=length(unique(sens.values$betaB_hcw2p[!is.na(sens.values$betaB_hcw2p)]))

sens.out = frac.h.sens(sets,sets2,sens.values,t,time.int,newhhfreq.vec = 0.5)
sens.out

# Run sensitivity analysis with different bacterial interference
t=730000
time.int<-365000
sets = length(sens.values[,1])
sets2=length(unique(sens.values$betaB_hcw2p[!is.na(sens.values$betaB_hcw2p)]))

sens.out.comp = frac.h.sens(sets,sets2,sens.values,t,time.int,newhhfreq.vec = 0.5)
sens.out.comp

# Check single admission numbers
sets = length(sens.values[,1])
sets2=length(unique(sens.values$betaB_hcw2p[!is.na(sens.values$betaB_hcw2p)]))
single.aR0 = single.adm.r0(sets,sets2,sens.values)
single.aR0

single.aR0$h_higher = ifelse(single.aR0$RA1.1>single.aR0$RB1.1,1,0)

# Save dataset
savename.mrsa <- paste0("./Output/",script.name.mrsa,".Rda")
save(matA.all.mrsa,matA1.mrsa,matA2and3.mrsa,matB.all.mrsa,matB1.mrsa,matB2and3.mrsa,
     matIA1.mrsa,matIA2and3.mrsa,matIB1.mrsa,matIB2and3.mrsa,
     matdIA1.mrsa,matdIA2and3.mrsa,matdIB1.mrsa,matdIB2and3.mrsa, file = savename.mrsa)

savename.mrsa.sens <- paste0("./Output/",script.name.mrsa,".sens",".Rda")
save(sens.out,file=savename.mrsa.sens)

savename.mrsa.admR0 <- paste0("./Output/",script.name.mrsa,".admR0",".Rda")
save(single.aR0,file=savename.mrsa.admR0)

savename.mrsa.sens.comp <- paste0("./Output/",script.name.mrsa,".sens.comp",".Rda")
save(sens.out.comp,file=savename.mrsa.sens.comp)
