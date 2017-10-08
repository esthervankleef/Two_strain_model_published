source("./R/model_functions.R")

library(lhs)
script.name = "hosp_com"
  
sens.values = read.csv("./Sensitivity_values/Values_sens_an.csv",sep=";")
sens.values.bi.mix = read.csv("./Sensitivity_values/Values_sens_bacinf_mix.csv",sep=";")

names(sens.values) = c("fh","betaA_hcw2p","betaA2","betaB_hcw2p","betaB2","R0A","R0B","baseA","baseB")

# Explore model output

# Parameter values
betaA1.hcw2p <- 0.187 #patient S to A transmission rate (from hcw to patient) in setting 1 etc
betaA1.p2hcw <- 1.87 #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
betaB1.hcw2p <- 0.1  #patient S to A transmission rate (from hcw to patient) in setting 1 etc
betaB1.p2hcw <- 1.0  #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
betaA2 <- 0.00212
betaB2 <- 0.0032
betaA3 <- 0.00212
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

# par(mfrow=c(1,1))
# manipulate(
#   plot.model(S.betaA1.hcw2p,S.betaA1.p2hcw ,S.betaA2,S.betaA3,S.betaB1.hcw2p, S.betaB1.p2hcw,S.betaB2,S.betaB3, S.hhfreq, S.time.int,S.newhhfreq,time=60000),
#   S.betaA1.hcw2p=slider(0,.5,0.187,step=.001),
#   S.betaA1.p2hcw=slider(0,5,1.87,step=.01), 
#   S.betaA2=slider(0,0.1,.00212,step=.00001),
#   S.betaA3=slider(0,0.1,.00212,step=.00001),
#   S.betaB1.hcw2p=slider(0,.4,.1,step=.001),
#   S.betaB1.p2hcw=slider(0,4,1.0,step=0.01),
#   S.betaB2=slider(0,0.02,0.0032,step=0.00001),
#   S.betaB3=slider(0,0.02,0.0032,step=0.00001),
#   S.newhhfreq=slider(0,0.99,0.5,step=0.01),
#   S.hhfreq=slider(0,0.99,.4,step=0.01),
#   S.time.int=slider(100,50000,45000,step=100)
# ) 
# 
# par(mfrow=c(1,1))
# manipulate(
#   plot.model.inc(S.betaA1.hcw2p,S.betaA1.p2hcw,S.betaA2,S.betaA3,S.betaB1.hcw2p, S.betaB1.p2hcw,S.betaB2,S.betaB3, S.hhfreq,S.time.int, S.newhhfreq,time=60000),
#   S.betaA1.hcw2p=slider(0,.5,0.187,step=.001),
#   S.betaA1.p2hcw=slider(0,5,1.87,step=.01), 
#   S.betaA2=slider(0,0.1,.00212,step=.00001),
#   S.betaA3=slider(0,0.1,.00212,step=.00001),
#   S.betaB1.hcw2p=slider(0,.4,.1,step=.001),
#   S.betaB1.p2hcw=slider(0,4,1.0,step=0.01),
#   S.betaB2=slider(0,0.02,0.0032,step=0.00001),
#   S.betaB3=slider(0,0.02,0.0032,step=0.00001),
#   S.newhhfreq=slider(0,0.99,0.5,step=0.01),
#   S.hhfreq=slider(0,0.99,.4,step=0.01),
#   S.time.int=slider(100,50000,45000,step=100)
# ) 

newhhfreq.vec<-c(.2,.25,.3,.35,.4,.45, 0.5, 0.55,.6) # 

l<-length(newhhfreq.vec)

matA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matB1<-matrix(rep(NA,73001*l),ncol=l)
matA2and3<-matrix(rep(NA,73001*l),ncol=l)
matB2and3<-matrix(rep(NA,73001*l),ncol=l)
matA.all<-matrix(rep(NA,73001*l),ncol=l)
matB.all<-matrix(rep(NA,73001*l),ncol=l)

matIA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matIB1<-matrix(rep(NA,73001*l),ncol=l)
matIA2and3<-matrix(rep(NA,73001*l),ncol=l)
matIB2and3<-matrix(rep(NA,73001*l),ncol=l)

matdIA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
matdIB1<-matrix(rep(NA,73001*l),ncol=l)
matdIA2and3<-matrix(rep(NA,73001*l),ncol=l)
matdIB2and3<-matrix(rep(NA,73001*l),ncol=l)

# STORE PREVALENCE & INCIDENCE - SA
for(i in 1:length(newhhfreq.vec)){
  print(newhhfreq.vec[i])
  newhhfreq<-newhhfreq.vec[i]
  time.int<-365000
  modpars <- c(betaA1.hcw2p, betaA1.p2hcw,betaB1.hcw2p,betaB1.p2hcw,betaA2,betaB2,betaA3,betaB3,betaAB1.hcw2p,betaBA1.hcw2p,betaAB2,betaBA2,betaAB3,betaBA3,muA1,muB1,muA2,muB2,muA3,muB3,mS1,mA1,mB1,rho, n.hcw,hhfreq,c, time.int,newhhfreq,f23,f32)   
  mod.init <-c(S1.t0,A1.t0,B1.t0,S2.t0,A2.t0,B2.t0,S3.t0,A3.t0,B3.t0,H.A1.t0,H.B1.t0,total.t0)
  mod.t <- seq(0,730000,by=10) 
  mod.sol <- lsoda(mod.init,mod.t,mod.dyn,modpars)
  N1<-S1.t0+A1.t0+B1.t0
  N2<-S2.t0+A2.t0+B2.t0
  N3<-S3.t0+A3.t0+B3.t0
  TIME <- mod.sol[,1] 
  A1 <- mod.sol[,3]/N1 
  B1 <- mod.sol[,4]/N1 
  A2and3<- (mod.sol[,6]+mod.sol[,9])/(N2+N3)
  B2and3<- (mod.sol[,7]+mod.sol[,10])/(N2+N3)
  A.all<-(mod.sol[,3]+mod.sol[,6]+mod.sol[,9])/(N1+N2+N3)
  B.all<-(mod.sol[,4]+mod.sol[,7]+mod.sol[,10])/(N1+N2+N3)
  
  matA1[,i]<-A1
  matB1[,i]<-B1
  matA2and3[,i]<-A2and3
  matB2and3[,i]<-B2and3
  matA.all[,i] <- A.all
  matB.all[,i] <- B.all
  
  # INCIDENCE
  IA1 <- mod.sol[,14] 
  IB1 <- mod.sol[,17] 
  IA2and3<- (mod.sol[,15]+mod.sol[,16])
  IB2and3<- (mod.sol[,18]+mod.sol[,19])
  matIA1[,i]<-IA1
  matIB1[,i]<-IB1
  matIA2and3[,i]<-IA2and3
  matIB2and3[,i]<-IB2and3
  
  # INCIDENCE PER DAY
  IdA1 <- mod.sol[,20] 
  IdB1 <- mod.sol[,23] 
  IdA2and3<- (mod.sol[,21]+mod.sol[,22])
  IdB2and3<- (mod.sol[,24]+mod.sol[,25])
  matdIA1[,i]<-IdA1
  matdIB1[,i]<-IdB1
  matdIA2and3[,i]<-IdA2and3
  matdIB2and3[,i]<-IdB2and3
} 

# Run model to find correct values for sensitivity analysis
time = 30000
time.int = time-365
lhs = randomLHS(1,2)
lhs[,1] = c(0.215)
lhs[,2] = c(0.00234)
out = data.frame(matrix(nrow=1, ncol=6))
R0.find = find.R0(lhs,out,time,time.int,"B",hhfreq = 0.4) 

# Run sensitivity analysis hospital adaptation 
t=730000
time.int<-365000
sets = length(sens.values[,1])
sets2=length(unique(sens.values$betaB_hcw2p[!is.na(sens.values$betaB_hcw2p)]))

sens.out = frac.h.sens(sets,sets2,sens.values,t,time.int,newhhfreq.vec = 0.5)
sens.out

# Run sensitivity analysis bacterial inference
sets = length(sens.values.bi.mix$bac_inf_half)
sets2=length(unique(sens.values$betaB_hcw2p[!is.na(sens.values$betaB_hcw2p)]))

sens.out = frac.h.sens(sets,sets2,sens.values,t,time.int,newhhfreq.vec = 0.5)
sens.out

# Run sensitivity analysis bacterial inference and mixing
sens.scen<-c("bi1","bi2","mix1","mix2") # 
l<-length(sens.scen)

sensA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
sensB1<-matrix(rep(NA,73001*l),ncol=l)
sensA2and3<-matrix(rep(NA,73001*l),ncol=l)
sensB2and3<-matrix(rep(NA,73001*l),ncol=l)
sensA.all<-matrix(rep(NA,73001*l),ncol=l)
sensB.all<-matrix(rep(NA,73001*l),ncol=l)

sensIA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
sensIB1<-matrix(rep(NA,73001*l),ncol=l)
sensIA2and3<-matrix(rep(NA,73001*l),ncol=l)
sensIB2and3<-matrix(rep(NA,73001*l),ncol=l)

sensdIA1<-matrix(rep(NA,73001*l),ncol=l) # 73001 is time
sensdIB1<-matrix(rep(NA,73001*l),ncol=l)
sensdIA2and3<-matrix(rep(NA,73001*l),ncol=l)
sensdIB2and3<-matrix(rep(NA,73001*l),ncol=l)

for(i in 1:length(sens.scen)){
  print(sens.scen[i])
  #time.int<-365000
  #modpars <- c(betaA1.hcw2p, betaA1.p2hcw,betaB1.hcw2p,betaB1.p2hcw,betaA2,betaB2,betaA3,betaB3,betaAB1.hcw2p,betaBA1.hcw2p,betaAB2,betaBA2,betaAB3,betaBA3,muA1,muB1,muA2,muB2,muA3,muB3,mS1,mA1,mB1,rho, n.hcw,hhfreq,c, time.int,newhhfreq,f23,f32)   
  modpars <- sens.values.bi.mix[,2+i]
  mod.init <-c(S1.t0,A1.t0,B1.t0,S2.t0,A2.t0,B2.t0,S3.t0,A3.t0,B3.t0,H.A1.t0,H.B1.t0,total.t0)
  mod.t <- seq(0,730000,by=10) 
  mod.sol <- lsoda(mod.init,mod.t,mod.dyn,modpars)
  N1<-S1.t0+A1.t0+B1.t0
  N2<-S2.t0+A2.t0+B2.t0
  N3<-S3.t0+A3.t0+B3.t0
  TIME <- mod.sol[,1] 
  A1 <- mod.sol[,3]/N1 
  B1 <- mod.sol[,4]/N1 
  A2and3<- (mod.sol[,6]+mod.sol[,9])/(N2+N3)
  B2and3<- (mod.sol[,7]+mod.sol[,10])/(N2+N3)
  A.all<-(mod.sol[,3]+mod.sol[,6]+mod.sol[,9])/(N1+N2+N3)
  B.all<-(mod.sol[,4]+mod.sol[,7]+mod.sol[,10])/(N1+N2+N3)
  
  sensA1[,i]<-A1
  sensB1[,i]<-B1
  sensA2and3[,i]<-A2and3
  sensB2and3[,i]<-B2and3
  sensA.all[,i] <- A.all
  sensB.all[,i] <- B.all
  
  # INCIDENCE
  IA1 <- mod.sol[,14] 
  IB1 <- mod.sol[,17] 
  IA2and3<- (mod.sol[,15]+mod.sol[,16])
  IB2and3<- (mod.sol[,18]+mod.sol[,19])
  sensIA1[,i]<-IA1
  sensIB1[,i]<-IB1
  sensIA2and3[,i]<-IA2and3
  sensIB2and3[,i]<-IB2and3
  
  # INCIDENCE PER DAY
  IdA1 <- mod.sol[,20] 
  IdB1 <- mod.sol[,23] 
  IdA2and3<- (mod.sol[,21]+mod.sol[,22])
  IdB2and3<- (mod.sol[,24]+mod.sol[,25])
  sensdIA1[,i]<-IdA1
  sensdIB1[,i]<-IdB1
  sensdIA2and3[,i]<-IdA2and3
  sensdIB2and3[,i]<-IdB2and3
} 

  
# Save dataset
savename <- paste0("./Output/",script.name,".Rda")
save(matA.all,matA1,matA2and3,matB.all,matB1,matB2and3,
     matIA1,matIA2and3,matIB1,matIB2and3,
     matdIA1,matdIA2and3,matdIB1,matdIB2and3, file = savename)

savename.sens <- paste0("./Output/",script.name,".sens",".Rda")
save(sens.out,file=savename.sens)

savename.sens.bi.mix <- paste0("./Output/",script.name,".sens_bi_mix",".Rda")
save(sensA.all,sensA1,sensA2and3,sensB.all,sensB1,sensB2and3,
     sensIA1,sensIA2and3,sensIB1,sensIB2and3,
     sensdIA1,sensdIA2and3,sensdIB1,sensdIB2and3, file = savename.sens.bi.mix)

