#######################################
# Two strain model
#######################################
rm(list=ls())
# Authors: B.S. Cooper with adjustments from E. van Kleef
# Date: 30 January 2017

library(deSolve)
require(manipulate)

# **********************************************************************************
# time.int is the the time of the hand hygiene intervention
# new.hhfreq is the hand hygiene level after the intervention.

modCDI.dyn <- function(t,var,par) { 
  # Basic model for two strains in competition allowing for a hand hygiene intervention
  # With transmission in the hospital and community 
  # And explicitly modelling HCW and hand-borne transmssion
  
  S1 <- var[1] # susceptible (not infected)  - setting 1(hospital)
  A1 <- var[2] # infected with A  - setting 1 (hospital)
  B1 <- var[3] # infected with B - setting 1 (hospital)
  S2 <- var[4] # susceptible (not infected)  - setting 2 (recently hospitalized)
  A2 <- var[5] # etc.
  B2 <- var[6]
  S3 <- var[7] # susceptible (not infected)  - setting 3 (not recently hospitalized)
  A3 <- var[8] #
  B3 <- var[9]
  H.A1<-var[10] #HCWs in pop 1 carrying A on their hands
  H.B1<-var[11] #HCWs in pop 2 carrying B on their hands
  Total<-var[12]
  
  # host vector model for hospital spread as we want to model hand hygiene  
  betaA1.hcw2p <- par[1] #patient S to A transmission rate (from hcw to patient) in setting 1 etc
  betaA1.p2hcw <- par[2] #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
  betaB1.hcw2p <- par[3] #patient S to A transmission rate (from hcw to patient) in setting 1 etc
  betaB1.p2hcw <- par[4] #hcw S to A transmission rate (from patient to hcw) in setting 1 etc
  
  betaA2 <- par[5]
  betaB2<- par[6]
  betaA3<- par[7]
  betaB3<- par[8]
  betaAB1.hcw2p <- par[9]  #A to B replace rate in setting 1 etc
  betaBA1.hcw2p <- par[10] 
  betaAB2 <- par[11]  # B to A replace rate in setting 1 
  betaBA2 <- par[12] 
  betaAB3 <- par[13] 
  betaBA3 <- par[14] 
  muA1<- par[15] # A clearance rate in setting 1 etc
  muB1<- par[16] 
  muA2<- par[17]
  muB2<- par[18]
  muA3<- par[19]
  muB3<- par[20]
  mS1<-par[21] # hospital discharge rate for those in S1
  mA1<-par[22] # etc
  mB1<-par[23] # etc
  rho<-par[24] #ratio of hospitalisation rate for those in pop2 to hosp rate for pop3 (generally >1)
  n.hcw<-par[25]  # number of HCWs
  hhfreq<-par[26] # hand hygiene frequency (from 0 to 99.99%) [and prob of hand hygiene after a contact]
  c<-par[27] # contact rate: number of contacts a patient on average gets per day 
  time.int<-par[28] # time of intervention
  new.hhfreq<-par[29] # new hand hygiene frequency after the intervention
  f23 <- par[30] # fraction of population 2 that mixes with 3 
  f32 <- par[31] # fraction of population 3 that mixes with 2 

    if(t>=time.int) hhfreq <- new.hhfreq
  
  N1<-S1+A1+B1
  N2<-S2+A2+B2
  N3<-S3+A3+B3
  
  # derive hand hygiene rate from hand hygiene frequency
  hh.rate <- (hhfreq*c*N1/n.hcw )/(1-hhfreq)
  
  hosp.discharge.rate<-mS1*S1 + mA1*A1 + mB1*B1
  pop3.hospitalisation.rate<-hosp.discharge.rate/(N3+rho*N2)
  hS3<-pop3.hospitalisation.rate
  hA3<-pop3.hospitalisation.rate
  hB3<-pop3.hospitalisation.rate
  hS2<-rho*pop3.hospitalisation.rate
  hA2<-rho*pop3.hospitalisation.rate
  hB2<-rho*pop3.hospitalisation.rate
  
  # the following holds if the population sizes are assumed to be stable
  mS2<- N3*pop3.hospitalisation.rate/N2 # rate those in S2 move to S3
  mA2<-mS2 #mA2<-par[29] #assumed to be the same as mS2 
  mB2<-mS2 #mB2<-par[30] #assumed to be the same as mS2   
  
  # Force of infection 
  foi_A2 <- betaA2*A2/N2 + betaA3*A3*f23/N2 # We assume an individual in population 2 is more likely to mix with an individual in population 2 than an individual in population 3
  foi_A3 <- betaA3*A3/N3 + betaA2*A2*f32/N3 # We assume an individual in population 3 is equally likely to mix with an individual in population 3 as an individual in population 2
  foi_B2 <- betaB2*B2/N2 + betaB3*B3*f23/N2
  foi_B3 <- betaB3*B3/N3 + betaB2*B2*f32/N3
  
  # Switch from A --> B and B --> A
  A2toB2 <- betaBA2*B2/N2 + betaBA3*B3*f23/N2 
  A3toB3 <- betaBA3*B3/N3 + betaBA2*B2*f32/N3 
  B2toA2 <- betaAB2*A2/N2 + betaAB3*A3*f23/N2 
  B3toA3 <- betaAB3*A3/N3 + betaAB2*A2*f32/N3 
  
  # Derivatives 
  dS1 <-   hS2*S2 + hS3*S3 - mS1*S1 - betaA1.hcw2p*S1*H.A1/n.hcw - betaB1.hcw2p*S1*H.B1/n.hcw + muA1*A1 + muB1*B1
  dS2 <- - hS2*S2 - mS2*S2 + mS1*S1 - foi_A2*S2 - foi_B2*S2 + muA2*A2 + muB2*B2
  dS3 <- - hS3*S3 + mS2*S2 - foi_A3*S3 - foi_B3*S3  + muA3*A3 + muB3*B3
  dA1 <-   hA2*A2 + hA3*A3 - mA1*A1 - muA1*A1 + betaA1.hcw2p*S1*H.A1/n.hcw + betaAB1.hcw2p*H.A1*B1/n.hcw - betaBA1.hcw2p*H.B1*A1/n.hcw 
  dA2 <- - hA2*A2 - mA2*A2 + mA1*A1 - muA2*A2 + foi_A2*S2 + B2toA2*B2 - A2toB2*A2  # Why would these recently hospitalised individuals not mix with the rest of the community?
  dA3 <- - hA3*A3 + mA2*A2 - muA3*A3 + foi_A3*S3 + B3toA3*B3 - A3toB3*A3
  dB1 <-   hB2*B2 + hB3*B3 - mB1*B1 - muB1*B1+ betaB1.hcw2p*S1*H.B1/n.hcw + betaBA1.hcw2p*H.B1*A1/n.hcw - betaAB1.hcw2p*H.A1*B1/n.hcw # when changing colonisation rate (reduce it), this one becomes negative, why?
  dB2 <- - hB2*B2 - mB2*B2 + mB1*B1 - muB2*B2 + foi_B2*S2 + A2toB2*A2 - B2toA2*B2
  dB3 <- - hB3*B3 + mB2*B2 - muB3*B3 + foi_B3*S3 + A3toB3*A3 - B3toA3*B3 
  dH.A1 <- betaA1.p2hcw*A1*(n.hcw - H.A1 -H.B1)/n.hcw - hh.rate*H.A1
  dH.B1 <- betaB1.p2hcw*B1*(n.hcw - H.A1 -H.B1)/n.hcw - hh.rate*H.B1
  dTotal <- dS1+dS2+dS3+dA1+dA2+dA3+dB1+dB2+dB3 
  
  Inc_A1 <- betaA1.hcw2p*S1*H.A1/n.hcw + betaAB1.hcw2p*H.A1*B1/n.hcw
  Inc_A2 <- foi_A2*S2 + B2toA2*B2  
  Inc_A3 <- foi_A3*S3 + B3toA3*B3
  Inc_B1 <- betaB1.hcw2p*S1*H.B1/n.hcw + betaBA1.hcw2p*A1*H.B1/n.hcw
  Inc_B2 <- foi_B2*S2 + A2toB2*A2
  Inc_B3 <- foi_B3*S3 + A3toB3*A3
  
  # Return values 
  list(c(dS1,dA1,dB1,dS2,dA2,dB2,dS3,dA3,dB3,dH.A1,dH.B1,dTotal),
       Inc_A1/N1*10000,Inc_A2/N2*10000,Inc_A3/N3*10000,Inc_B1/N1*10000,Inc_B2/N2*10000,Inc_B3/N3*10000,
       Inc_A1,Inc_A2,Inc_A3,Inc_B1,Inc_B2,Inc_B3)
}

calcR0.fqr<-function(S1.t0,A1.t0, B1.t0,S2.t0,A2.t0, B2.t0,S3.t0,A3.t0, B3.t0,mS1,mA1,mB1,rho,muA1,muA2,muA3,betaA1.hcw2p, betaA1.p2hcw,betaA2,betaA3,n.hcw, hhfreq, c,f23,f32){
  # calcs below work out R0 and next generation matrix 
  N1<-sum(S1.t0+A1.t0+B1.t0)
  N2<-sum(S2.t0+A2.t0+B2.t0)
  N3<-sum(S3.t0+A3.t0+B3.t0)
  
  # derive hand hygiene rate from hand hygiene frequency
  hh.rate <- (hhfreq*c*N1/n.hcw )/(1-hhfreq)
  
  hosp.discharge.rate<-mS1*S1.t0 + mA1*A1.t0 + mB1*B1.t0
  pop3.hospitalisation.rate<-hosp.discharge.rate/(N3+rho*N2)
  hS3<-pop3.hospitalisation.rate
  hA3<-pop3.hospitalisation.rate
  hB3<-pop3.hospitalisation.rate
  hS2<-rho*pop3.hospitalisation.rate
  hA2<-rho*pop3.hospitalisation.rate
  hB2<-rho*pop3.hospitalisation.rate
  
  # the following holds if the population sizes are assumed to be stable
  mS2<- N3*pop3.hospitalisation.rate/N2 # rate those in S2 move to S3
  mA2<-mS2 #mA2<-par[29] #assumed to be the same as mS2
  mB2<-mS2 #mB2<-par[30] #assumed to be the same as mS2   
  
  p12A <-  mA1/(mA1 + muA1) # prob person with A in pop 1 moves to pop 2 still with A
  p21A <-  hA2/(hA2+ muA2 + mA2) # prob person with A in pop 2 moves to pop 1 still with A
  p23A <- mA2/(hA2+ muA2 + mA2)    # etc
  p31A <- hA3/(hA3 + muA3)
  R1As <-  (N1-1)*betaA1.p2hcw*betaA1.hcw2p/(n.hcw*hh.rate*(muA1+mA1)) #betaA1/(muA1 + mA1) # single episode repro number for A in pop 1
  R2As <- betaA2/(muA2 + mA2 + hA2) # single episode repro number for A in pop 2
  R3As <- betaA3/(muA3 + hA3) # single episode repro number for A in pop 3
  R23As <- betaA2*f32/(muA2 + mA2 + hA2) # including transmission from A in pop2 to pop3
  R32As <- betaA3*f23/(muA3 + hA3) # including transmission from A in pop3 to pop2
  
  R0A1<- (R1As + p12A*(R2As+R23As) + p12A*p23A*(R3As +R32As))/ (1-(p12A*p21A +p12A*p23A*p31A)) # R0 starting in pop1
  R0A2<- R2As + R23As + p21A*R0A1 + p23A*(R3As +R32As) + p23A*p31A*R0A1  # R0 starting in pop2
  R0A3<- R3As + R32As + p31A*R0A1  # R0 starting in pop3
  
  R0A1.1 <-(R1As )/ (1-(p12A*p21A +p12A*p23A*p31A)) #component of R0A1 due to transmission in pop 1
  R0A1.2<-(p12A*R2As + p12A*p23A*R32As)/ (1-(p12A*p21A +p12A*p23A*p31A)) #component of R0A1 due to transmission in pop 2
  R0A1.3<- (p12A*p23A*R3As +p12A*R23As)/ (1-(p12A*p21A +p12A*p23A*p31A)) # R0 starting in pop 3
  
  R11A<- R0A1.1 # RijA gives expected number of secondary cases in pop j resulting from one initial case in pop i (aacounting for readmissions)
  R12A<- R0A1.2
  R13A<- R0A1.3
  R21A <- p21A*R0A1.1 + p23A*p31A*R0A1.1
  R22A <- R2As + p23A*R32As + p21A*R0A1.2 + p23A*p31A*R0A1.2
  R23A <- p23A*R3As + R23As + p21A*R0A1.3 + p23A*p31A*R0A1.3  
  R31A <- p31A*R0A1.1 
  R32A <- R32As + p31A*R0A1.2  # including transmission from A in pop3 to pop2
  R33A <- R3As + p31A*R0A1.3
  
  ngm<-matrix(c(R11A,R12A,R13A, R21A,R22A,R23A,R31A,R32A,R33A), nrow=3, byrow=F)
  R0<-max(as.double(eigen(ngm)$values))
  return(list(R0=R0, ngm=ngm))
}

calcR0B.fqr<-function(S1.t0,A1.t0, B1.t0,S2.t0,A2.t0, B2.t0,S3.t0,A3.t0, B3.t0,mS1,mA1,mB1,rho,muB1,muB2,muB3,betaB1.hcw2p, betaB1.p2hcw,betaB2,betaB3,n.hcw, hhfreq, c,f23,f32){
  # calcs below work out R0 and next generation matrix 
  N1<-sum(S1.t0+A1.t0+B1.t0)
  N2<-sum(S2.t0+A2.t0+B2.t0)
  N3<-sum(S3.t0+A3.t0+B3.t0)
  
  # derive hand hygiene rate from hand hygiene frequency
  hh.rate <- (hhfreq*c*N1/n.hcw )/(1-hhfreq)
  
  hosp.discharge.rate<-mS1*S1.t0 + mA1*A1.t0 + mB1*B1.t0
  pop3.hospitalisation.rate<-hosp.discharge.rate/(N3+rho*N2)
  hS3<-pop3.hospitalisation.rate
  hA3<-pop3.hospitalisation.rate
  hB3<-pop3.hospitalisation.rate
  hS2<-rho*pop3.hospitalisation.rate
  hA2<-rho*pop3.hospitalisation.rate
  hB2<-rho*pop3.hospitalisation.rate
  
  # the following holds if the population sizes are assumed to be stable
  mS2<- N3*pop3.hospitalisation.rate/N2 # rate those in S2 move to S3
  mA2<-mS2 #mA2<-par[29] #assumed to be the same as mS2
  mB2<-mS2 #mB2<-par[30] #assumed to be the same as mS2   
  
  p12B <-  mB1/(mB1 + muB1) # prob person with B in pop 1 moves to pop 2 still with B
  p21B <-  hB2/(hB2+ muB2 + mB2) # prob person with B in pop 2 moves to pop 1 still with B
  p23B <- mB2/(hB2+ muB2 + mB2)    # etc
  p31B <- hB3/(hB3 + muB3)
  R1Bs <-  (N1-1)*betaB1.p2hcw*betaB1.hcw2p/(n.hcw*hh.rate*(muB1+mB1)) #single episode repro number for B in pop 1
  R2Bs <- betaB2/(muB2 + mB2 + hB2) # single episode repro number for B in pop 2
  R3Bs <- betaB3/(muB3 + hB3) # single episode repro number for B in pop 3
  R23Bs <- betaB2*f32/(muB2 + mB2 + hB2) # including transmission from B in pop2 to pop3
  R32Bs <- betaB3*f23/(muB3 + hB3) # including transmission from B in pop3 to pop2
  
  R0B1<- (R1Bs + p12B*(R2Bs+R23Bs) + p12B*p23B*(R3Bs +R32Bs))/ (1-(p12B*p21B +p12B*p23B*p31B)) # R0 starting in pop1
  R0B2<- R2Bs + R23Bs + p21B*R0B1 + p23B*(R3Bs +R32Bs) + p23B*p31B*R0B1  # R0 starting in pop2
  R0B3<- R3Bs + R32Bs + p31B*R0B1  # R0 starting in pop3
  
  R0B1.1 <-(R1Bs )/ (1-(p12B*p21B +p12B*p23B*p31B)) #component of R0B1 due to transmission in pop 1
  R0B1.2<-(p12B*R2Bs + p12B*p23B*R32Bs)/ (1-(p12B*p21B +p12B*p23B*p31B)) #component of R0B1 due to transmission in pop 2
  R0B1.3<- (p12B*p23B*R3Bs +p12B*R23Bs)/ (1-(p12B*p21B +p12B*p23B*p31B)) # R0 starting in pop 3
  
  R11B<- R0B1.1 # RijB gives expected number of secondary cases in pop j resulting from one initial case in pop i (aacounting for readmissions)
  R12B<- R0B1.2
  R13B<- R0B1.3
  R21B <- p21B*R0B1.1 + p23B*p31B*R0B1.1
  R22B <- R2Bs + p23B*R32Bs + p21B*R0B1.2 + p23B*p31B*R0B1.2
  R23B <- p23B*R3Bs + R23Bs + p21B*R0B1.3 + p23B*p31B*R0B1.3  
  R31B <- p31B*R0B1.1 
  R32B <- R32Bs + p31B*R0B1.2  # HERE - including transmission from B in pop3 to pop2
  R33B <- R3Bs + p31B*R0B1.3
  
  
  ngm<-matrix(c(R11B,R12B,R13B, R21B,R22B,R23B,R31B,R32B,R33B), nrow=3, byrow=F)
  R0<-max(as.double(eigen(ngm)$values))
  return(list(R0=R0, ngm=ngm))
}

