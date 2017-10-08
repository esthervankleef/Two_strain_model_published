##########################################
# FIGURES COMMUNITY VS HOSPITAL MODEL
##########################################

# Author: E. van Kleef
# Date: 7 February 2017

rm(list=ls())
load("./Output/hosp_com.Rda")
load("./Output/hosp_com_sens.Rda")
load("./Output/hosp_com.sens_bi_mix.Rda")

source("./R/im_scale.R")
library(ggplot2)

#######################################################
# FIGURES VanKleef et al (2017) Wellcome Open Research
#######################################################
col1<-grey(0.3); lty1<-3
col2<-grey(0.8); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.3); lty4<-1
lw1<-3 #default line weight

# FIGURE 2 - Reduction IRR at baseline
rqpois = function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

time.int <- 365000
Iint <- time.int/10
Ipost_int <- Iint+365

period = 365

random.seed=60
n=1000
size = 5
std <- function(x) sd(x)/sqrt(length(x))
prob=0.1

pirrA.h.post = rbinom(n,rqpois(n,mu=round(sum(matdIA1[Iint:Ipost_int,7])),theta=size),prob=prob)
pirrA.h.pre = rbinom(n,rqpois(n,mu=round(sum(matdIA1[(Iint-365):Iint,7])),theta=size),prob=prob)
irrA.h = pirrA.h.post/pirrA.h.pre

pirrA.c.post = rbinom(n,rqpois(n,mu=round(sum(matdIA2and3[Iint:Ipost_int,7])),theta=size),prob=prob)
pirrA.c.pre = rbinom(n,rqpois(n,mu=round(sum(matdIA2and3[(Iint-365):Iint,7])),theta=size),prob=prob)
irrA.c = pirrA.c.post/pirrA.c.pre
irrA.all = (pirrA.h.post+pirrA.c.post)/(pirrA.h.pre+pirrA.c.pre)

pirrB.h.post = rbinom(n,rqpois(n,mu=round(sum(matdIB1[Iint:Ipost_int,7])),theta =size),prob=prob)
pirrB.h.pre = rbinom(n,rqpois(n,mu=round(sum(matdIB1[(Iint-365):Iint,7])),theta=size),prob=prob)
irrB.h = pirrB.h.post/pirrB.h.pre

pirrB.c.post = rbinom(n,rqpois(n,mu=round(sum(matdIB2and3[Iint:Ipost_int,7])),theta=size),prob=prob)
pirrB.c.pre = rbinom(n,rqpois(n,mu=round(sum(matdIB2and3[(Iint-365):Iint,7])),theta=size),prob=prob)
irrB.c = pirrB.c.post/pirrB.c.pre
irrB.all = (pirrB.h.post+pirrB.c.post)/(pirrB.h.pre+pirrB.c.pre)


irr.data = data.frame(var=c(1:6), irr.median=c(median(irrA.h,na.rm=T),median(irrA.c,na.rm=T),median(irrA.all,na.rm=T),median(irrB.h,na.rm=T),median(irrB.c,na.rm=T),median(irrB.all,na.rm=T)),
                      q1 = c(quantile(irrA.h,0.05,na.rm=T),quantile(irrA.c,0.05,na.rm=T),quantile(irrA.all,0.05,na.rm=T),quantile(irrB.h,0.05,na.rm=T),quantile(irrB.c,0.05,na.rm=T),quantile(irrB.all,0.05,na.rm=T)),
                      q3 = c(quantile(irrA.h,0.95,na.rm=T),quantile(irrA.c,0.95,na.rm=T),quantile(irrA.all,0.95,na.rm=T),quantile(irrB.h,0.95,na.rm=T),quantile(irrB.c,0.95,na.rm=T),quantile(irrB.all,0.95,na.rm=T)),
                      type=c("With hospital link A","No hospital link A","All A","With hospital link B","No hospital link B","All B"),
                      path=c("A","A","A","B","B","B"))

irr.data$type=factor(irr.data$type,levels=c("With hospital link A","With hospital link B","No hospital link A","No hospital link B", "All A", "All B"))
irr.data

irr.data.long = data.frame(var=c(rep(1,n),rep(NA,n),rep(2,n),rep(NA,n),rep(3,n),rep(NA,n),rep(4,n),rep(NA,n),rep(5,n),rep(NA,n),rep(6,n),rep(NA,n),rep(NA,n)), irr=c(irrA.h,rep(NA,n),irrA.c,rep(NA,n),irrB.h,rep(NA,n),irrB.c,rep(NA,n),
                                                                                                                                                                     irrA.all,rep(NA,n),irrB.all,rep(NA,n),rep(NA,n)),
                           type=c(rep("With hospital link A",n),rep("With hospital link FQ",n),rep("No hospital link A",n),rep("No hospital link FQ",n),
                                  rep("With hospital link B",n),rep("With hospital link FS",n),rep("No hospital link B",n),rep("No hospital link FS",n),
                                  rep("Total A",n),rep("Total FQ CDI",n),rep("Total B",n),rep("Total FS CDI",n),rep(NA,n)),
                           path1=c(rep("A",n),rep(NA,n),rep("A",n),rep(NA,n),rep("B",n),rep(NA,n),rep("B",n),rep(NA,n),rep("A",n),rep(NA,n),rep("B",n),rep(NA,n),rep(NA,n)),
                           path2=c(rep(NA,n),rep("A",n),rep(NA,n),rep("A",n),rep(NA,n),rep("B",n),rep(NA,n),rep("B",n),rep(NA,n),rep("A",n),rep(NA,n),rep("B",n),rep(NA,n)),
                           cdi = c(rep(NA,n),rep(0.21,n),rep(NA,n),rep(0.45,n),rep(NA,n),rep(0.87,n),rep(NA,n), rep(1.14,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n)),
                           cdiq1 = c(rep(NA,n),rep(0.13,n),rep(NA,n),rep(0.29,n),rep(NA,n),rep(0.67,n),rep(NA,n), rep(0.92,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n)),
                           cdiq3 = c(rep(NA,n),rep(0.34,n),rep(NA,n),rep(0.71,n),rep(NA,n),rep(1.13,n),rep(NA,n), rep(1.42,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n),rep(NA,n)))

irr.data.long$type=factor(irr.data.long$type,levels=c("With hospital link FQ","With hospital link A","With hospital link FS","With hospital link B",
                                                      "No hospital link FQ","No hospital link A","No hospital link FS","No hospital link B",
                                                      "Total FQ CDI","Total A","Total FS CDI","Total B"))
head(irr.data.long)
tail(irr.data.long)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

print(ggplot(irr.data.long,aes(x=type,y=irr,fill=path1,group=var))+geom_vline(xintercept=4.7,col="grey",alpha=0.2)+geom_vline(xintercept=8.9,col="grey",alpha=0.2)+
        geom_violin()+geom_hline(yintercept=1,lty=2)+
        theme_bw()+ylab("IRR (per year)")+stat_summary(fun.data=data_summary,col="black",shape=16,size=0.6)+
        scale_y_continuous(limits=c(0, 1.6),breaks=seq(0,1.6,0.2),labels=seq(0,1.6,0.2))+theme(axis.title.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                                                                               axis.ticks.x=element_blank(),
                                                                                               axis.text.x = element_text(hjust = 0.01),
                                                                                               axis.text=element_text(size=10),
                                                                                               legend.position="bottom")+scale_x_discrete("", labels = c("","Hospital-linked","","","",
                                                                                                                                                         "Community-linked","","","","",
                                                                                                                                                         "Overall","","","",""))+
        scale_fill_manual("", values=c("darkslategray3",col3),labels=c("Infections with resistant strain (simulations)",
                                                                       "Infections with sensitive strain (simulations)")) +
        guides(fill=guide_legend(nrow=4,byrow=TRUE),col=guide_legend(nrow=4,byrow=TRUE),shape=F)+
        geom_point(aes(x=type,y=cdi,col=path2,shape=path2),size=2.7)+
        geom_segment(aes(x = 1, y = 0.13, xend = 1, yend = 0.34),col="darkslategray3")+
        geom_segment(aes(x = 5, y = 0.29, xend = 5, yend = 0.71),col="darkslategray3")+
        geom_segment(aes(x = 3, y = 0.67, xend = 3, yend = 1.13),col=col3)+
        geom_segment(aes(x = 7, y = 0.92, xend = 7, yend = 1.42),col=col3)+
        scale_colour_manual("",values=c("darkslategray3",col3,"white"),
                            labels=c("Infections with resistant strain (data from [4])",
                                     "Infections with sensitive strain (data from [4])",
                                     ""))+
        scale_shape_manual("",values=c(16,16),
                           labels=c("Infections with resistant strain (data from [4])",
                                    "Infections with sensitive strain (data from [4])",
                                    "")))


# FIGURE 3 - incidence and prevalence at baseline
par(mfrow=c(2,2),oma=c(3,1,4,1),mar=c(1,5,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-100

plot(matIA1[xmin:xmax,7],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=3,lwd=lw1, col=col3,bty='n',ylab="Incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
#lines(c(time.int/10-xmin,time.int+10/10-xmin),rep(ymax-10,length(c(time.int/10-xmin,time.int+10/10-xmin)),lty=1,lwd=lw1,col="red")
lines(matIA1[xmin:xmax,7],lty=3, lwd=lw1,col=col3)
lines(matIA2and3[xmin:xmax,7],lty=2, lwd=lw1,col=col2)
lines(matIA1[xmin:xmax,7]+matIA2and3[xmin:xmax,7],lty=1, lwd=lw1,col=col1)
all.eq = round((matIA1[xmin:xmax,7]+matIA2and3[xmin:xmax,7])[length(matIA1[xmin:xmax,7]+matIA2and3[xmin:xmax,7])],1)
all.s = round((matIA1[xmin:xmax,7]+matIA2and3[xmin:xmax,7])[1],1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,20),labels=seq(0,ymax+10,20))
legendtext<-c("All", "Hospital","Community") 
legend("topleft",legend=legendtext, col=c(col1,col3,col2),lty=c(1,3,2),lwd=lw1,bty='n',cex=0.8)

par(mar=c(1,1,1,4))
plot(matIB1[xmin:xmax,7],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=3,lwd=lw1, col=col3,bty='n', xlab="Time (years)", ylab="Incidence (acquisitions/10,000 person days)",main="")
lines(matIB1[xmin:xmax,7],lty=3, lwd=lw1,col=col3)
lines(matIB2and3[xmin:xmax,7],lty=2, lwd=lw1,col=col2)
lines(matIB1[xmin:xmax,7]+matIB2and3[xmin:xmax,7],lty=1, lwd=lw1,col=col1)
legendtext<-c("All", "Hospital","Community") 
legend("topleft",legend=legendtext, col=c(col1,col3,col2),lty=c(1,3,2),lwd=lw1,bty='n',cex=0.8)
all.eq = round((matIB1[xmin:xmax,7]+matIB2and3[xmin:xmax,7])[length(matIB1[xmin:xmax,7]+matIB2and3[xmin:xmax,7])],1)
all.s = round((matIB1[xmin:xmax,7]+matIB2and3[xmin:xmax,7])[1],1)
axis(side=2,at=seq(0,ymax+10,20),labels=seq(0,ymax+10,20))
axis(side=1,at=(0:years)*365,labels=0:years)

xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-0.5
col1<-grey(0.3); lty1<-3
col2<-grey(0.8); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.3); lty4<-1
lw1<-3 #default line weight

par(mar=c(1,5,1,1))
plot(matA1[xmin:xmax,7],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=3,lwd=lw1, col=col3,bty='n',ylab="Prevalence\n(% of patients colonized)",xlab="Time (years)",main="")
lines(matA1[xmin:xmax,7],lty=3, lwd=lw1,col=col3)
lines(matA2and3[xmin:xmax,7],lty=2, lwd=lw1,col=col2)
lines(matA.all[xmin:xmax,7],lty=1, lwd=lw1,col=col1)
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))

par(mar=c(1,1,1,4))
plot(matB1[xmin:xmax,7],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=3,lwd=lw1, col=col3,bty='n', xlab="Time (years)", ylab="Incidence (per 10 000 beddays)",main="")
lines(matB1[xmin:xmax,7],lty=3, lwd=lw1,col=col3)
lines(matB2and3[xmin:xmax,7],lty=2, lwd=lw1,col=col2)
lines(matB.all[xmin:xmax,7],lty=1, lwd=lw1,col=col1)
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))


mtext("Resistant strain",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive strain",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)
mtext("Time (years)",side=1,at=.29,line=1,outer=TRUE)
mtext("Time (years)",side=1,at=.72,line=1,outer=TRUE)

# FIGURE 4 - Varying degree of hospital and community adaptation
head(sens.out)
irr.sens = sens.out[,c(1:4)]

# If annual incidence is <1 before intervention (i.e. no coexistence), set incidence to 0
nocoex = ifelse(sens.out$preA1+sens.out$preA2and3<=1|sens.out$preB1+sens.out$preB2and3<=1,0,1)
irr.sens = sens.out[c("postA1","postA2and3","postB1", "postB2and3")]/sens.out[,c("preA1","preA2and3","preB1", "preB2and3")]
irr.sens = cbind(sens.out[,c(1:2)],irr.sens)
irr.sens[which(nocoex==0),c("postA1","postA2and3","postB1", "postB2and3")] = NA

# Calculate irr
irr.Ah = matrix(NA,nrow=length(unique(irr.sens$fhA)),ncol=length(unique(irr.sens$fhB))) # rows = fha, cols = fhb
irr = list(irr.Ah, irr.Ac=irr.Ah,irr.Bh=irr.Ah,irr.Bc=irr.Ah)

for(m in 1:length(irr)){
  for(i in unique(irr.sens$fhA)){
    f = which(unique(irr.sens$fhA)==i)
    rows = which(irr.sens$fhA==i)
    irr[[m]][f,] = round(irr.sens[rows,2+m],2)
  }
}

# create discrete levels (i.e. IRR <1. IRR = 1, IRR >1)
irr.dis = list(irr.Ah, irr.Ac=irr.Ah,irr.Bh=irr.Ah,irr.Bc=irr.Ah)
for(i in 1:length(irr.dis)){
  irr.dis[[i]] = ifelse(round(irr[[i]],1)<1,0,
                        ifelse(round(irr[[i]],1)==1,1,2)) 
}

breaks = seq(0,1.4,length=101)
pA = c(0.005,0.025,0.05*c(1:12))*100

pB = c(0.005,0.025,0.05*c(1:3))*100

parOld = par()
par(oma=c(1,2,5,3))
par(mar=c(5,6,3,1)) # Numbers refer to c(bottom, left, top, right)
graphics::layout(cbind(matrix(1:4,ncol=2),rep(5,2)), widths=c(2,2,0.3), heights=c(2,2))
#layout.show(5)

cols_neg <- colorRampPalette(c("darkblue","turquoise3", "yellow"))(n=length(breaks[breaks<1]))
cols_pos <-colorRampPalette(c("darkgoldenrod1", "red", "black"))(n = length(breaks[breaks>=1])-1)

cols1 = c(cols_neg,cols_pos)
cex1=1.5
cex2=1.2
image(x=pB, y=pA, z=t(irr[[1]][1:length(pA),1:length(pB)]),breaks=breaks,col=cols1, ylab ="% Resistant pathogen hospital-acquired" ,xlab = "% Sensitive pathogen hospital-acquired", cex.lab=1.5,cex.axis=1.5)
mtext("Resistant strain", side=3, line=4, adj=0.5,cex = cex1)
mtext("Hospital", side=3, line=0.5, adj=0.5,cex = cex2)

image(x=pB, y=pA, z=t(irr[[2]][1:length(pA),1:length(pB)]),breaks=breaks,col=cols1, ylab = "% Resistant pathogen hospital-acquired" ,xlab = "% Sensitive pathogen hospital-acquired", cex.lab=1.5,cex.axis=1.5)
mtext("Community", side=3, line=0.5, adj=0.5,cex =  cex2)

image(x=pB, y=pA, z=t(irr[[3]][1:length(pA),1:length(pB)]),breaks=breaks,col=cols1, ylab = "% Resistant pathogen hospital-acquired" ,xlab = "% Sensitive pathogen hospital-acquired", cex.lab=1.5,cex.axis=1.5)
mtext("Sensitive strain", side=3, line=4, adj=0.5,cex = cex1)
mtext("Hospital", side=3, line=0.5, adj=0.5,cex =  cex2)

image(x=pB, y=pA, z=t(irr[[4]][1:length(pA),1:length(pB)]),breaks=breaks,col=cols1,  ylab ="% Resistant pathogen hospital-acquired" ,xlab = "% Sensitive pathogen hospital-acquired", cex.lab=1.5,cex.axis=1.5)
mtext("Community", side=3, line=0.5, adj=0.5,cex =  cex2)

par(mar=c(5,0,1,2))

image.scale(c(0,0.1), breaks=breaks,col=cols1,horiz=FALSE, yaxt="n", xaxt="n", xlab="", ylab="",cex=cex2)
axis(4)
mtext("IRR (per year)", side=4, line=2.5,cex=cex2)
box()


# FIGURE S1 - Incidence under different levels of hand hygiene compliance
par(mfrow=c(2,2),oma=c(3,2,4,1),mar=c(1,5,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-100
col1<-grey(0.2); lty1<-3
col2<-grey(0.4); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.8); lty4<-1
lw1<-3 #default line weight

plot(matIA1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Hospital incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(matIA1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIA1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIA1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,25),labels=seq(0,ymax+10,25))
legendtext<-c("5% compliance gain","10% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

par(mar=c(1,1,1,4))
plot(matIB1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',xlab="Time (years)",main="")
lines(matIB1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIB1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIB1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,25),labels=seq(0,ymax+10,25))
legendtext<-c("5% compliance gain","10% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

ymax<-50
par(mar=c(4,5,1,1))
plot(matIA2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Community incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(matIA2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIA2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIA2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))

par(mar=c(4,1,1,4))
plot(matIB2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Community incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(matIB2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIB2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIB2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))
mtext("Resistant strain",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive strain",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)


# FIGURE S2 - Prevalence under different levels of hand hygiene compliance
par(mfrow=c(2,2),oma=c(3,1,4,1),mar=c(1,5,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-0.4
col1<-grey(0.2); lty1<-3
col2<-grey(0.4); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.8); lty4<-1
lw1<-3 #default line weight

plot(matA1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Hospital prevalence\n(% of patients colonized)",xlab="Time (years)",main="")
lines(matA1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matA1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matA1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))
legendtext<-c("5% compliance gain","10% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

par(mar=c(1,1,1,4))
plot(matB1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',xlab="Time (years)",main="")
lines(matB1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matB1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matB1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))
legendtext<-c("5% compliance gain","10% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

par(mar=c(4,5,1,1))
plot(matA2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Community prevalence\n(% of patients colonized)",xlab="Time (years)",main="")
lines(matA2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matA2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matA2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))

par(mar=c(4,1,1,4))
plot(matB2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Incidence (x10 000 person-days)",xlab="Time (years)",main="")
lines(matB2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matB2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matB2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))
mtext("Resistant strain",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive strain",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)

# FIGURE S3: change in bacterial inference and mixing
par(mfrow=c(2,2),oma=c(3,2,4,1),mar=c(1,5,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-200
col0<-grey(0.7); lty0<-3
col1<-"firebrick1"; lty1<-1
col2<-"lightcoral"; lty2<-1
col3<-"dodgerblue4"; lty3<-1
col4<-"dodgerblue"; lty4<-1
lw1<-3 #default line weight

plot(sensIA1[xmin:xmax,1],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=lty0,lwd=lw1, col=col0,bty='n',ylab="Hospital incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(sensIA1[xmin:xmax,2], lty=lty2, lwd=lw1,col=col4)
lines(sensIA1[xmin:xmax,3], lty=lty3, lwd=lw1,col=col3)
lines(sensIA1[xmin:xmax,4], lty=lty4,lwd=lw1,col=col2)
lines(sensIA1[xmin:xmax,5], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,25),labels=seq(0,ymax+10,25))
legendtext<-c("Baseline", "Bacterial inference (w = 0.5)","Bacterial inference (w = 1)" ,
              "Homogeneous mixing \n(f23 = N2/N3; f32 = N3/N2)", "No mixing (f23 = f32 = 0)") 
legend("topleft",legend=legendtext, col=c(col0,col4,col3, col2,col1),lty=c(lty0,lty1,lty2,lty3,lty4),lwd=lw1,bty='n')

par(mar=c(1,1,1,4))
plot(sensIB1[xmin:xmax,1],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=lty0,lwd=lw1, col=col0,bty='n',xlab="Time (years)",main="")
lines(sensIB1[xmin:xmax,2], lty=lty2, lwd=lw1,col=col4)
lines(sensIB1[xmin:xmax,3], lty=lty3, lwd=lw1,col=col3)
lines(sensIB1[xmin:xmax,4], lty=lty4,lwd=lw1,col=col2)
lines(sensIB1[xmin:xmax,5], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,25),labels=seq(0,ymax+10,25))
legendtext<-c("Baseline", "Bacterial inference (w = 0.5)","Bacterial inference (w = 1)" ,
              "Homogeneous mixing \n(f23 = N2/N3; f32 = N3/N2)", "No mixing (f23 = f32 = 0)") 
legend("topleft",legend=legendtext, col=c(col0,col4,col3,col2,col1),lty=c(lty0,lty1,lty2,lty3,lty4),lwd=lw1,bty='n')

ymax<-50
par(mar=c(4,5,1,1))
plot(sensIA2and3[xmin:xmax,1],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=lty0,lwd=lw1, col=col0,bty='n',ylab="Community incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(sensIA2and3[xmin:xmax,2], lty=lty2, lwd=lw1,col=col4)
lines(sensIA2and3[xmin:xmax,3], lty=lty3, lwd=lw1,col=col3)
lines(sensIA2and3[xmin:xmax,4], lty=lty4,lwd=lw1,col=col2)
lines(sensIA2and3[xmin:xmax,5], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))

par(mar=c(4,1,1,4))
plot(sensIB2and3[xmin:xmax,1],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=lty0,lwd=lw1, col=col0,bty='n',ylab="Community incidence \n(acquisitions/10,000 person days)",xlab="Time (years)",main="")
lines(sensIB2and3[xmin:xmax,2], lty=lty2, lwd=lw1,col=col4)
lines(sensIB2and3[xmin:xmax,3], lty=lty3, lwd=lw1,col=col3)
lines(sensIB2and3[xmin:xmax,4], lty=lty4,lwd=lw1,col=col2)
lines(sensIB2and3[xmin:xmax,5], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))
mtext("Resistant strain",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive strain",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)




