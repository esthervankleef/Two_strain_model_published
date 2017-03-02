#############################################
# FIGURES CDI model - Lancet ID correspondence
#############################################


# Author: E. van Kleef
# Date: 7 February 2017

rm(list=ls())
load("./Output/hosp_com_cdi.Rda")

library(ggplot2)

###################################
# PLOT CORRESPONDENCE LANCET ID
col1<-grey(0.3); lty1<-3
col2<-grey(0.8); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.3); lty4<-1

rqpois = function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

time.int <- 365000
Iint <- time.int/10
Ipost_int <- Iint+365

period = 365

# Incidence is per day so *365 = Incidence per year
random.seed=60
n=1000
size = 15
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


irr.data = data.frame(var=c(1:8), irr.median=c(median(irrA.h,na.rm=T),0.21,median(irrA.c,na.rm=T),0.45,median(irrB.h,na.rm=T),0.87,median(irrB.c,na.rm=T),1.14),
                      q1 = c(quantile(irrA.h,0.05,na.rm=T),0.13,quantile(irrA.c,0.05,na.rm=T),0.29,quantile(irrB.h,0.05,na.rm=T),0.67,quantile(irrB.c,0.05,na.rm=T),0.92),
                      q3 = c(quantile(irrA.h,0.95,na.rm=T),0.34,quantile(irrA.c,0.95,na.rm=T),0.71,quantile(irrB.h,0.95,na.rm=T),1.13,quantile(irrB.c,0.95,na.rm=T),1.42),
                      type=c("With hospital link A - model","With hospital link A - data", "No hospital link A - model","No hospital link A - data","With hospital link B - model","With hospital link B - data","No hospital link B - model","No hospital link B - data"),
                      path=c("A - m","A - d","A - m","A - d","B - m","B - d","B - m","B - d"),
                      path2=c(rep("A",4), rep("B",4)), model = c(1,0,1,0,1,0,1,0))

irr.data$type=factor(irr.data$type,levels=c("With hospital link A - model","With hospital link A - data","With hospital link B - model","With hospital link B - data",
                                            "No hospital link A - model","No hospital link A - data","No hospital link B - model", "No hospital link B - data"))
irr.data

irr.data.long = data.frame(var=c(rep(1,n),rep(NA,n),rep(2,n),rep(NA,n),rep(3,n),rep(NA,n),rep(4,n),rep(NA,n)), irr=c(irrA.h,rep(NA,n),irrA.c,rep(NA,n),irrB.h,rep(NA,n),irrB.c,rep(NA,n)),
                           type=c(rep("With hospital link A",n),rep("With hospital link FQ",n),rep("No hospital link A",n),rep("No hospital link FQ",n),
                                  rep("With hospital link B",n),rep("With hospital link FS",n),rep("No hospital link B",n),rep("No hospital link FS",n)),
                           path=c(rep("A",n),rep(NA,n),rep("A",n),rep(NA,n),rep("B",n),rep(NA,n),rep("B",n),rep(NA,n)),
                           cdi = c(rep(NA,n),rep(0.21,n),rep(NA,n),rep(0.45,n),rep(NA,n),rep(0.87,n),rep(NA,n), rep(1.14,n)))
irr.data.long$type=factor(irr.data.long$type,levels=c("With hospital link A","With hospital link FQ","With hospital link B","With hospital link FS",
                                                      "No hospital link A","No hospital link FQ","No hospital link B","No hospital link FS"))
head(irr.data.long)
tail(irr.data.long)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


print(ggplot(irr.data[!irr.data$type%in%c("All A","All B"),],aes(x=type,y=irr.median,col=path2,lty=factor(model),shape=factor(model)))+
        scale_y_log10(breaks = c(0,0.25,0.5,1,2,4),labels = c(0,0.25,0.5,1,2,4))+
        geom_point(size=4)+geom_hline(yintercept=1,lty=2)+
        theme_bw()+geom_errorbar(aes(ymin=q1, ymax=q3),width=0,size=0.8)+ylab("IRR per year")+geom_vline(xintercept=4.5,col="grey",alpha=0.3)+
        theme(axis.title.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                                                                            axis.ticks.x=element_blank(),
                                                                                            axis.text.x = element_text(hjust = 0.01),
                                                                                            axis.text=element_text(size=10),
                                                                                            legend.position="bottom")+scale_x_discrete("", labels = c("","Hospital link","","",
                                                                                                                                                      "","No hospital link","","","",""))+
        scale_color_manual("", values=c("darkred","darkblue"),labels=c("Change in fluoroquinolone-resistant secondary cases",
                                                                       "Change in fluoroquinolone-sensitive secondary cases")) +
        scale_linetype_manual("",values=c(1,12))+
        scale_shape_manual("",values=c(20,18,20,18),labels=c("Dingle et al (2017)","Model output"))+
        guides(color=guide_legend(nrow=2,byrow=TRUE),
               shape=guide_legend(nrow=2,byrow=TRUE),linetype=F))

#################################
# PLOT PREVALENCE
par(mfrow=c(2,2),oma=c(3,1,4,1),mar=c(1,4,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-0.4
col1<-grey(0.2); lty1<-3
col2<-grey(0.4); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.8); lty4<-1
lw1<-3 #default line weight

plot(matA1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Hospital prevalence (%)",xlab="Time (years)",main="")
lines(matA1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matA1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matA1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax,0.1),labels=seq(0,ymax*100,10))
legendtext<-c("5% compliance gain","7.5% compliance gain" ,
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
legendtext<-c("5% compliance gain","7.5% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

par(mar=c(4,4,1,1))
plot(matA2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Community prevalence (%)",xlab="Time (years)",main="")
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
mtext("Resistant pathogen",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive pathogen",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)

# PLOT INCIDENCE
par(mfrow=c(2,2),oma=c(3,1,4,1),mar=c(1,4,1,1))

time.int = 365000
xmin<-time.int/10-365
xmax<-xmin+365*5
ymax<-120
col1<-grey(0.2); lty1<-3
col2<-grey(0.4); lty2<-1
col3<-grey(0.6); lty3<-1
col4<-grey(0.8); lty4<-1
lw1<-3 #default line weight

plot(matIA1[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Hospital incidence (x10 000 person-days)",xlab="Time (years)",main="")
lines(matIA1[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIA1[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIA1[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,25),labels=seq(0,ymax+10,25))
legendtext<-c("5% compliance gain","7.5% compliance gain" ,
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
legendtext<-c("5% compliance gain","7.5% compliance gain" ,
              "15% compliance gain", "20% compliance gain") 
legend("topleft",legend=legendtext, col=c(col4,col3, col2,col1),lty=c(1,lty2,lty3,lty4),lwd=lw1,bty='n')

ymax<-50
par(mar=c(4,4,1,1))
plot(matIA2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Community incidence (x10 000 person-days)",xlab="Time (years)",main="")
lines(matIA2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIA2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIA2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))

par(mar=c(4,1,1,4))
plot(matIB2and3[xmin:xmax,6],type='l',ylim=c(0,ymax),xaxt='n',yaxt='n',lty=1,lwd=lw1, col=col4,bty='n',ylab="Incidence (x10 000 person-days)",xlab="Time (years)",main="")
lines(matIB2and3[xmin:xmax,7], lty=lty2, lwd=lw1,col=col3)
lines(matIB2and3[xmin:xmax,8], lty=lty3, lwd=lw1,col=col2)
lines(matIB2and3[xmin:xmax,9], lty=lty4,lwd=lw1,col=col1)
days<-xmax-xmin
years<-round(days/365)
axis(side=1,at=(0:years)*365,labels=0:years)
axis(side=2,at=seq(0,ymax+10,10),labels=seq(0,ymax+10,10))
mtext("Resistant pathogen",side=3,at=.3,line=2,outer=TRUE)
mtext("Adapted to hospital",side=3,at=.3,line=1,outer=TRUE,cex=.75)
mtext("Sensitive pathogen",side=3,at=.75,line=2,outer=TRUE)
mtext("Adapted to community",side=3,at=.75,line=1,outer=TRUE,cex=.75)

