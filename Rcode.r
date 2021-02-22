

#Clear the work space
rm(list=ls())

#get and set wrokign directory
getwd()
setwd(getwd())


dat <- read.csv(file="InterLACE4sample.csv", header=TRUE, sep=",")

#R libraries
library(mfp)
library(ggplot2) #plots
library(ggpubr) #plots
library(lattice) #plots
library(cowplot) #plots
library(Hmisc) #data handling

d1 <- subset(dat, study==1); d2 <- subset(dat, study==2); d3 <- subset(dat, study==3); d4 <- subset(dat, study==4)
f1 <- mfp(meno ~ fp(bmi, df=4, select=NA, alpha=NA, scale=F) + factor(smoke)+factor(educ)+factor(child)+age, family=gaussian, data=d1)
summary(f1)
f1

f2 <- mfp(meno ~ fp(bmi, df=4, select=NA, alpha=NA, scale=F) + factor(smoke)+factor(educ)+factor(child)+age, family=gaussian, data=d2)
summary(f2)
f2

f3 <- mfp(meno ~ fp(bmi, df=4, select=NA, alpha=NA, scale=F) + factor(smoke)+factor(educ)+factor(child)+age, family=gaussian, data=d3)
summary(f3)
f3

f4 <- mfp(meno ~ fp(bmi, df=4, select=NA, alpha=NA, scale=F) + factor(smoke)+factor(educ)+factor(child)+age, family=gaussian, data=d4)
summary(f4)
f4

s1p <- predict(f1, se.fit=T, newdata=dat)
y1 <- s1p$fit; se1 <- s1p$se.fit; v1 <- se1^2
lcl1 <- y1-1.96*se1
ucl1 <- y1+1.96*se1

s2p <- predict(f2, se.fit=T, newdata=dat)
y2 <- s2p$fit; se2 <- s2p$se.fit; v2 <- se2^2
lcl2 <- y2-1.96*se2
ucl2 <- y2+1.96*se2

s3p <- predict(f3, se.fit=T, newdata=dat)
y3 <- s3p$fit; se3 <- s3p$se.fit; v3 <- se3^2
lcl3 <- y3-1.96*se3
ucl3 <- y3+1.96*se3

s4p <- predict(f4, se.fit=T, newdata=dat)
y4 <- s4p$fit; se4 <- s4p$se.fit; v4 <- se4^2
lcl4 <- y4-1.96*se4
ucl4 <- y4+1.96*se4

#Pooling the functional forms across studies
suminv = 1/v1 + 1/v2 + 1/v3 + 1/v4

#Standardised fixed effect weights
w1 = (1/v1)/suminv
w2 = (1/v2)/suminv
w3 = (1/v3)/suminv
w4 = (1/v4)/suminv

#Overall fixed effect estimate and the variance
phiFE = w1*y1 + w2*y2 + w3*y3 + w4*y4
varphiFE = 1/suminv
sephiFE = sqrt(varphiFE)
lclFE = phiFE - 1.96*sephiFE
uclFE = phiFE + 1.96*sephiFE

####Calculation random effect weights#########

Q <- ((y1 - phiFE)^2*(1/v1)) + ((y2 - phiFE)^2*(1/v2)) + ((y3 - phiFE)^2*(1/v3)) + ((y4 - phiFE)^2*(1/v4))
denom = suminv - ((1/v1)^2 + (1/v2)^2 + (1/v3)^2 + (1/v4)^2)/suminv


#S Squared
#tausq = max(0, ((Q - (4-1))/denom))
tausq1 = ((Q - (4-1))/denom)
tausq <- ifelse(tausq1<0,0,tausq1)

#random-effect weights
wran1 = 1/(v1 + tausq)
wran2 = 1/(v2 + tausq)
wran3 = 1/(v3 + tausq)
wran4 = 1/(v4 + tausq)
wransum = wran1 + wran2 + wran3 + wran4

w1std = wran1/wransum
w2std = wran2/wransum
w3std = wran3/wransum
w4std = wran4/wransum

#Overall random-effect estimate and the variance
phiRE =  w1std*y1 + w2std*y2 + w3std*y3 + w4std*y4
varphiRE = 1/wransum
sephiRE <- sqrt(varphiRE)
lclRE = phiRE - 1.96*sephiRE
uclRE = phiRE + 1.96*sephiRE

#Plot of the estimated function for each study 
  
  #Plot of the estimated function for each study - smoother applied - Include CI
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7))
  cord = coord_cartesian(xlim=c(10,50), ylim=c(35,65))
  
  d1 <- subset(dat, study==1); d3 <- subset(dat, study==3); d4 <- subset(dat, study==4)
  y1 <- subset(y1, dat$study==1); y2 <- subset(y2, dat$study==2); y3 <- subset(y3, dat$study==3); y4 <- subset(y4, dat$study==4)
  
  
  pr1plot <- ggplot(d1, aes(x=bmi, y=y1)) + geom_point(aes(y=meno), colour="lightgray", size=1.5) + 
    geom_smooth(aes(x=bmi, y=y1), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord
  
  pr2plot <- ggplot(d2, aes(x=bmi, y=y2)) + geom_point(aes(y=meno), colour="lightgray", size=1.5) + 
    geom_smooth(aes(x=bmi, y=y2), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord
  
  pr3plot <- ggplot(d3, aes(x=bmi, y=y3)) + geom_point(aes(y=meno), colour="lightgray", size=1.5) + 
    geom_smooth(aes(x=bmi, y=y3), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme + cord
  
  pr4plot <- ggplot(d4, aes(x=bmi, y=y4)) + geom_point(aes(y=meno), colour="lightgray", size=1.5) + 
    geom_smooth(aes(x=bmi, y=y4), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme + cord
  
  
  ggdraw(plot=NULL, xlim=c(0,4), ylim=c(0,4.1), clip="off") + 
    draw_plot(pr1plot, x=0, y=2, width=2, height=2) +
    draw_plot(pr2plot, x=2, y=2, width=2, height=2) +
    
    draw_plot(pr3plot, x=0, y=0, width=2, height=2) + 
    draw_plot(pr4plot, x=2, y=0, width=2, height=2) +
  
  draw_plot_label(label= c("Study 1","Study 2", "Study 3", "Study 4"), size=8, x=c(0.08,2.08,0.08,2.08), y=c(4.05,4.05,2.05,2.05))
    

  ggsave2("SmthEstFunStdy.jpg", width=15, height=20, units=c("cm"), dpi=350)
  
  

  #Plot of the estimated function for each study - smoother applied - Include CI
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7))
  cord = coord_cartesian(xlim=c(10,50), ylim=c(45,55))
  
  d1 <- subset(dat, study==1); d3 <- subset(dat, study==3); d4 <- subset(dat, study==4)
  y1 <- subset(y1, dat$study==1); y2 <- subset(y2, dat$study==2); y3 <- subset(y3, dat$study==3); y4 <- subset(y4, dat$study==4)
  lcl1 <- subset(lcl1, dat$study==1); lcl2 <- subset(lcl2, dat$study==2); lcl3 <- subset(lcl3, dat$study==3); lcl4 <- subset(lcl4, dat$study==4)
  ucl1 <- subset(ucl1, dat$study==1); ucl2 <- subset(ucl2, dat$study==2); ucl3 <- subset(ucl3, dat$study==3); ucl4 <- subset(ucl4, dat$study==4)
  
  
  pr1plot <- ggplot(d1, aes(x=bmi, y=y1)) + 
    geom_smooth(aes(x=bmi, y=y1), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lcl1), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=ucl1), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord
  
  pr2plot <- ggplot(d2, aes(x=bmi, y=y2)) + 
    geom_smooth(aes(x=bmi, y=y2), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lcl2), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=ucl2), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord
  
  pr3plot <- ggplot(d3, aes(x=bmi, y=y3)) + 
    geom_smooth(aes(x=bmi, y=y3), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lcl3), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=ucl3), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme + cord
  
  pr4plot <- ggplot(d4, aes(x=bmi, y=y4)) + 
    geom_smooth(aes(x=bmi, y=y4), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lcl4), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=ucl4), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme + cord

  
  ggdraw(plot=NULL, xlim=c(0,4), ylim=c(0,4.1), clip="off") + 
    draw_plot(pr1plot, x=0, y=2, width=2, height=2) +
    draw_plot(pr2plot, x=2, y=2, width=2, height=2) +
    
    draw_plot(pr3plot, x=0, y=0, width=2, height=2) + 
    draw_plot(pr4plot, x=2, y=0, width=2, height=2) +
    
    draw_plot_label(label= c("Study 1","Study 2", "Study 3", "Study 4"), size=8, x=c(0.08,2.08,0.08,2.08), y=c(4.05,4.05,2.05,2.05))
  
  
  ggsave2("SmthEstFunStdy_ci.jpg", width=15, height=20, units=c("cm"), dpi=350)
  
  
  
  #Plot the overall estimated function smoother applied  - fixed effect
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7))
  cord = coord_cartesian(xlim=c(10,50), ylim=c(45,55))
  
  fixprplot <- ggplot(dat, aes(x=bmi, y=phiFE)) + #geom_point(aes(y=meno), colour="lightgray", size=1.5) +
    geom_smooth(aes(x=bmi, y=phiFE), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lclFE), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=uclFE), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme + cord


  ggdraw(plot=NULL, xlim=c(0,2), ylim=c(0,2), clip="off") + 
    draw_plot(fixprplot, x=0, y=0, width=2, height=2) +
  
    
    ggsave2("FixedEstimate_ci.jpg", width=10, height=10, units=c("cm"), dpi=350)
  
  
  
  #Plot the overall estimated function smoother applied  - random effect
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7))
  cord = coord_cartesian(xlim=c(10,50), ylim=c(45,55))
  
  randprplot <- ggplot(dat, aes(x=bmi, y=phiRE)) + #geom_point(aes(y=meno), colour="lightgray", size=1.5) +
    geom_smooth(aes(x=bmi, y=phiRE), method="loess", se=F, size=0.8, colour="black") + 
    geom_smooth(aes(x=bmi, y=lclRE), method="loess", se=F, lty=5, size=0.8, colour="black") +
    geom_smooth(aes(x=bmi, y=uclRE), method="loess", se=F, lty=5, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord
  
  
  ggdraw(plot=NULL, xlim=c(0,2), ylim=c(0,2), clip="off") + 
    draw_plot(randprplot, x=0, y=0, width=2, height=2) +
    
    
    ggsave2("RandEstimate_ci.jpg", width=10, height=10, units=c("cm"), dpi=350)
  
  
  
  
  #Plot the overall estimated function smoother applied  - fixed effect weights
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7), 
    legend.position = "bottom")
    
  cord = coord_cartesian(xlim=c(10,50), ylim=c(0,1))
  
  wtprplot <- ggplot(dat, aes(x=bmi, y=suminv)) +
    geom_smooth(aes(x=bmi, y=w1, colour="Study 1"), method="loess", se=F, lty=1, size=0.8) + 
    geom_smooth(aes(x=bmi, y=w2, colour="Study 2"), method="loess", se=F, lty=2, size=0.8) +
    geom_smooth(aes(x=bmi, y=w3, colour="Study 3"), method="loess", se=F, lty=3, size=0.8) + 
    geom_smooth(aes(x=bmi, y=w4, colour="Study 4"), method="loess", se=F, lty=4, size=0.8) + labs(y="Fixed effect weights",x="BMI") + 
    scale_colour_manual(name="", values=c("red", "black", "blue", "green")) + mytheme + cord
  
  
  ggdraw(plot=NULL, xlim=c(0,2), ylim=c(0,2), clip="off") + 
    draw_plot(wtprplot, x=0, y=0, width=2, height=2) +
    
    
    ggsave2("FixedEffectWeights.jpg", width=10, height=10, units=c("cm"), dpi=350)
  
  
  
  
  #Plot the overall estimated function smoother applied  - random effect weights
  mytheme = theme(
    axis.title.x=element_text(size=7), 
    axis.title.y=element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7), 
    legend.position = "bottom")
  
  cord = coord_cartesian(xlim=c(10,50), ylim=c(0,1))
  
  wtprplot <- ggplot(dat, aes(x=bmi, y=wransum)) +
    geom_smooth(aes(x=bmi, y=w1std, colour="Study 1"), method="loess", se=F, lty=1, size=0.8) + 
    geom_smooth(aes(x=bmi, y=w2std, colour="Study 2"), method="loess", se=F, lty=2, size=0.8) +
    geom_smooth(aes(x=bmi, y=w3std, colour="Study 3"), method="loess", se=F, lty=3, size=0.8) + 
    geom_smooth(aes(x=bmi, y=w4std, colour="Study 4"), method="loess", se=F, lty=4, size=0.8) + labs(y="Random effect weights",x="BMI") + 
    scale_colour_manual(name="", values=c("red", "black", "blue", "green")) + mytheme + cord
  
  
  ggdraw(plot=NULL, xlim=c(0,2), ylim=c(0,2), clip="off") + 
    draw_plot(wtprplot, x=0, y=0, width=2, height=2) +
    
    
    ggsave2("RandomEffectWeights.jpg", width=10, height=10, units=c("cm"), dpi=350)
  
  

  