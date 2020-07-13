
###########################################################################
#Purpose: Pseudo-code for pooling study-specific non-linear relationships
#using fractional polynomial: a step-by-step guide.
#We use three studies in this code; study 1, 2, and 3
#R-version Rx64 3.5.2 or later*/
##########################################################################

#Clear the work space
rm(list=ls())

#set seed
set.seed(999)

#directory
getwd() #Get the directory
setwd("C:/xxx/yyy") #Change the working directory where you saved the data

#install the following R-packages
#ggplot2, ggpubr, lattice, cowplot, Hmisc

#library
library(ggplot2) #plots
library(ggpubr) #plots
library(lattice) #plots
library(cowplot) #plots
library(Hmisc) #data handling

#Read the original data for the paper
dat <- read.csv(file="hypodat.csv", header=TRUE, sep=",")
names(dat); dim(dat); head(dat)

#data format
dat$study <- factor(dat$study)
dat$smoking <- factor(dat$smoking)
dat$nochild <- factor(dat$nochild)
dat$education <- factor(dat$education)

#Step 1: Exploratory analysis---------------------------------------------------------
#Data description
summary(dat$bmi)
summary(dat$menopauseage)
summary(dat$age)
table(dat$smoking)
table(dat$nochild)
table(dat$education)
#Plots
#BMI.........
#Density plot all in one
densityplot(~bmi, data = dat,groups = study,  plot.points = F, 
auto.key = list(space="bottom", columns=2))
#Density plot separate
densityplot(~bmi | study, data=dat, col="blue")
#Histogram separate
histogram(~bmi | study, data=dat, breaks=seq(from=0,to=70,by=0.6))
#Box plot separate
bwplot(bmi~study, data = dat)
#Density plot all
densityplot(~bmi, data=dat, col="blue")
#Box plot all
bwplot(~bmi, data = dat)
#Menopause Age.......
#Density plot all in one
densityplot(~menopauseage, data = dat,groups = study,  plot.points = F, 
auto.key = list(space="bottom", columns=2))
#Density plot separate
densityplot(~menopauseage | study, data=dat, col="blue")
#Histogram
histogram(~menopauseage | study, data=dat, breaks=seq(from=0,to=90,by=0.6))
#Box plot
bwplot(menopauseage~study, data = dat)
#Density plot all studies
densityplot(~menopauseage, data=dat, col="blue")
#Box plot all studies
bwplot(~menopauseage, data = dat)

#Step 2: Assessing the linearity of the relationship---------------------------------------
#Each study
xyplot(menopauseage ~ bmi| study, data = dat, type =c("p", "smooth"), auto.key = TRUE)
#All studies
xyplot(menopauseage ~ bmi, data = dat, type =c("p", "smooth"), auto.key = TRUE)

#Step 3: Modelling the effect of confounding variables on the explanatory variable---------
d1 <- subset(dat, study==1); d2 <- subset(dat, study==2); d3 <- subset(dat, study==3)
modeld1 <- lm(bmi ~ smoking + nochild + education + age, data=d1); summary(modeld1)
modeld2 <- lm(bmi ~ smoking + nochild + education + age, data=d2); summary(modeld2)
modeld3 <- lm(bmi ~ smoking + nochild + education + age, data=d3); summary(modeld3)
s3pred <- predict.lm(modeld3); s2pred <- predict.lm(modeld2); s1pred <- predict.lm(modeld1)
predall <- c(s3pred, s2pred, s1pred)
length(s3pred); length(predall)
dat1 <- cbind(dat,predall)
names(dat1); head(dat1)
#Confounder variable model index, 
#i.e. the residuals, is used to adjust for confounders in our main model
dat1$cvmi <- dat1$bmi - dat1$predall
d1 <- subset(dat1, study==1); d2 <- subset(dat1, study==2); d3 <- subset(dat1, study==3)
dim(d1); dim(d2); dim(d3)
#Step 4: Selecting the best study-specific fractional polynomials-
#-for the association between the outcome and exposure variables------------------------

#Fit up to fourth order polynomial with powers in {-2, -1, -0.5, 0, 0.5, 1, 2, 3}, if a power repeats use log for the second,
#See the article reference number 7, and https://doi.org/10.1016/j.csda.2005.07.015. Alternatively can use the R package mfp if preferred.

#First order, select the one with highest likelihood. You also need to centre bmi at its mean and divide by 10 to do a reasonable model fit;
#For example, use bmi/10 - mean(bmi/10) instead of bmi
o1p1 <- lm(formula = menopauseage ~ I(bmi^(-2)) + cvmi, data = d1)
o1p2 <- lm(formula = menopauseage ~ I(bmi^(-1)) + cvmi, data = d1)
logLik(o1p1); logLik(o1p2)
#continue the process...for illustration for order 1 we chose the power 3 - FP1(3)
o1p8 <- lm(formula = menopauseage ~ I(bmi^(3)) + cvmi, data = d1)

#Repeat the process for order 2
o2p11 <- lm(formula = menopauseage ~ I(bmi^(-2)) + log(I(bmi^(-2))) + cvmi, data = d1)
o2p12 <- lm(formula = menopauseage ~ I(bmi^(-2)) + I(bmi^(-1)) + cvmi, data = d1)
logLik(o2p11); logLik(o2p12)
#continue the process...For illustration of this code for order 2 we chose the power (3,3) - FP2(3,3)
o2p88 <- lm(formula = menopauseage ~ I(bmi^(3)) + log(I(bmi^(3))) + cvmi, data = d1)
#Repeat the process for order 3
o3p112 <- lm(formula = menopauseage ~ I(bmi^(-2)) + log(I(bmi^(-2))) + I(bmi^(-1)) + cvmi, data = d1)
o3p118 <- lm(formula = menopauseage ~ I(bmi^(-2)) + log(I(bmi^(-2))) + I(bmi^(3)) + cvmi, data = d1)
logLik(o3p112); logLik(o3p118)
#By continuing the process we chose "FP3(3 3 3)" and "FP4(2 2 3 3)" for 3rd and 4th order FP
o3p888 <- lm(formula = menopauseage ~ I(bmi^(3)) + log(I(bmi^(3))) + log(I(bmi^(3))) + cvmi, data = d1)
o4p7788 <- lm(formula = menopauseage ~ I(bmi^(2)) + log(I(bmi^(2))) + I(bmi^(3)) + log(I(bmi^(3))) + cvmi, data = d1)

#Plot predicted value of all, FP1 to FP4 to select the best model
mytheme = theme(
  axis.title.x=element_text(size=7), 
  axis.title.y=element_text(size=7),
  axis.text.x = element_text(size=7),
  axis.text.y = element_text(size=7))
cord = coord_cartesian(xlim=c(0,70), ylim=c(39,65))

ggplot(d1, aes(x = bmi, y = menopauseage)) + geom_point()  + geom_line(aes(y=predict(o1p8), color="FP1(3)"), size=1) + 
  geom_line(aes(y=predict(o2p88), color="FP2(3,3)"), size=1) + 
  geom_line(aes(y=predict(o3p888), color="FP3(3,3,3)"), size=1) + 
  geom_line(aes(y=predict(o4p7788), color="FP4(2,2,3,3)"), size=1) +
  mytheme + cord

#For the illustration purpose of this pseudo code the model choice for study 1 is FP2(3 3);

#Repeat the above process for study 2 (this includes the model selection process above and the plot)
#For illustration purpose of this code the model choice for study 2 is FP1(3);
#Similarly for study 3 the model is FP3(0.5 3 3)

#Step 5: Estimating the functional forms for the association between the outcome and exposure variables for each study, adjusting for confounders
#Study 1 - FP2(3 3)
s1fit <- lm(formula = menopauseage ~ I(bmi^(3)) + log(I(bmi^(3))) + cvmi, data = d1)
s1p <- predict(s1fit, se.fit=T, newdata=dat1)
s1pval <- s1p$fit; s1pse <- s1p$se.fit; v1pse <- s1pse^2
#Study 2 - FP1(3)
s2fit <- lm(formula = menopauseage ~ I(bmi^(3)) + cvmi, data = d2)
s2p <- predict(s2fit, se.fit=T, newdata=dat2)
s2pval <- s2p$fit; s2pse <- s2p$se.fit; v2pse <- s2pse^2
#Study 3 - FP3(0.5,3,3)
s3fit <- lm(formula = menopauseage ~ I(bmi^(0.5)) + I(bmi^(3)) + log(I(bmi^(3))) + cvmi, data = d3)
s3p <- predict(s3fit, se.fit=T, newdata=dat3)
s3pval <- s3p$fit; s3pse <- s3p$se.fit; v3pse <- s3pse^2

#Step 6: Pooling the functional forms across studies
#Fixed effect weights
inv_vme1 <- ifelse(s1pse==0,NA,(v1pse)^(-1))
inv_vme2 <- ifelse(s2pse==0,NA,(v2pse)^(-1))
inv_vme3 <- ifelse(s3pse==0,NA,(v3pse)^(-1))

#sum of fixed effect weights
sumfixw <- inv_vme1 + inv_vme2 + inv_vme3

#Standardised fixed effect weights
sw1 <- inv_vme1/sumfixw
sw2 <- inv_vme2/sumfixw
sw3 <- inv_vme3/sumfixw

#Overall fixed effect estimate and the variance
fixme <- sw1*s1pval + sw2*s2pval + sw3*s3pval
vfixme <- ifelse(sumfixw==0,NA,(sumfixw)^(-1))

####Calculation random effect weights#########

Q <- inv_vme1*((s1pval-fixme)^2) + inv_vme2*((s2pval-fixme)^2) + inv_vme3*((s3pval-fixme)^2)

#Sum of inverse variance squared
suminvsq <- inv_vme1^2 + inv_vme2^2 + inv_vme3^2

#Tau Square
tausq <- ifelse(sumfixw==0,0,max(0,((Q-3+1)/(sumfixw - (suminvsq/sumfixw)))))

#sum of random-effect weights
sumrandw = (1/(inv_vme1 + tausq)) + (1/(inv_vme2 + tausq)) + (1/(inv_vme3 + tausq))

#Standardised random effect weight
srw1 = ((inv_vme1 + tausq)^(-1))/sumrandw
srw2 = ((inv_vme2 + tausq)^(-1))/sumrandw
srw3 = ((inv_vme3 + tausq)^(-1))/sumrandw

#Overall random-effect estimate and the variance
ranme = srw1*s1pval + srw2*s2pval + srw3*s3pval
varrandme = ifelse(sumrandw==0,.,(sumrandw)^(-1))
  
###########Plots

#Plot of the estimated function for each study 
datplot <- cbind(dat1,s1pval,s2pval,s3pval); dim(datplot)
d1plot <- subset(datplot, study==1); d2plot <- subset(datplot, study==2); d3plot <- subset(datplot, study==3)

mytheme = theme(
  axis.title.x=element_text(size=7), 
  axis.title.y=element_text(size=7),
  axis.text.x = element_text(size=7),
  axis.text.y = element_text(size=7))
cord = coord_cartesian(xlim=c(0,60), ylim=c(35,65))

pr1plot <- ggplot(d1plot, aes(x=bmi, y=s1pval)) + geom_point(aes(y=menopauseage), colour="lightgray", size=1.5) + 
  geom_smooth(aes(x=bmi, y=s1pval), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord

pr2plot <- ggplot(d2plot, aes(x=bmi, y=s2pval)) + geom_point(aes(y=menopauseage), colour="lightgray", size=1.5) + 
  geom_smooth(aes(x=bmi, y=s2pval), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord

pr3plot <- ggplot(d3plot, aes(x=bmi, y=s3pval)) + geom_point(aes(y=menopauseage), colour="lightgray", size=1.5) + 
  geom_smooth(aes(x=bmi, y=s3pval), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="") + mytheme + cord

ggdraw(plot=NULL, xlim=c(0,6), ylim=c(0,8), clip="off") + 
  draw_plot(pr1plot, x=0, y=6, width=2, height=2) +
  draw_plot(pr2plot, x=2, y=6, width=2, height=2) +
  draw_plot(pr3plot, x=4, y=6, width=2, height=2) + 
  
draw_plot_label(label= c("Study 1","Study 2", "Study 3"), size=6, x=c(0.25,2.27,4.30), y=c(8.005,8.005,8.005))


#######Plot the overall estimated function and the between study variance 
mytheme = theme(
  axis.title.x=element_text(size=7), 
  axis.title.y=element_text(size=7),
  axis.text.x = element_text(size=7),
  axis.text.y = element_text(size=7))

prallplot <- ggplot(datplot, aes(x=bmi)) + geom_point(aes(y=menopauseage), colour="lightgray", size=1.5) + 
  geom_smooth(aes(x=bmi, y=ranme), method="loess", se=F, size=0.8, colour="black") + labs(y="Menopause Age",x="BMI") + mytheme

ggdraw(plot=NULL, xlim=c(0,4), ylim=c(0,2), clip="off") + 
  draw_plot(prallplot, x=0, y=0, width=2, height=2) 


