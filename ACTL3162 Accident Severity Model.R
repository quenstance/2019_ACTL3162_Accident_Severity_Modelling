#Set working directory
setwd("C:/Users/Quenstance Lau/OneDrive - UNSW/ACTL3162 General Insurance/Assignment")

#Load packages
library(fitdistrplus)
library(actuar)
library(moments)
library(rSymPy) 
library(rJava)
library(rootSolve)

#Import data
data<-read.csv("LossData.csv",header=FALSE)

########################################TASK 1############################################

#Examine data
str(data)
summary(data)
sd(data$V1) #8.552255
kurtosis(data$V1) #6.13269 
skewness(data$V1) #1.374291 

#Create directory to store plots
dir.create("figs")

#Visual data exploration
png(filename = "figs/Fig1.png",width=1000,height=600)
plotdist(data$V1,histo=TRUE,demp=TRUE,col="light blue")
mtext("Fig.1: Histogram and CDF plots of an empirical distribution",side=1,line=4) 
dev.off()

png(filename = "figs/Fig1.png",width=1000,height=600)
descdist(data$V1,boot=1000)
mtext("Fig.2: Skewness-Kurtosis plot ",side=1,line=4.1) 
dev.off()

#Fit right-tailed distributions
w<-fitdist(data$V1,"weibull")
g<-fitdist(data$V1,"gamma")
ln<-fitdist(data$V1,"lnorm")
b<-fitdist(data$V1,"burr", start=list(shape1=0.3,shape2=1,rate=1)) #guess start value 

#Summary Statistics for the fitted distributions
summary(w)
summary(g)
summary(ln)
summary(b)

#Plot 4 goodness-of-fit plots to visualize the fit
plot.legend <- c("Weibull", "lognormal","gamma","burr")

png(filename = "figs/Fig3.png",width=1000,height=1000)
denscomp(list(w,ln,g,b), legendtext = plot.legend)
mtext("Fig.3: Theorectical and Empirical Density Plot",side=1,line=4) 
dev.off()

png(filename = "figs/Fig4.png",width=1000,height=1000)
cdfcomp(list(w,ln,g,b), legendtext = plot.legend)
mtext("Fig.4: Theorectical and Empirical CDF",side=1,line=4) 
dev.off()

png(filename = "figs/Fig5.png",width=1000,height=1000)
qqcomp(list(w,ln,g,b), legendtext = plot.legend)
mtext("Fig.5: Q-Q plot of empirical quantiles (y-axis) against the theoretical quantiles (x-axis)",side=1,line=4) 
dev.off()

png(filename = "figs/Fig6.png",width=1000,height=1000)
ppcomp(list(w,ln,g,b), legendtext = plot.legend)
mtext("Fig.6: P-P plot of empirical (y-axis) against the theoretical quantiles (x-axis)",side=1,line=4) 
dev.off()

png(filename = "figs/Fig7.png",width=1000,height=1000)
cdfcomp(list(w,ln,g,b), legendtext = plot.legend, xlogscale = TRUE,ylogscale = TRUE)
mtext("Fig.7: Theorectical and Empirical CDF on a logscale",side=1,line=4) 
dev.off()

#GOODNESS OF FIT STATISTICS
gofstat(list(w,ln,g,b),fitnames=plot.legend)

#UNCERTAINTY IN BURR PARAMETERS ESTIMATE
b.boot<-bootdist(b,niter = 1001)
b.boot.summary<-summary(b.boot)
b.boot.summary

##Correlation plot of Burr Parameters estimate
png(filename = "figs/Fig8.png",width=1000,height=1000)
plot(b.boot,col="light blue")
mtext("Fig.8: Correlation plot of Burr parameters bootstrapped estimate",side=1,line=4)
dev.off()

##Extract the median Burr parameters estimate
b.shape1<-b.boot.summary$CI[1]
b.shape2<-b.boot.summary$CI[2]
b.rate<-b.boot.summary$CI[3]

##Update the new Burr distributions with the median Burr parameters estimate
b.new<-fitdist(data$V1,"burr", start=list(shape1=b.shape1,shape2=b.shape2,rate=b.rate))
summary(b.new)
gofstat(b.new)

#Fit a Pareto distribution
png(filename = "figs/Fig9.png",width=1000,height=1000)
pa<-fitdist(data$V1,"pareto", start=list(shape=0.3,scale=1))
plot.legend <- c("gamma","burr","pareto")
par(mfrow=c(2,2))
denscomp(list(g,b,pa), legendtext = plot.legend)
cdfcomp(list(g,b,pa), legendtext = plot.legend)
qqcomp(list(g,b,pa), legendtext = plot.legend)
ppcomp(list(g,b,pa), legendtext = plot.legend)
dev.off()

########################################TASK 2############################################
theta<-0.375

#Define SymPy variables
sympy("var('x')")
sympy("var('r')")

#Find the MGF
sympy("integrate(exp(x*r)*(exp(-x)/3+2*x*exp(-x)/3),(x,0,oo))")
  #[1] "9/(9 - 27*r + 27*r**2 - 9*r**3) - 12*r/(9 - 27*r + 27*r**2 - 9*r**3) + 3*r**2/(9 - 27*r + 27*r**2 - 9*r**3)"
#Simplify(function (r) 9/(9 - 27*r + 27*r**2 - 9*r**3) - 12*r/(9 - 27*r + 27*r**2 - 9*r**3) + 3*r**2/(9 - 27*r + 27*r**2 - 9*r**3))

#Find expected value
sympy("integrate(x*(exp(-x)/3 +2*x*exp(-x)/3),(x, 0, oo))")
  #[1] "5/3"

#PART A
#Create equation to solve for Lundberg coefficient
eqR.a<-function(r){
 (9-12*r+3*r^2)/(9-27*r+27*r^2-9*r^3)-1-((1+theta)*r*5/3)
}

#Graph it to see where the roots are located
png(filename = "figs/Fig10.png",width=1000,height=1000)
curve(eqR.a,xlim=c(0,0.3),col="blue")
abline(h=0,lty=1)
mtext("Fig.10: PART A- Locate roots to solve for Lundberg coefficients",side=1,line=4)
dev.off()

#Solve for R and set R>0
uniroot.all(eqR.a,c(0,1),tol=1e-15)
#[1] 0.000000 0.200000

#PART B: Solve for Mx(alpha*R)-1-1.135E[X]R
#PART Bi: alpha = a = 0.84
eqR.b<-function(r,a){
(9-12*(a*r)+3*(a*r)^2)/(9-27*(a*r)+27*(a*r)^2-9*(a*r)^3)-1-((((1+theta)-1.5*(1-a))*5/3)*r)
}

eqR.b<-function(r,a){
  (9-12*(a*r)+3*(a*r)^2)/(9-27*(a*r)+27*(a*r)^2-9*(a*r)^3)-1-((((1.375)-1.5*(1-a))*5/3)*r)
}

fR<-function(a){
  uniroot.all(eqR.b,c(0.0001,1),tol=1e-15,a=a)
}

fR(0.84) #[1]  0.2265746

fR.check<-function(a){
  uniroot(eqR.b,lower=0.0001,upper=1,a=a)$root 
}
  #NOTE: using lower = 0, will give root=0, and no other values

fR.check(0.84) #[1] 0.2265764

#PART Bii: finding the optimal alpha and corresponding R
optimise(fR,interval=c(0.01,0.999),maximum=TRUE)
  #$maximum [1] 0.457578 $objective [1] 0.2934871

#PART C
#Find E[h(x)]
sympy("integrate((x-3)*(exp(-x)/3+2*x*exp(-x)/3),(x,3,oo)") #11*exp(-3)/3

#Find new MGF
sympy("integrate(exp(r*x)*(exp(-x)/3 + 2*x*exp(-x)/3), (x, 0, 3))") 
  #(exp(3*r)*(7*r-9)+(3-r)*exp(3))/(3*exp(3)*(r-1)^2)
sympy("integrate(exp(3*r)*(exp(-x)/3 + 2*x*exp(-x)/3), (x,3,oo))")
  # 3*exp(-3)*exp(3*r)

eqR.c<-function(r){
  (exp(3*r)*(7*r-9)+(3-r)*exp(3))/(3*exp(3)*(r-1)^2)+3*exp(-3)*exp(3*r)-1-(1.375*5/3-1.5*(11/3)*exp(-3))*r
}

#Graph it to see where the roots are located
png(filename = "figs/Fig11.png",width=1000,height=1000)
curve(eqR.c,xlim=c(0,0.3),col="blue")
abline(h=0,lty=1)
mtext("Fig.11: Part C-Locate roots to solve for Lundberg coefficients",side=1,line=4)
dev.off()

#Solve for R and set R>0
R<-uniroot.all(eqR.c,c(0,1))
  #[1] 0.0000000 0.2687631
R<-max(R)
  #[1] 0.2687631

#PART D
c0<-1
d<-3
lowerbound<-exp(-R*(c0+d))
upperbound<-exp(-R*(c0))
lowerbound #[1] 0.3412833
upperbound #[1] 0.7643262
