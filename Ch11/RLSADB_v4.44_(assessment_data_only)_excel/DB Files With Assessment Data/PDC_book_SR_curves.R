# SR plots for PDC book, using RAM legacy database

# Must run from working directory containing DBdata.Rdata
library(ggplot2)
library(gridExtra)
library(cowplot)


# load data:
load("DBdata.RData")

Y <- seq(1872,2018) # total length of timeseries

# Get YT Rockfish
YT.ssb <- ssb.data$YTROCKNPCOAST
YT.r <-r.data$YTROCKNPCOAST

# what years are these for?
#timeseries_years_views$SSB[timeseries_years_views$stockid=="YTROCKNPCOAST"]
YT.y = c(1967,2005)
YT.steep <- c(0.25,0.35, 0.65, 0.95) # overall posterior mean (low information)
YT.R0 <- mean(YT.r[Y > YT.y[1] & Y < (YT.y[1]+10)])
#YT.phi <- 2.403

# make a dataframe
YT <- data.frame(YT.ssb,YT.r)

# BH curve
YT.x <- seq(0,max(YT.ssb,na.rm=TRUE),length.out=1e5)
YT.phi <- max(YT.x)/YT.R0
YT.BH = matrix(data=NA,nrow=length(YT.x),ncol=length(YT.steep))
for (i in 1:length(YT.steep)){
YT.BH[,i] <- 0.8*YT.R0*YT.steep[i]*YT.x/(0.2*YT.phi*YT.R0*(1-YT.steep[i]) +
         (YT.steep[i] - 0.2)*YT.x) }
YT.dummy <- data.frame(YT.x,YT.BH)

# Get POP
POP.ssb <- ssb.data$POPERCHPCOAST
POP.r <-r.data$POPERCHPCOAST

# what years are these for?
#timeseries_years_views$SSB[timeseries_years_views$stockid=="YTROCKNPCOAST"]
POP.y = c(1940,2011)
POP.steep <- c(0.25,0.35, 0.65, 0.95) 
POP.R0 <- mean(POP.r[Y > POP.y[1] & Y < (POP.y[1]+10)])

# make a dataframe
POP <- data.frame(POP.ssb,POP.r)

# BH curve
POP.x <- seq(0,max(POP.ssb,na.rm=TRUE),length.out=1e5)
POP.phi <- max(POP.x)/POP.R0
POP.BH = matrix(data=NA,nrow=length(POP.x),ncol=length(POP.steep))
for (i in 1:length(POP.steep)){
POP.BH[,i] <- 0.8*YT.R0*POP.steep[i]*POP.x/(0.2*POP.phi*POP.R0*(1-POP.steep[i]) +
                                   (POP.steep[i] - 0.2)*POP.x) }
POP.dummy <- data.frame(POP.x,POP.BH)


# Get Chilipepper data
CP.ssb <- ssb.data$CHILISPCOAST
CP.f <- f.data$CHILISPCOAST
CP.Fmsy <- bioparams_values_views$Fmsy[bioparams_values_views$stockid=="CHILISPCOAST"]
CP.SSBmsy <- bioparams_values_views$SSBmsy[bioparams_values_views$stockid=="CHILISPCOAST"]
CP.SSBlim <- bioparams_values_views$SSBlim[bioparams_values_views$stockid=="CHILISPCOAST"]

# NEED TO FIX THIS. NEED YEARS, PLUS DATA
CP.f <- CP.f[72:218]
CP.f <- 1-CP.f
CP <- data.frame(Y,CP.ssb,CP.f)
CP2 <- data.frame(c(Y[1],Y[147]),rep(CP.Fmsy,2),rep(CP.SSBmsy,2),rep(CP.SSBlim,2))
names(CP2)<-list('Y','Fmsy','SSBmsy','SSBlim')

# Now do some plotting
YT.gg <- ggplot(YT,aes(x=YT.ssb,y=YT.r))+
  geom_point(size=2,shape=21,color='black',fill=NA)+
  geom_path(size=0.5)+
  geom_line(data=YT.dummy,aes(x=YT.x,y=X1),size=1,color='gray')+
  geom_line(data=YT.dummy,aes(x=YT.x,y=X2),size=1,color='gray')+
  geom_line(data=YT.dummy,aes(x=YT.x,y=X3),size=1,color='black')+
  geom_line(data=YT.dummy,aes(x=YT.x,y=X4),size=1,color='gray')+
  xlab('Spawning stock biomass (kg)')+
  ylab('Recruitment (No. fish)')+
  theme_classic()

POP.gg <- ggplot(POP,aes(x=POP.ssb,y=POP.r))+
  geom_point(size=2,shape=21,color='black',fill=NA)+
  geom_path(size=0.5)+
  geom_line(data=POP.dummy,aes(x=POP.x,y=X1),size=1,color='gray')+
  geom_line(data=POP.dummy,aes(x=POP.x,y=X2),size=1,color='black')+
  geom_line(data=POP.dummy,aes(x=POP.x,y=X3),size=1,color='gray')+
  geom_line(data=POP.dummy,aes(x=POP.x,y=X4),size=1,color='gray')+
  xlab('Spawning stock biomass (kg)')+
  ylab('Recruitment (No. fish)')+
  theme_classic()



quartz(width=3.5,height=6)
grid.arrange(YT.gg,POP.gg, ncol=1, nrow=2)

# Plot Chilipepper data

CP.p1 <- ggplot(CP,aes(x=Y,y=CP.ssb/1e9))+
  geom_point(size=1,shape=21,color='black',fill=NA)+
  geom_line(size=0.5)+
  geom_hline(yintercept=CP.SSBmsy/1e9,size=0.5,lty=2,color='gray')+
  geom_hline(yintercept=CP.SSBlim/1e9,size=0.5,lty=2,color='gray')+
  ylab('Spawning stock biomass (10^9 kg)')+
  ylim(c(0,7.1))+
  scale_x_continuous(breaks=seq(1900,2010,10),limits=c(1930,2020))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black",size=8))+
  theme(axis.text.y = element_text(color="black",size=8))

CP.p2 <- ggplot(CP,aes(x=Y,y=CP.f))+
  geom_point(size=1,shape=21,color='black',fill=NA)+
  geom_line(size=0.5)+
  geom_abline(data=CP2,slope=0,intercept=CP2$Fmsy,size=0.5,lty=2,color='gray')+
  ylab('Fishing mortality (y-1)')+
  ylim(c(0,1))+
  scale_x_continuous(breaks=seq(1900,2010,10),limits=c(1930,2020),trans=)+
  theme_classic()+
  theme(axis.text.x = element_text(color="black",size=8))+
  theme(axis.text.y = element_text(color="black",size=8))

Y2 = as.character(Y)
Y2 = Y2[Y>=1930]
CP.p3 <- ggplot(CP[CP$Y>=1930,],aes(y=CP.f,x=CP.ssb))+
  geom_point(size=1,shape=21,color='black',fill=NA)+
  geom_text(aes(y=CP.f,x=CP.ssb,label=ifelse(Y %% 10 == 0,Y2,'')),hjust=-0.1,size=3)+
  geom_path(size=0.5)+
  geom_hline(yintercept=CP2$Fmsy,size=0.5,lty=2,color='gray')+
  geom_vline(xintercept=CP2$SSBmsy,size=0.5,lty=2,color='gray')+
  geom_vline(xintercept=CP2$SSBlim,size=0.5,lty=2,color='gray')+
  ylab('Fishing mortality (y-1)')+
  xlab('Spawning stock biomass (kg)')+
  ylim(c(0,1))+
  xlim(c(0,7.1e9))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black",size=8))+
  theme(axis.text.y = element_text(color="black",size=8))
  
  
quartz(width=3.5,height=6)
plot_grid(CP.p1, CP.p2, CP.p3, align = "v", nrow = 3)

  

