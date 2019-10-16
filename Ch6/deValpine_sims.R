# de Valpine simulations
# recreate de Valpline's (2009) analysis showing you can get the same stage distribution with different growth distributions

require(varDev)
# this package is available at Ecological Archives E090-205-S1
# http://www.esapubs.org/archive/ecol/E090/205/suppl-1.htm

# First, create the matrix: ( the 3-stage cactus model)
VDM <- VD.model(3,marginal.durations=list(VD.dist("geomp1",list(prob=2/69)),
               VD.dist("geomp1",list(prob=14/98)),VD.dist("geomp1",list(prob=0.03))),
               marginal.death.times=list(VD.dist("geomp1",list(prob=0.31)),
               VD.dist("geomp1",list(prob=0.02)),VD.dist("infinite")))

# Now find combos with the same SSD, but different parameters:
#Means = seq(log10(7),log10(40),length.out =10)
n = 20
Shapes <- 10^(seq(-1.5,1,length.out=n))
Rates <- rep(NA,length(Shapes))
for (s in 1:length(Shapes)){
  Rates[s] = VD.solve.fraction.maturing.at.SASD.independent(VD.dist("gamma",
          list(shape = Shapes[s], rate = 1)),VD.dist("geomp1",list(prob=0.02)),
          r=-0.00024,fraction.maturing = 14/98,param="rate",range=c(1e-04,30),
          max.age=200)
}


Means = Shapes/Rates
SDs = sqrt( Shapes/(Rates)^2)

quartz(width=5,height=5)
plot(Means,SDs,log="xy",type='l',cex.axis=0.75)
points(Means[15],SDs[15])
points(Means[5],SDs[5])
quartz.save('Gamma_distributions.pdf',type='pdf')

# Plot the two extreme examples of the age distribution:
VDM1 <- VD.model(3,marginal.durations=list(VD.dist("geomp1",list(prob=2/69)),
                VD.dist("gamma",list(shape = Shapes[15],rate=Rates[15])),VD.dist("geomp1",list(prob=0.03))),
                marginal.death.times=list(VD.dist("geomp1",list(prob=0.31)),
                VD.dist("geomp1",list(prob=0.02)),VD.dist("infinite")))

VDM2 <- VD.model(3,marginal.durations=list(VD.dist("geomp1",list(prob=2/69)),
                                           VD.dist("gamma",list(shape = Shapes[5],rate=Rates[5])),VD.dist("geomp1",list(prob=0.03))),
                 marginal.death.times=list(VD.dist("geomp1",list(prob=0.31)),
                                           VD.dist("geomp1",list(prob=0.02)),VD.dist("infinite")))


# Create VDS objects & calculate r for each example:
VDS1 <- VD.run(VDM1)
dev.table1<- compile.dev.table(VDS1)
mean.fec1 <- calc.average.surv.rep.by.age(dev.table1,F=0.56)
r1 <- VD.solve.euler(mean.fec1)

VDS2 <- VD.run(VDM2)
dev.table2 <- compile.dev.table(VDS2)
mean.fec2 <- calc.average.surv.rep.by.age(dev.table2,F=0.56)
r2 <- VD.solve.euler(mean.fec2)


quartz(width=5,height=5)
SASD1 <- make.SASD(VDS1,r1)
plot(SASD1,age.vector=1:50)
quartz.save('SSD1.pdf',type='pdf')

quartz(width=5,height=5)
SASD2 <- make.SASD(VDS2,r2)
plot(SASD2,age.vector=1:50)
quartz.save('SSD2.pdf',type='pdf')

