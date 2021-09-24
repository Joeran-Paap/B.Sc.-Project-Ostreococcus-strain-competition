#N2PZD model with both strains
#MID-TEMP scenario

#clear workspace
rm(list = ls())


#load packages
library(tidyverse)# for general plotting and data handling
library(rTPC)# for accessing functions for the thermal response curve
library(deSolve)# for numerical integration



#set time components
secperday <- 86400     #seconds per day [seconds]
deltat    <- 6000      #delta t in [seconds]
ye <- 50               #timespan of temperature increase [years]
  








#set N2PZD parameters as vector
parameters <- c(secperday=86400,
                kN13 = 0.5, #[half saturation constant strain 13.1/ -]
                kN21 = 0.35, #[half saturation constant strain 21.1/ -]
                gaP13 = 0.05/secperday, #[mortality rate strain 13.1/sec^-1]
                gaP21 = 0.05/secperday, #[mortality rate strain 21.1/sec^-1]
                gpmax13 = ((johnsonlewin_1946(temp=2.701e+01, r0=8.393e+06, e=3.786e-01, eh=4.000e+01, topt=2.701e+01))/secperday), #[maximum gp rate strain 13.1/sec^-1]
                gpmax21 = ((spain_1982(temp=15, a = -3.509e-02, b= 4.663e-25, c=2.000e+00, r0=2.170e+00))/secperday), #[maximum gp rate strain 21.1/sec^-1]
                rmax13 = (spain_1982(temp=26.192, a=4.493e-02, b=3.909e-25, c=2.000e+00, r0=7.273e-01))/secperday, #[maximum r rate strain 13.1/sec^-1]
                rmax21 = ((spain_1982(temp=15, a=-2.374e-02, b=4.429e-25, c=2.000e+00, r0=1.512e+00))/secperday), #[maximum r rate strain 21.1/sec^-1]
                gmax = 1/secperday, #[maximum grazing rate/sec^-1]
                gaGra = 0.25/secperday, #[mortality rate grazer/sec^-1]
                kP = 0.5, #[half saturation konstant grazer/-]
                beta = 0.7, #[assimilation coefficient/-]
                tau = 0.05/secperday, #[remineralization rate/sec^-1]
                deltatemp = (0.018*ye)/(ye*360*secperday), #[temperature increase per time step/°C sec^-1]]
                alpha = 0.03/secperday, #[slope of light limiting function]
                deltat = 6000) #[delta t/sec^-1]


#set intitial values for state variables
state      <- c(light=200,temp=12.5, nut=4,detri=0,phyto13=3,phyto21=3,gra=2)
#set time scale
times      <- seq(from=0, to=111.75*360*secperday, by = deltat)


## N2PZD, set up model structure for the model
NPZD <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dlight <- ((2*pi*150*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    if(t<=30.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    if (t>30.75*360*secperday & t<=80.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))+deltatemp
    if(t>80.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    
    p13 <- gpmax13*(johnsonlewin_1946(temp=temp, r0=8.393e+06, e=3.786e-01, eh=4.000e+01, topt=2.701e+01)/(gpmax13*secperday))*(nut/(kN13+nut))*tanh((alpha*light)/gpmax13)
    r13 <- rmax13*(spain_1982(temp=temp, a=4.493e-02, b=3.909e-25, c=2.000e+00, r0=7.273e-01)/(rmax13*secperday))
    
    p21 <- gpmax21*(spain_1982(temp=temp,a = -3.509e-02, b= 4.663e-25, c=2.000e+00, r0=2.170e+00)/(gpmax21*secperday))*(nut/(kN21+nut))*tanh((alpha*light)/gpmax21)
    r21 <- rmax21*(spain_1982(temp=temp, a=-2.374e-02, b=4.429e-25, c=2.000e+00, r0=1.512e+00)/(rmax21*secperday))
    
    
    growth13 <- (p13-r13)*phyto13
    growth21 <- (p21-r21)*phyto21
    
    grazing <- gmax*((phyto13+phyto21)/(kP+(phyto13+phyto21)))*gra
    a <- (phyto13/(phyto13+phyto21))

    dnut <- ( - growth13 - growth21 + tau*detri)

    ddetri <- (+gaP13*phyto13+gaP21*phyto21+(1-beta)*grazing+ gaGra*gra-tau*detri)

    dphyto13 <- (growth13-gaP13*phyto13- a*grazing)
    dphyto21 <- (growth21-gaP21*phyto21- (1-a)*grazing)

    dgra <- (beta*grazing - gaGra*gra)

    return(list(c(dlight,dtemp, dnut, ddetri,dphyto13, dphyto21, dgra)))
  }) 
}



#runge kuttha 4th order integration 
out.rk4 <- rk4(y = state, times = times, func = NPZD, parms = parameters)
out.rk4.df <- as.data.frame(out.rk4)


#save the df
write.csv(out.rk4.df, "/home/zmaw/u301231/bsc_project/skripts/out.rk4.df.temp_midtemp.csv")

#calculate time in days and years
out.rk4.df <- out.rk4.df %>% mutate(days = time/secperday, years = time/secperday/360)

#plot full time scale
ggplot(out.rk4.df, aes(x=time/secperday/360, y=phyto13))+
  geom_line()+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=phyto21), col="green")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=detri), col="grey")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=nut), col="lightblue")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=gra), col="red")+
  theme_bw()+
  labs(x="time[days]", y="concentration", title="NPZD RK4 - temp increase=0.4/30years")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=temp), col="purple")

#plot first half - no temp increase
ggplot(out.rk4.df[1:36001,], aes(x=time/secperday, y=phyto13))+
  geom_line()+
  geom_line(data=out.rk4.df[1:36001,], aes(x=time/secperday, y=phyto21), col="green")+
  geom_line(data=out.rk4.df[1:36001,], aes(x=time/secperday, y=detri), col="grey")+
  geom_line(data=out.rk4.df[1:36001,], aes(x=time/secperday, y=nut), col="lightblue")+
  geom_line(data=out.rk4.df[1:36001,], aes(x=time/secperday, y=gra), col="red")+
  theme_bw()+
  labs(x="time[days]", y="concentration", title="NPZD RK4")+
  geom_line(data=out.rk4.df[1:36001,], aes(x=time/secperday, y=temp), col="purple")

#plot second half - with temp increase
ggplot(out.rk4.df[36001:86761,], aes(x=time/secperday, y=phyto13))+
  geom_line()+
  geom_line(data=out.rk4.df[36001:86761,], aes(x=time/secperday, y=phyto21), col="green")+
  geom_line(data=out.rk4.df[36001:86761,], aes(x=time/secperday, y=detri), col="grey")+
  geom_line(data=out.rk4.df[36001:86761,], aes(x=time/secperday, y=nut), col="lightblue")+
  geom_line(data=out.rk4.df[36001:86761,], aes(x=time/secperday, y=gra), col="red")+
  theme_bw()+
  labs(x="time[days]", y="concentration", title="NPZD RK4 - temp increase")+
  geom_line(data=out.rk4.df[36001:86761,], aes(x=time/secperday, y=temp), col="purple")
