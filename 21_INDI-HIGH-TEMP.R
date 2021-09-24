# NPZD model 21_INDI
#HIGH-TEMP scenario

#clear workspace
rm(list = ls())

#load packages
library(tidyverse)# for general plotting and data handling
library(rTPC)# for accessing functions for the thermal response curve
library(deSolve)# for numerical integration

#switch for saving the df, binary: TRUE or FALSE
write.df = TRUE


#set time components
secperday <- 86400     #seconds per day [seconds]
deltat    <- 6000      #delta t in [seconds]
ye <- 50               #timespan of temperature increase [years]
 










#set NPZD parameters as vector
parameters <- c(secperday=86400,
                kN21 = 0.35, #[half saturation constant/ -]
                gaP21 = 0.05/secperday, #[mortality rate strain 21.1/sec^-1]
                gpmax21 = ((spain_1982(temp=15, a = -3.509e-02, b= 4.663e-25, c=2.000e+00, r0=2.170e+00))/secperday), #[maximum gp rate/sec^-1]
                rmax21 = ((spain_1982(temp=15, a=-2.374e-02, b=4.429e-25, c=2.000e+00, r0=1.512e+00))/secperday), #[maximum r rate/sec^-1]
                gmax = 1/secperday, #[maximum grazing rate/sec^-1]
                gaGra = 0.25/secperday, #[mortality rate grazer/sec^-1]
                kP = 0.5, #[half saturation konstant grazer/-]
                beta = 0.7, #[assimilation coefficient/-]
                tau = 0.05/secperday, #[remineralization rate/sec^-1]
                deltatemp = (0.025*ye)/(ye*360*secperday), #[temperature increase per time step/°C sec^-1]]
                alpha = 0.03/secperday,  #[slope of light limiting function]
                deltat = 6000) #[delta t/sec^-1]

#set intitial values for state variables
state      <- c(light=200,temp=12.5, nut=4,detri=0,phyto21=3,gra=2)
#set time scale
times      <- seq(from=0, to=111.75*360*secperday, by = deltat)


## 13_INDI, set up model structure for the model
NPZD <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dlight <- ((2*pi*150*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    if(t<=30.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    if (t>30.75*360*secperday & t<=80.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))+deltatemp
    if(t>80.75*360*secperday) dtemp <- ((2*pi*7.5*cos((2*pi*t)/(360*secperday)))/(360*secperday))
    
    
    p21 <- gpmax21*(spain_1982(temp=temp,a = -3.509e-02, b= 4.663e-25, c=2.000e+00, r0=2.170e+00)/(gpmax21*secperday))*(nut/(kN21+nut))*tanh((alpha*light)/gpmax21)
    r21 <- rmax21*(spain_1982(temp=temp, a=-2.374e-02, b=4.429e-25, c=2.000e+00, r0=1.512e+00)/(rmax21*secperday))
    
    

    growth21 <- (p21-r21)*phyto21
    grazing <- gmax*((phyto21)/(kP+(phyto21)))*gra
    
    dnut <- (  - growth21 + tau*detri)
    
    ddetri <- (gaP21*phyto21+(1-beta)*grazing+ gaGra*gra-tau*detri)
    
    dphyto21 <- (growth21-gaP21*phyto21- grazing)
    
    dgra <- (beta*grazing - gaGra*gra)
    
    return(list(c(dlight,dtemp, dnut, ddetri, dphyto21, dgra)))
  }) 
}






# calculating model via runge kuttha 4th order integration
out.rk4 <- rk4(y = state, times = times, func = NPZD, parms = parameters)
#store model output as df
out.rk4.df <- as.data.frame(out.rk4)

#save model output as data
if(write.df == TRUE) write.csv(out.rk4.df, "directory/filename.csv")

#plot temp and light
ggplot(out.rk4.df, aes(x=time/secperday/360, y=temp))+
  geom_line()+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=light), col="yellow")
#plot full model
ggplot(out.rk4.df, aes(x=time/secperday/360, y=phyto21))+
  geom_line(col="green")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=detri), col="grey")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=nut), col="lightblue")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=gra), col="red")+
  theme_bw()+
  labs(x="time[days]", y="concentration", title="NPZD 21.1")+
  geom_line(data=out.rk4.df, aes(x=time/secperday/360, y=temp), col="purple")
