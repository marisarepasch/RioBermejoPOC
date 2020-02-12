# Model for source-to-sink carbon cycling in the Rio Bermejo
# Only modeling the aged fluvial POC pool during downstream transit
# objective is to replicate the measured values from the river sediment
# This is a 5-pool series linear model
# Each pool represents the stock of aged POC and its radiocarbon content at each of five distances along the channel

#install.packages("devtools")
#install.packages("FME", repos="http://R-Forge.R-project.org")
#devtools::install_github("MPIBGC-TEE/SoilR-exp/pkg")

library(SoilR)
library(forecast)
library(FME)
source('C:/Users/repasch/Dropbox/R/SoilR_Bermejo/SeriesPools14.R')
#source('/Users/marisa/Dropbox/R/SoilR_Bermejo/SeriesPools14.R')

## Prepare your dataset (actual measured values for river samples)
ber_Cdat = read.csv('C:/Users/repasch/Dropbox/Bulk_OC_data/ber_Cdat_indiv.csv', header = TRUE, sep = ",", quote = "\"",
                       dec = ".", fill = TRUE, comment.char = "")
# ber_Cdat = read.csv('/Users/marisa/Dropbox/Bulk_OC_data/ber_Cdat_indiv.csv', header = TRUE, sep = ",", quote = "\"",
#                    dec = ".", fill = TRUE, comment.char = "")
ber_Cdat_DI = read.csv('C:/Users/repasch/Dropbox/Bulk_OC_data/ber_Cdat_DI.csv', header = TRUE, sep = ",", quote = "\"",
                    dec = ".", fill = TRUE, comment.char = "")

ber_C14dat = read.csv('C:/Users/repasch/Dropbox/Bulk_OC_data/ber_C14dat_indiv.csv', header = TRUE, sep = ",", quote = "\"",
                    dec = ".", fill = TRUE, comment.char = "")
# ber_C14dat = read.csv('/Users/marisa/Dropbox/Bulk_OC_data/ber_C14dat_indiv.csv', header = TRUE, sep = ",", quote = "\"",
#                    dec = ".", fill = TRUE, comment.char = "")
ber_C14dat_DI = read.csv('C:/Users/repasch/Dropbox/Bulk_OC_data/ber_C14dat_DI.csv', header = TRUE, sep = ",", quote = "\"",
                      dec = ".", fill = TRUE, comment.char = "")

dist = c(0,135,422,866,1221)

# berdata = data.frame(
#   time = c(-6503,-6503+350,-6503+840,-6503+5330,-6503+8520),
#   Ct = c(AgedPOCdata[,2]),
#   Ct_SD = c(AgedPOCdata[,3]),
#   d14C = c(AgedPOCdata[,4]), #kg
#   d14C_SD = c(AgedPOCdata[,5]) #kg
# )

# plot measured data to visualize
plot(ber_Cdat_DI[,2],ber_Cdat_DI[,3],type="p",xlab="Transit time (yr)",
     ylab="C stock (kg)")
arrows(ber_Cdat_DI[,2],ber_Cdat_DI[,3]-ber_Cdat_DI[,4],ber_Cdat_DI[,2],ber_Cdat_DI[,3]+
         ber_Cdat_DI[,4],angle=90, length=0.1,code=3)

plot(ber_C14dat_DI[,2],ber_C14dat_DI[,3],type="p",xlab="Transit time (yr)",
     ylab="d14C",  ylim=c(-800,100))
arrows(ber_C14dat_DI[,2],ber_C14dat_DI[,3]-ber_C14dat_DI[,4],ber_C14dat_DI[,2],ber_C14dat_DI[,3]+
         ber_C14dat_DI[,4],angle=90, length=0.1,code=3)

## Define input variables
yr = seq(-6383,2017) # years over which the model will run (8500 yrs)
C0 = c(87138,0,0,0,0) #kg (initial aged POC flux at confluence)
F0_Delta14C = c(-167.91,-167.91,-167.91,-167.91,-167.91) #per mille (initial aged POC d14C at confluence (modeled value))
topsoil = c(33300,4050,24220,29104,6667) #total litter C stock avaialble (kg/yr)
litter = c(55500,16200,6055,7276,167) #total topsoil C stock available (kg/yr)

## Create atmospheric 14C curve for the last 50 kyr
# Bind IntCal13 (prebomb) with Hua (2013)
Fatm = bind.C14curves(prebomb=IntCal13,postbomb=Hua2013$SHZone12,time.scale="AD")
yrs=seq(1966,2009.5,by=1/4) # A series of years by quarters
sz12=spline(Hua2013$SHZone12[,c(1,4)],xout=yrs) #Spline interpolation of the NH_Zone 2 dataset at a quaterly basis
shz12=ts((sz12$y-1)*1000,start=1966,freq=4) #Transformation into a time-series object
m=ets(shz12) #Fits an exponential smoothing state space model to the time series
f2=forecast(m,h=11*4) #Uses the fitted model to forecast 11 years into the future
# bind forecasted 14C curve with Fatm composite curve
Fatm_2=data.frame(Year=c(Fatm[-dim(Fatm)[1],1],
                     seq(tsp(f2$mean)[1],tsp(f2$mean)[2], by=1/tsp(f2$mean)[3])),
              Delta14C=c(Fatm[-dim(Fatm)[1],2],as.numeric(f2$mean)))
plot(Fatm_2,type="l",xlim=c(-6383,2020))

# set up the model
berPOCfunc = function(pars) {
  mdl=SeriesPools14(
    t = yr,
    m.pools = 5,
    ki = pars[1:5],
    Tij = c(inipars[6],inipars[7],inipars[7],inipars[7]), #rep(pars[6], 4)
    C0 = C0,
    F0_Delta14C = F0_Delta14C,
    In = c(pars[8]*litter[1]+pars[13]*topsoil[1],pars[9]*litter[2]+pars[14]*topsoil[2],pars[10]*litter[3]+pars[15]*topsoil[3],
    pars[11]*litter[4]+pars[16]*topsoil[4],pars[12]*litter[5]+pars[17]*topsoil[5]),
    xi = 1,
    inputFc=Fatm_2,
    lambda = -0.0001209681,
    lag = 100,
    solver = deSolve.lsoda.wrapper,
    pass = FALSE
  )
  Ct = getC(mdl)
  C14t = getF14(mdl)
  CO2t = getAccumulatedRelease(mdl)
  ind = which(yr==2017)
  return(data.frame(time=c(-6503,-6503+350,-6503+840,-6503+5330,-6503+8520),C_kg=Ct[ind,],d14C=C14t[ind,],CO2=CO2t[ind,]))
  #return(data.frame(time=t,Ct=Ct,C14t=C14t,CO2t=CO2t))
}

# set the initial parameters for decomposition rates (ki), transfer proportions between pools (Ti), 
# and the fraction of leaf litter transferred to aged POC pool (gam1) and fraction of topsoil
# transferred to aged POC pool (gam2)
inipars = c(k1=1/5000,k2=1/3595,k3=1/1000,k4=1/2740,k5=1/3260,T1=0.75,T2=0.9,gamL1=0.3,
            gamL2=0.4,gamL3=0.6,gamL4=0.2,gamL5=0.3,gamS1=0.3,gamS2=0.3,gamS3=0.6,gamS4=0.3,gamS5=0.6)

test_2017_17pars = berPOCfunc(inipars)
head(test_2017_17pars)

plot(ber_Cdat_DI[,2],ber_Cdat_DI[,3],type="p",xlab="Year",
     ylab="C (kg)",ylim=c(1.5E7,4.0E8))
arrows(ber_Cdat_DI[,2],ber_Cdat_DI[,3]-ber_Cdat_DI[,4],ber_Cdat_DI[,2],ber_Cdat_DI[,3]+
         ber_Cdat_DI[,4],angle=90, length=0.1,code=3)
lines(test_2017_17pars[,1],test_2017_17pars[,2], type="l", lty=1, col='red', xlab = 'year', ylab = 'OC load (kg)')

plot(ber_C14dat_DI[,2],ber_C14dat_DI[,3],type="p",xlab="Year",
     ylab="d14C (per mille)",ylim=c(-800,100))
arrows(ber_C14dat_DI[,2],ber_C14dat_DI[,3]-ber_C14dat_DI[,4],ber_C14dat_DI[,2],ber_C14dat_DI[,3]+
         ber_C14dat_DI[,4],angle=90, length=0.1,code=3)
lines(test_2017_17pars[,1],test_2017_17pars[,3], type="l", lty=1, col='red')

CO2_cum = cumsum(test_2017_17pars[,4])
plot(dist,CO2_cum)

# here I re-organize the measured/observed data for the cost function:
C_obs = data.frame(time=ber_Cdat[,2],C_kg=ber_Cdat[,3],C_kg_SD=ber_Cdat[,4])
d14C_obs = data.frame(time=ber_C14dat[,2],d14C=ber_C14dat[,3],d14C_SD=ber_C14dat[,4])

C_obs_DI = data.frame(time=ber_Cdat_DI[,2],C_kg=ber_Cdat_DI[,3],C_kg_SD=ber_Cdat_DI[,4])
d14C_obs_DI = data.frame(time=ber_C14dat_DI[,2],d14C=ber_C14dat_DI[,3],d14C_SD=ber_C14dat_DI[,4])

# Set up the cost function...maybe my observed data are not organized correctly?
berCostFun = function(pars){
  modelOutput=berPOCfunc(inipars)
  cost1 = modCost(model=modelOutput, obs = C_obs, err = 'C_kg_SD') #Soil stocks for fitting
  cost2 = modCost(model=modelOutput, obs = d14C_obs, cost = cost1, err = 'd14C_SD') #14C data for fitting, with soil stocks cost included
  return(cost2)
}

# Fit the model parameters to the data, with least-error estimation
berFit_C_14C = modFit(f=berCostFun, p=inipars, method='Nelder-Mead', upper=1, 
                      lower=0)

print(berFit_C_14C$par) #Best fit parameters
print(berFit_C_14C$var_ms) #Some error stats

# Re-run fitted model with parameters fit above
fitMod = SeriesPools14(
    t=yr,
    m.pools=5,
    ki=berFit_C_14C$par[1:5],
    Tij=c(berFit_C_14C$par[6],berFit_C_14C$par[7],berFit_C_14C$par[7],berFit_C_14C$par[7]), #rep(pars[6], 4)
    C0=C0,
    F0_Delta14C=F0_Delta14C,
    In=c(berFit_C_14C$par[8]*litter[1]+berFit_C_14C$par[13]*topsoil[1],berFit_C_14C$par[9]*litter[2]+berFit_C_14C$par[14]*topsoil[2],berFit_C_14C$par[10]*litter[3]+berFit_C_14C$par[15]*topsoil[3],
         berFit_C_14C$par[11]*litter[4]+berFit_C_14C$par[16]*topsoil[4],berFit_C_14C$par[12]*litter[5]+berFit_C_14C$par[17]*topsoil[5]),
    xi = 1,
    inputFc=Fatm_2,
    lambda = -0.0001209681,
    lag = 100,
    solver = deSolve.lsoda.wrapper,
    pass = FALSE
  )
fit_Ct = getC(fitMod)
fit_C14t = getF14(fitMod)
fit_CO2t = getAccumulatedRelease(fitMod)
ind = which(yr==2017)
fit_C_2017 = fit_Ct[ind,]
fit_C14_2017 = fit_C14t[ind,]
fit_CO2_2017 = fit_CO2t[ind,]

plot(yr,fit_CO2t[,1], type="l", lty=1, col=1:5, xlab = 'year', ylab = 'Cumulative CO2 release (kgC)', 
     legend('topleft', legend = colnames(fit_CO2t)[1:5], col = 1:5, lty = 1))

#Plot the results of the model run with the best fit parameters
plot(ber_C14dat[,1],ber_C14dat[,3],type="p",xlab="Distance (km)",
     ylab="d14C",  ylim=c(-800,100))
arrows(ber_C14dat[,1],ber_C14dat[,3]-ber_C14dat[,4],ber_C14dat[,1],ber_C14dat[,3]+
         ber_C14dat[,4],angle=90, length=0.1,code=3)
lines(dist,fit_C14_2017)

plot(ber_Cdat_DI[,1],ber_Cdat_DI[,3],type="p",xlab="Distance (km)",
     ylab="C (kg)")
arrows(ber_Cdat_DI[,1],ber_Cdat_DI[,3]-ber_Cdat_DI[,4],ber_Cdat_DI[,1],ber_Cdat_DI[,3]+
         ber_Cdat_DI[,4],angle=90, length=0.1,code=3)
lines(dist,fit_C_2017,col='red')

plot(dist,fit_CO2_2017,col='blue')
lines(dist,fit_CO2_2017,col='blue')

save.image("C:/Users/repasch/Dropbox/R/SoilR_Bermejo/bermejo_5pool_17pars_DI.RData")

# Construct matrices to calculate System Age and Transit Time. "A" matrix should be a square with your -k values on diagonals,
# transfer coefficients on the left-side of k's, between the pools. Remember: a21 is "to, from"; to pool 2, from pool 1.
kbest <- c(berFit_C_14C$par[1], berFit_C_14C$par[2], berFit_C_14C$par[3], berFit_C_14C$par[4], berFit_C_14C$par[5])
a21best <- winchFit$par[3] * winchFit$par[1]
a32best <- winchFit$par[4] * winchFit$par[2]
a43best <- 
a54best <- 
A2=diag(-kbest) #This is some sloppy, carry-over code, haha don't judge
u = matrix(c(inputs, 0), ncol = 1)
A2[2,1] = a21best
A2[1,2] = a12best
SA = systemAge(A=A2, u=u3) #System Age
TT = transitTime(A=A2, u=u3) #Transit Time

SA$meanSystemAge
SA$meanPoolAge
SA$systemAgeDensity

TT$meanTransitTime
TT$transitTimeDensity

# Calculate Bayesian parameter sensitivity
niter <- 11000 #Number of iterations
burnin <- 7500
updatecov <- 100 #Somewhat arbitrary
outLength <- 1000

# This step takes a while!
bayes_fit <- modMCMC(f=berCostFun, p=c(berFit_C_14C$par), niter = niter, updatecov = updatecov,
                     var0 = winchFit$var_ms_unweighted, outputlength = outLength,
                     upper=c(10,10,1,1), lower=c(0,0,0,0), burninlength = burnin)

# Check stability of parameters
plot(bayes_fit)

# Check co-linearity
pairs(bayes_fit)

# Statistics on accepted parameter sets
pred_uncert <- sensRange(winchFunc, parInput = bayes_fit$pars, parRange = c(0,1), num = outLength)




