# REVISED Dorleta Garcia 2022-11-28

# install.packages(c("ggplot2"))
# install.packages(c("FLa4a","FLash","FLXSA","FLBRP","ggplotFL"), repos="http://flr-project.org/R")
#### Loads all necessary packages #### 
library(FLa4a)
library(FLash)
library(FLXSA)
library(FLBRP)
library(ggplotFL)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
####  CONDITIONING THE OPERATING MODEL #### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-

# Read in stock assessment data ------------------------------------------------
data(ple4)
data(ple4.index)
stk <- ple4 
idx <- FLIndices(idx=ple4.index)

# Set up the iteration and projection window parameters (years) ----------------
it <- 20                    # iterations
y0 <- range(stk)["minyear"] # initial data year
dy <- range(stk)["maxyear"] # final data year
iy <- dy+1                  # initial year of projection (also intermediate year)
ny <- 12    # number of years to project from initial year
fy <- dy+ny # final year
nsqy <- 3   # number of years to compute status quo metrics


# ** Fit stock assessment model a4a ** -----------------------------------------

# Set up the catchability submodel with a smoothing spline (setting up a 'list' allows for more than one index) 
qmod <- list(~s(age, k=6))

# Set up the fishing mortality submodel as a tensor spline, which allows age and year to interact
fmod <- ~te(replace(age, age>9,9), year, k=c(6,8))

# Set up the MCMC parameters
mcsave <- 100
mcmc   <- it * mcsave

# Fit the model
fit <- sca(stk, idx, fmodel = fmod, qmodel = qmod, fit = "MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))

# Update the FLStock object
stk <- stk + fit
plot(stk)

# Reduce to keep one iteration only for reference points
stk0 <- qapply(stk, iterMedians)


# ** Fit the stock-recruit model ** --------------------------------------------
srbh  <- fmle(as.FLSR(stk,  model="bevholt"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
srbh0  <- fmle(as.FLSR(stk0,  model="bevholt"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
plot(srbh0)

# Segmented regression to 'estimate' Blim.
srsegreg0 <- fmle(as.FLSR(stk0, model="segreg"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf)) 

plot(srsegreg0)

# Generate stock-recruit residuals for the projection period
srbh.res <- rnorm(it, FLQuant(0, dimnames=list(year=iy:fy)), mean(c(apply(residuals(srbh), 6, sd))))

                                                                                    
## ** Calculate reference points ** --------------------------------------------                                                                     

# Calculate the reference points
brp  <- brp(FLBRP(stk0, srbh0))

refpts(brp)

Fmsy <- c(refpts(brp)["msy","harvest"])
msy  <- c(refpts(brp)["msy","yield"])
Bmsy <- c(refpts(brp)["msy","ssb"])
Blim <- srsegreg0@params$b[drop=T]
Bpa  <- Blim*exp(1.654*sqrt(var(ssb(stk)/mean(ssb(stk))))) 
  
# ** set up the operating model for the projection window ** -------------------
#Prepare the FLStock object for projections
stk <- stf(stk, fy-dy, nsqy, nsqy)

stk@stock.wt
stk@m

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
## SET UP OBSERVATION ERROR MODEL ELEMENTS -------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-

# Estimate the index catchabilities from the a4a  (without simulation)
# Set up the FLIndices object and populate it (note, FLIndices potentially has more than one index, hence the for loop)

idcs <- FLIndices()
for (i in 1:length(idx)){
  
  #   Set up FLQuants and calculate mean and sd for catchability
  lst        <- mcf(list(index(idx[[i]]), stock.n(stk0))) # make FLQuants same dimensions
  idx.lq     <- log(lst[[1]]/lst[[2]]) # log catchability of index
  idx.qmu    <- idx.qsig <- stock.n(iter(stk,1)) # create quants
  idx.qmu[]  <- yearMeans(idx.lq) # allocate same mean-at-age to every year
  idx.qsig[] <- sqrt(yearVars(idx.lq)) # allocate same sd-at-age to every year
  
  #   Build index catchability based on lognormal distribution with mean and sd calculated above
  idx.q    <- rlnorm(it, idx.qmu, idx.qsig)
  idx_temp <- idx.q * stock.n(stk)
  idx_temp <- FLIndex(index=idx_temp, index.q=idx.q) # generate initial index
  
  range(idx_temp)[c("startf", "endf")] <- c(0, 0) # timing of index (as proportion of year)
  idcs[[i]] <- idx_temp
}
names(idcs) <- names(idx)
idx<-idcs[1]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# SET UP MSE LOOP ---------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-

# Needed Functions:
# * Observation error model
# * XSA assessment model
# * Control object for projections

# * Observation error model ** -------------------------------------------------
o <- function(stk, idx, assessmentYear, dataYears) {
  # dataYears is a position vector, not the years themselves
  stk.tmp <- stk[, dataYears]
  # add small amount to avoid zeros
  catch.n(stk.tmp) <- catch.n(stk.tmp) + 0.1
   # Generate the indices - just data years
  idx.tmp <- lapply(idx, function(x) x[,dataYears])
  # Generate observed index
  for (i in 1:length(idx)) index(idx[[i]])[, assessmentYear] <-  stock.n(stk)[, assessmentYear]*index.q(idx[[i]])[, assessmentYear]
  
  return(list(stk=stk.tmp, idx=idx.tmp, idx.om=idx))
}

# * XSA assessment model ** ----------------------------------------------------
xsa <- function(stk, idx){
  # Set up XSA settings
  control  <- FLXSA.control(tol = 1e-09, maxit=99, min.nse=0.3, fse=2.0,
                            rage = -1, qage = range(stk)["max"]-1, shk.n = TRUE, shk.f = TRUE,
                            shk.yrs = 5, shk.ages= 5, window = 100, tsrange = 99, tspower = 0)
  # Fit XSA
  fit <- FLXSA(stk, idx, control)
  # convergence diagnostic (quick and dirty)
  maxit <- c("maxit" = fit@control@maxit)
  # Update stk
  stk   <- transform(stk, harvest = harvest(fit), stock.n = stock.n(fit))
  return(list(stk = stk, converge = maxit))
}

# * Control object for projections ** ------------------------------------------
getCtrl <- function(values, quantity, years, it){
  dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- unlist(values)
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA)) 
  ctrl@trgtArray <- arr0
  ctrl
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# MSE initialisation -----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-

vy <- ac(iy:fy)
TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(dy,vy), iter=1:it))
TAC[,ac(dy)] <- 140000
TAC[,ac(iy)] <- TAC[,ac(dy)] #assume same TAC in the first intermediate year
ctrl <- getCtrl(c(TAC[,ac(iy)]), "catch", iy, it)

# Set up the operating model FLStock object
stk.om <- fwd(stk, control=ctrl, sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# Start the MSE loop ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
set.seed(231) # set seed to ensure comparability between different runs

for(i in vy[-length(vy)]){
  
  # set up simulations parameters
  ay <- an(i)
  cat(i, " < ")
  flush.console()
  vy0 <- 1:(ay-y0)              # data years (positions vector) - one less than current year
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector) - one less than current year
  
  # apply observation error
  oem    <- o(stk.om, idx, i, vy0) 
  stk.mp <- oem$stk
  idx.mp <- oem$idx
  idx    <- oem$idx.om
  
  # perform assessment
  out.assess <- xsa(stk.mp, idx.mp)
  stk.mp     <- out.assess$stk
  
  # apply ICES MSY-like Rule to obtain Ftrgt (note this is not the ICES MSY rule, but is similar)
  flag  <- ssb(stk.mp)[,ac(ay-1)]<Bpa
  Ftrgt <- ifelse(flag,ssb(stk.mp)[,ac(ay-1)]*Fmsy/Bpa,Fmsy) 
  
  # project the perceived stock to get the TAC for ay+1
  fsq.mp    <- yearMeans(fbar(stk.mp)[,sqy]) # Use status quo years defined above
  ctrl      <- getCtrl(c(fsq.mp, Ftrgt), "f", c(ay, ay+1), it)
  stk.mp    <- stf(stk.mp, 2)
  gmean_rec <- c(exp(yearMeans(log(rec(stk.mp)))))
  stk.mp    <- fwd(stk.mp, control=ctrl, sr=list(model="mean", params = FLPar(gmean_rec,iter=it)))
  TAC[,ac(ay+1)] <- catch(stk.mp)[,ac(ay+1)]
  
  # apply the TAC to the operating model stock and project the population one year forward.
  ctrl   <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
  stk.om <- fwd(stk.om, control=ctrl,sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE) 
}

save(stk.om, file = 'C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# PERFORMANCE STATISTICS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# Some example performance statstics, but first isolate the projection period
stk.tmp <- window(stk.om,start=iy)

# annual probability of being below Blim 
risky.Bmsy <- iterSums((ssb(stk.tmp)/Bmsy)<1)/it
risky.Bpa  <- iterSums((ssb(stk.tmp)/Bpa)<1)/it
risky.Blim <- iterSums((ssb(stk.tmp)/Blim)<1)/it

risky.Bmsy 
risky.Bpa
risky.Blim

# mean probability of being below Bmsy in the first half of the projection period
mean(risky.Bmsy [,1:floor(length(risky.Bmsy)/2)])

# ...and second half 
mean(risky.Bmsy [,(floor(length(risky.Bmsy)/2)+1):length(risky.Bmsy)])

# plot of SSB relative to Bmsy
boxplot(data~year,data=as.data.frame(ssb(stk.tmp)),main="SSB")
abline(h=Bmsy,col="red")

plot(FLStocks(stk.om=stk.om, stk.mp=stk.mp)) + theme(legend.position="top") + geom_vline(aes(xintercept=2009))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-
# EXERCISES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-

# 1. Change the stock-recruitment relationship in the OM for a ricker model.
# 2. In the OM change the natural mortality, make it variable in the range [0.05,0.2]
# 3. With natural mortality variable in the OM, in the observation of the stock use M = 0.1
# 4. In the observation model add a log-normal error with median 1 and CV = 30% in the observation of the index.
# 5. In the implementation of the advice add an overshoot of the TAC using a uniform distribution with range [1,1.3] 


load('C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk_exe1.RData')
load('C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk_exe2.RData')
load('C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk_exe3.RData')
load('C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk_exe4.RData')
load('C:/users/dgarcia/Dropbox/FLBEIA_CursoIEO/Tutorials/02_MSE_with_FLR/stk_exe5.RData')


plot(FLStocks(stk.om = stk.om, 
           stk.om.sr = stk.om.sr, 
       stk.om.im.err =  stk.om.im.err, 
            stk.om.m = stk.om.m,
         stk.om.merr = stk.om.merr, 
   stk.om.obs.err.id = stk.om.obs.err.id)) + theme(legend.position="top") + geom_vline(aes(xintercept=2018))








## Not run: 
library(FLBEIA)
library(FLAssess)          # required to use the IcesHCR. Not available for win64
library(FLash)             # required to use the IcesHCR. Not available for win64
library(ggplot2)  


#---------------------------------------------------------------- 
# Example with 1 stock, 1 Fleets, 1 seasons and 1 iteration: one
#----------------------------------------------------------------

# Load the data to run FLBEIA in a one stock one fleet example using the HCR used by ICES 
# in the MSY framework. 
data(one) 

# The names and the class of the objects needed to run FLBEIA.
# sapply(ls(), function(x) class(get(x)))

# In this scenario a single, age-structured, stock is exploited by a single fleet with a 
# unique metier. 
# The fleet catches yearly exactly the adviced TAC and there is no exit-entry of vessels 
# in the fishery.  
# The stock abundance and exploitation level is observed without error in the observation 
# model.
# There is no assessment model and the TAC advice is used through the HCR used by ICES 
# in the MSY framework.  




s0 <- FLBEIA(biols = oneBio,   # FLBiols object with one FLBiol element for stk1.
             SRs = oneSR,    # A list with one FLSRsim object for stk1.
             BDs = NULL,     # No Biomass Dynamic populations in this case.
             fleets = oneFl,    # FLFleets object with on fleet.
             covars = oneCv,    # covars not used
             indices = NULL,     # indices not used 
             advice = oneAdv,   # A list with two elements 'TAC' and 'quota.share'
             main.ctrl = oneMainC, # A list with one element to define the start and end of 
             # the simulation.
             biols.ctrl = oneBioC,  # A list with one element to select the model to simulate 
             # the stock dynamics.
             fleets.ctrl = oneFlC,   # A list with several elements to select fleet dynamic models 
             # and store additional parameters.
             covars.ctrl = oneCvC,   # covars control not used 
             obs.ctrl = oneObsC,  # A list with one element to define how the stock observed 
             # ("PerfectObs").
             assess.ctrl = oneAssC,  # A list with one element to define how the stock assessment 
             # model used ("NoAssessment").
             advice.ctrl = oneAdvC)  # A list with one element to define how the TAC advice is 
# obtained ("IcesHCR").

# Names of the object returned by FLBEIA
names(s0)

# The default plot for FLBiol defined in FLCore
plot(s0$biols[[1]])

# Extract reference points for summaries
s0_brps <- extractBRP(oneAdvC, stkn = names(oneBio)) 

# Create summary data frames (biological, economic, and catch)
proj.yr     <- 2013 
s0_sum      <- bioSum(s0, brp = s0_brps)
s0$fleets$fl1 <- setUnitsNA(s0$fleets$fl1) # set units to NA to avoid errors in fltSum
s0_flt      <- fltSum(s0)
s0_fltStk   <- fltStkSum(s0)


# Create several plots and save them in the working directory using 'pdf' format and 
# 's0' suffix in the name.


plotFLBiols(s0$biols, pdfnm='s0', ss = 'all')
plotFLFleets(s0$fleets, pdfnm='s0', ss = 'all')
plotEco(s0, pdfnm='s0')
plotfltStkSum(s0, pdfnm='s0')

#------------------------------------------------------------ 
# Example with several iterations: oneIters
#------------------------------------------------------------

# Load the same data set as before but with 3 iterations.
# Run FLBEIA and plot the results

data(oneIt)

s1 <- FLBEIA(biols = oneItBio,   # FLBiols object with one FLBiol element for stk1.
             SRs = oneItSR,    # A list with one FLSRsim object for stk1.
             BDs = NULL,       # No Biomass Dynamic populations in this case.
             fleets = oneItFl,    # FLFleets object with on fleet.
             covars = oneItCv,    # covars not used
             indices = NULL,       # indices not used 
             advice = oneItAdv,   # A list with two elements 'TAC' and 'quota.share'
             main.ctrl = oneItMainC, # A list with one element to define the start and end of 
             # the simulation.
             biols.ctrl = oneItBioC,  # A list with one element to select the model to simulate 
             # the stock dynamics.
             fleets.ctrl = oneItFlC,   # A list with several elements to select fleet dynamic 
             # models and store additional parameters.
             covars.ctrl = oneItCvC,   # covars control not used 
             obs.ctrl = oneItObsC,  # A list with one element to define how the stock observed 
             # ("PerfectObs").
             assess.ctrl = oneItAssC,  # A list with one element to define how the stock 
             # assessment model used ("NoAssessment").
             advice.ctrl = oneItAdvC)  # A list with one element to define how the TAC advice is 
# obtained ("IcesHCR").

# Names of the object returned by FLBEIA
names(s1)

# The default plot for FLBiol defined in FLCore
plot(s1$biols[[1]])

# Extract reference points for summaries
s1_brps <- extractBRP(oneItAdvC, stkn = names(oneItBio)) 

# Create summary data frames (biological, economic, and catch)
proj.yr     <- 2013 
s1_bio     <- bioSum(s1, brp = s1_brps)
s1$fleets$fl1 <- setUnitsNA(s1$fleets$fl1) # set units to NA to avoid errors in fltSum
s1_flt     <- fltSum(s1)
s1_fltStk  <- fltStkSum(s1)

s1_bioQ    <- bioSumQ(s1_bio)
s1_fltQ    <- fltSumQ(s1_flt)
s1_fltStkQ <- fltStkSumQ(s1_fltStk)

s1b_bio     <- bioSum(s1, long = FALSE)
s1b_flt     <- fltSum(s1, long = FALSE)
s1b_fltStk  <- fltStkSum(s1, long = FALSE)

s1b_fltQ    <- bioSumQ(s1b_bio)
s1b_fltQ    <- fltSumQ(s1b_flt)
s1b_fltStkQ <- fltStkSumQ(s1b_fltStk)

# Create several plots and save them in the working directory using 'pdf' format and 
# 's1' suffix in the name.

plotFLBiols(s1$biols, pdfnm='s1', ss = 'all')
plotFLFleets(s1$fleets, pdfnm='s1', ss = 'all')
plotEco(s1, pdfnm='s1')
plotfltStkSum(s1, pdfnm='s1') 


#------------------------------------------------------------------ 
# Example with 2 stock, 2 Fleets, 4 seasons and 1 iteration: multi
#------------------------------------------------------------------

# Load the multi data set. This dataset has 2 stocks, one stk1 is 
# age structured and the second one stk2 is aggregated in biomass.

data(multi)

# Run FLBEIA.

s2 <- FLBEIA(biols = multiBio,   # FLBiols object with 2 FLBiol element for stk1.
             SRs = multiSR,    # A list with 1 FLSRsim object for stk1.
             BDs = multiBD,    # A list with 1 FLBDSim object for stk2.
             fleets = multiFl,    # FLFleets object with on fleet.
             covars = multiCv,    # covars not used
             indices = NULL,       # indices not used 
             advice = multiAdv,   # A list with two elements 'TAC' and 'quota.share'
             main.ctrl = multiMainC, # A list with one element to define the start and end 
             # of the simulation.
             biols.ctrl = multiBioC,  # A list with one element to select the model to simulate 
             # the stock dynamics.
             fleets.ctrl = multiFlC,   # A list with several elements to select fleet dynamic 
             # models and store additional parameters.
             covars.ctrl = multiCvC,   # covars control not used 
             obs.ctrl = multiObsC,  # A list with one element to define how the stock observed 
             # ("PerfectObs").
             assess.ctrl = multiAssC,  # A list with one element to define how the stock 
             # assessment model used ("NoAssessment").
             advice.ctrl = multiAdvC)  # A list with one element to define how the TAC advice is 
# obtained ("IcesHCR").

# Names of the object returned by FLBEIA
names(s2)

# The default plot for FLBiol defined in FLCore
plot(s2$biols[[1]])

# Extract reference points for summaries
s2_brps <- extractBRP(multiAdvC, stkn = names(multiBio)) 

# Create summary data frames (biological, economic, and catch)

s2_sum      <- bioSum(s2, brp = s2_brps)
for (fl in names(s2$fleets))  # set units to NA to avoid errors in fltSum
  s2$fleets[[fl]] <- setUnitsNA(s2$fleets[[fl]])
s2_flt      <- fltSum(s2)

s2b_flt     <- fltSum(s2, byyear = FALSE)

s2_fltStk   <- fltStkSum(s2)

# Create several plots and save them in the working directory using 'pdf' format and 
# 's2' suffix in the name.

plotFLBiols(s2$biols, pdfnm='s2', ss = 2)
plotFLFleets(s2$fleets, pdfnm='s2', ss = 2)
plotEco(s2, pdfnm='s2')
plotfltStkSum(s2, pdfnm='s2')

## End(Not run)

















