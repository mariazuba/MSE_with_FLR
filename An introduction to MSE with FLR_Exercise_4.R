# 4. In the observation model add a log-normal error with median 1 and CV = 30% in the observation of the index.
# Loads all necessary packages
library(FLa4a)
library(FLash)
library(FLXSA)
library(FLBRP)
library(ggplotFL)

#---------------------------------------------------------------------------------
## CONDITIONING THE OPERATING MODEL
#---------------------------------------------------------------------------------

# Read in stock assessment data
data(ple4)
data(ple4.index)
stk <- ple4 
idx <- FLIndices(idx=ple4.index)

# Set up the iteration and projection window parameters
it <- 20                    # iterations
y0 <- range(stk)["minyear"] # initial data year
dy <- range(stk)["maxyear"] # final data year
iy <- dy+1                  # initial year of projection (also intermediate year)
ny <- 12    # number of years to project from initial year
fy <- dy+ny # final year
nsqy <- 3   # number of years to compute status quo metrics


## ** Fit stock assessment model a4a **

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

# Reduce to keep one iteration only for reference points
stk0 <- qapply(stk, iterMedians)


# ** Fit the stock-recruit model **
srbh  <- fmle(as.FLSR(stk,  model="bevholt"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
srbh0 <- fmle(as.FLSR(stk0, model="bevholt"), method="L-BFGS-B", lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf)) 

# Generate stock-recruit residuals for the projection period
srbh.res <- rnorm(it, FLQuant(0, dimnames=list(year=iy:fy)), mean(c(apply(residuals(srbh), 6, sd))))

                                                                                    
## ** Calculate reference points **                                                                                   bh), 6, sd))))

# Calculate the reference points
brp  <- brp(FLBRP(stk0, srbh0))

Fmsy <- c(refpts(brp)["msy","harvest"])
msy  <- c(refpts(brp)["msy","yield"])
Bmsy <- c(refpts(brp)["msy","ssb"])
Bpa  <- 0.5*Bmsy
Blim <- Bpa/1.4

#** set up the operating model for the projection window **
#Prepare the FLStock object for projections
stk <- stf(stk, fy-dy, nsqy, nsqy)


#---------------------------------------------------------------------------------
## SET UP OBSERVATION ERROR MODEL ELEMENTS
#---------------------------------------------------------------------------------

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


#---------------------------------------------------------------------------------
# SET UP MSE LOOP
#---------------------------------------------------------------------------------

# Needed Functions:
# * Observation error model
# * XSA assessment model
# * Control object for projections

# * Observation error model
o <- function(stk, idx, assessmentYear, dataYears, obs.error.id) {
  # dataYears is a position vector, not the years themselves
  stk.tmp <- stk[, dataYears]
  # add small amount to avoid zeros
  catch.n(stk.tmp) <- catch.n(stk.tmp) + 0.1
  # Generate the indices - just data years
  idx.tmp <- lapply(idx, function(x) x[,dataYears])
  # Generate observed index
  for (i in 1:length(idx)) index(idx[[i]])[, assessmentYear] <-  stock.n(stk)[, assessmentYear]*
                                                                 index.q(idx[[i]])[, assessmentYear]*
                                                                 obs.error.id[, assessmentYear]
  return(list(stk=stk.tmp, idx=idx.tmp, idx.om=idx))
}

# * XSA assessment model
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

# * Control object for projections
getCtrl <- function(values, quantity, years, it){
  dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- unlist(values)
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA)) 
  ctrl@trgtArray <- arr0
  ctrl
}

#---------------------------------------------------------------------------------
# MSE initialisation
#---------------------------------------------------------------------------------

vy <- ac(iy:fy)
TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(dy,vy), iter=1:it))
TAC[,ac(dy)] <- catch(stk)[,ac(dy)]
TAC[,ac(iy)] <- TAC[,ac(dy)] #assume same TAC in the first intermediate year
ctrl <- getCtrl(c(TAC[,ac(iy)]), "catch", iy, it)

# Set up the operating model FLStock object
stk.om.obs.err.id <- fwd(stk, control=ctrl, sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE)

# * Observation error model
obs.error.id <- FLQuant(rnorm(prod(dim(stk.om.obs.err.id@m)),1,0.3), dimnames = dimnames(stk.om.obs.err.id@m))
obs.error.id[obs.error.id <= 0.1] <- 0.1

#---------------------------------------------------------------------------------
# Start the MSE loop
#---------------------------------------------------------------------------------
set.seed(231) # set seed to ensure comparability between different runs

for(i in vy[-length(vy)]){
  
  # set up simulations parameters
  ay <- an(i)
  cat(i, " > ")
  flush.console()
  vy0 <- 1:(ay-y0)              # data years (positions vector) - one less than current year
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector) - one less than current year
  
  # apply observation error
  oem    <- o(stk.om.obs.err.id, idx, i, vy0, obs.error.id) 
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
  stk.om.obs.err.id <- fwd(stk.om.obs.err.id, control=ctrl,sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE) 
}

#---------------------------------------------------------------------------------
# PERFORMANCE STATISTICS
#---------------------------------------------------------------------------------
# Some example performance statstics, but first isolate the projection period
stk.tmp <- window(stk.om.obs.err.id,start=iy)

# annual probability of being below Blim 
risky <- iterSums((ssb(stk.tmp)/Blim)<1)/it

# mean probability of being below Blim in the first half of the projection period
mean(risky[,1:trunc(length(risky)/2)])

# ...and second half 
mean(risky[,(trunc(length(risky)/2)+1):length(risky)])

# plot of SSB relative to Bmsy
boxplot(data~year,data=as.data.frame(ssb(stk.tmp)),main="SSB")
abline(h=Bmsy,col="red")

plot(FLStocks(stk.om.obs.err.id=stk.om.obs.err.id, stk.mp=stk.mp)) + theme(legend.position="top") + geom_vline(aes(xintercept=2009))


save(stk.om.obs.err.id, file = 'C:/use/Dropbox/FLBEIA Course/Tutorials/02_MSE_with_FLR/stk_exe4.RData')

#---------------------------------------------------------------------------------
# EXERCISES
#---------------------------------------------------------------------------------

# 1. Change the stock-recruitment relationship in the OM for a ricker model.
# 2. In the OM change the natural mortality, make it variable in the range [0.05,0.2]
# 3. Witn natural mortality variable in the OM, in the observation of the stock use M = 0.1
# 4. In the observation model add a log-normal error with median 1 and CV = 30% in the observation of the index.
# 5. In the implementation of the advice add an overshoot of the TAC using a uniform distribution with range [1,1.3] 



























