
suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(lubridate)))
suppressWarnings(suppressPackageStartupMessages(library(bbmle)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim))
                 
# our model of growth (using differential euqutions)
simulate_regrowth_single_strain = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
  Jg =Vh*G/(Kh+G) 
  Gdot   = -Jg*N 
  Ndot   = a*Jg*N
return(list(c(Gdot, Ndot)))
  })
}

## this is an example of how to run an ODE with defined parameters and initial values
#pars = c(a   = 0.008,    # /day, rate of ingestion
#         Vh  = 126.0,    # /day, growth rate of prey
#         Kh  = 90.0)     # mmol/m3, carrying capacity
#
#yini  <- c(G = 111, N = 0.001)
#times <- seq(0, 20, by = 1)
#out   <- ode(yini, times, simulate_regrowth_single_strain, pars)
#plot(out[, "time"], out[, "N"])


# FITTING PARAMETERS TO FRESH CELLS
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
simulate_regrowth_with_lag = function(a, Vh, Kh, tlag, N0, Nprop, H0) {
  timesLag = data %>% filter(time_hours < tlag) %>% pull(time_hours)
  timesGrowth = data %>% filter(time_hours >= tlag) %>% mutate(time_hours = time_hours - tlag) %>% pull(time_hours)
  aliveN0 = Nprop*N0
  deadN0 = (1-Nprop)*N0
  pars <- c(a = a, Vh = Vh, Kh = Kh)
  inits = c(G=H0, N=aliveN0)
  simulation <- as.data.frame(ode(inits, timesGrowth, simulate_regrowth_single_strain, pars)) %>% mutate(N = N + deadN0)
  simulatedN <- c(rep(N0, length(timesLag)),  simulation$N)
  return(simulatedN)
}

# we will use that function to fit the metabolic parameters Vh and Kh which are assumed to be the same for fresh and starved cells. Here we assume tlag is known
# it returns sum of squared errors so we want to find Kh and Vh such that they minimise this function
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
sumLeastSquaresFitGrowth = function(param, a, tlag, N0, H0) {
  Vh = param[1]
  Kh = param[2]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop=1, H0)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}



# we will use that function to fit the metabolic parameters Vh and Kh which are assumed to be the same for fresh and starved cells and parameter tlag
# it returns sum of squared errors so we want to find Kh, Vh and tlag such that they minimise this function
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
#Parameters a,N0 and H0 are taken from the global environment
sumLeastSquaresFitGrowthToDeoptim = function(param) {
  Vh = param[1]
  Kh = param[2]
  tlag = param[3]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop=1, H0)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}

# we will use that function to fit No assuming lag is already known,
# it returns sum of squared errors so we want to find N0 such that they minimise this function
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
sumLeastSquaresFitNprop = function(param, a, Vh,Kh, N0, H0, tlag) {
  Nprop = param[1]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}

# we will use that function to fit No and lag at the same time.
# it returns sum of squared errors so we want to find N0 and lag such that they minimise this function
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
sumLeastSquaresFitNprop2 = function(param, a, Vh,Kh, N0, H0) {
  Nprop = param[1]
  tlag = param[2]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}

# Same as above but a, Vh, Kh, tlag,N0 and H0 are taken from the global environment
# This is because function DEoptim can only handle a function of one input = param
sumLeastSquaresFitNprop3 = function(param) {
  Nprop = param[1]
  tlag = param[2]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}


# This is a function used for MLE fitting
#nll = function(a, Vh, Kh, mu, sigma, tlag, N0, Nprop, H0) {
#  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
#  difference = data$biomass - simulatedN
#  logf = suppressWarnings(dnorm(difference[-1],mu,sigma,log = TRUE))
#  ll = sum(logf)
#  return(-ll)
#}


# molar mass of glucose: 180g / mol = 0.180 g per mmol
# 2% glucose i.e. 20g per 1L = 20[g/L] * 1/0.18 [mmol/g] = 20/0.18 [mmol/L]
H0 = 111
# we will rescale number of cells to miligrams so that we have a reasnable order of magnitude
# assume 1mg = 33*10^6 (MacLean et.al. 2010)
mg_count=33*10^6
max_time = 16
freshCells = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/w0_fresh_cells/W0_fresh cells-summary.txt", header = TRUE, row.names = NULL)
freshCells = freshCells %>% 
  filter(time_hours <= max_time) %>%
  mutate(biomass = meanBiomass/mg_count)
# find initial biomass (dead + alive)
N0 = freshCells$biomass[1]
# find final biomass
Nend = freshCells$biomass[nrow(freshCells)]
# find the parameter a i.e. how many grams of proteins can be created per mmol of glucose
a = mean((Nend-N0)/(H0))


## This is an example of how to fit an MLE estimator
#free = c(Vh=50, Kh=50, sigma=1)
#fixed = c(a=a, mu = 0, tlag = 2, N0=N0, Nprop=1, H0=H0)
## Use MLE
#data = freshCells
#fit <- mle2(nll, start = as.list(free), fixed = as.list(fixed), method = "Nelder-Mead")
# Use optim
#fit2 = optim(par = c(as.numeric(free['Vh']), as.numeric(free['Kh'])), fn = sumLeastSquaresFitGrowth,  a=a, tlag=2, N0=N0, H0=H0,
#            method = "L-BFGS-B") #, lower = c(0,0), upper=c(500,500)
#fitted_params = coef(fit)
#freshCells$predicted =simulate_regrowth_with_lag(a=as.numeric(fitted_params['a']), 
#                                                 Vh=as.numeric(fitted_params['Vh']), 
#                                                 Kh=as.numeric(fitted_params['Kh']), 
#                                                 tlag=as.numeric(fixed['tlag']), 
#                                                 N0=as.numeric(fixed['N0']),
#                                                 Nprop=as.numeric(fitted_params['Nprop']),
#                                                 H0=as.numeric(fitted_params['H0']))
#fit.ci <- confint(fit)
#
#out = freshCells %>%
#  tidyr::gather(key="type", value="biomass", biomass, predicted)
#ggplot(out) + 
#  geom_point(aes(time_hours,biomass, color = type))


# fit params to fresh cells using Deoptim
# restrict the growth curve to 16h onl, otherwise we get unreliable parameters. Probably because there is a lot of noise between h12 and h24
data = freshCells
deoptim_out = DEoptim(fn = sumLeastSquaresFitGrowthToDeoptim, lower = c(0,0,0), upper=c(300,300,3))

# Q cells
Vh=deoptim_out$optim$bestmem[1] %>% as.numeric()#fitted_params['Vh']
Kh=deoptim_out$optim$bestmem[2] %>% as.numeric()#fitted_params['VK']
freshLag = deoptim_out$optim$bestmem[3] %>% as.numeric()

#single_fit_params = data.frame(N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0))
#multi_fit_params = data.frame(N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0))
deoptim_fit_params = data.frame(N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0))
#Vh=106;Kh=106;a=0.0084 % fitted to freshcells when lag = 0
#Vh=165.63; Kh=127.81; a=0.0079; # fitted to fresh cells by Matlab when lag = 2.5
minLag = 0.5
H0 = 111
# we will rescale number of cells to miligrams so that we have a reasnable order of magnitude
# assume 1mg = 33*10^6 (MacLean et.al. 2010)
mg_count=33*10^6
#pm <- expand.grid(seq(0.1,0.9,0.1), seq(0.1,10.1,1))
#pm <- as.matrix(pm)
#colnames(pm) = NULL
types = c("Q", "NQ", "S")
for (type in types) {
#timeslag = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/LAG_Bogna.txt"), header = TRUE);
max_time = ifelse(type == "NQ" & i > 2, 24, 16)
for (i in 1:6) {
dataQ = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/W", i, "_",type, " H2O_summary.txt"), header = TRUE);
dataQ = dataQ %>%
  filter(time_hours <= max_time) %>%
  mutate(biomass = Biomass/mg_count)
N0 = dataQ$biomass[1]
data = dataQ
#Nprop0=0.5;
#tlag0=2;
#fit_single = optim(par = c(Nprop0, tlag0), 
#            fn = sumLeastSquaresFitNprop2,  a=a, Vh=Vh, Kh = Kh,N0 = dataQ$biomass[1],H0=H0,
#            method = "L-BFGS-B", lower = c(0,0), upper=c(1,10))
#single_fit_params = single_fit_params %>% rbind(data.frame(N0= fit_single$par[1], lag =  fit_single$par[2], week = i, type = type))


#fit_multi = optimr::multistart(par =pm, , 
#                               fn = sumLeastSquaresFitNprop3, method = "L-BFGS-B", lower = c(0.01,0.01), upper=c(0.99,10))
#bsst_fit=fit_multi[fit_multi$value == min(fit_multi$value),]
#multi_fit_params = multi_fit_params %>% rbind(data.frame(N0= bsst_fit$p1, lag =  bsst_fit$p2, week = i, type = type))


deoptim_out = DEoptim(fn = sumLeastSquaresFitNprop3, lower = c(0,minLag), upper=c(1,10))
deoptim_fit_params = deoptim_fit_params %>% rbind(data.frame(N0= deoptim_out$optim$bestmem[1], lag =  deoptim_out$optim$bestmem[2], week = i, type = type))
dataQ_new =dataQ
dataQ_new$predicted = simulate_regrowth_with_lag(a, Vh, Kh,deoptim_out$optim$bestmem[2] %>% as.numeric(),N0, deoptim_out$optim$bestmem[1] %>% as.numeric(), H0)
dataQ_new = dataQ_new %>%
  tidyr::gather(key="type", value="biomass", biomass, predicted)
p1=ggplot(dataQ_new) + 
  geom_point(aes(time_hours,biomass, color = type)) +
  ggtitle(paste0("Pop: ", type, " week: ", i))
print(p1)
}
}

ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = N0, col = type)) + 
  geom_line(aes(y = N0, col = type)) +
  ylim(c(0,1))


ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = lag, col = type)) + 
  geom_line(aes(y = lag, col = type)) +
  ylim(c(0,10))


# DeOPtim gives same results like multistart but is quicker and doesn't require setting arbitrary start points
