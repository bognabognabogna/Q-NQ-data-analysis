# set directory o the directory of this script!


suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(lubridate)))
suppressWarnings(suppressPackageStartupMessages(library(bbmle)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim)))
source('Q-NQ-helper-functions.R')                 


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


deoptim_fit_params = data.frame(N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0), bestval = numeric(0))
#Vh=106;Kh=106;a=0.0084 % fitted to freshcells when lag = 0
#Vh=165.63; Kh=127.81; a=0.0079; # fitted to fresh cells by Matlab when lag = 2.5
minLag = 0.5
types = c("Q", "NQ", "S")
for (type in types) {
#timeslag = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/LAG_Bogna.txt"), header = TRUE);
  for (i in 1:6) {
    max_time = ifelse(type == "NQ" & i > 2, 24, 16)
    dataQ = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/W", i, "_",type, " H2O_summary.txt"), header = TRUE);
    dataQ = dataQ %>%
      filter(time_hours <= max_time) %>%
      mutate(biomass = Biomass/mg_count)
    N0 = dataQ$biomass[1]
    data = dataQ
    deoptim_out = DEoptim(fn = sumLeastSquaresFitNprop3, lower = c(0,minLag), upper=c(1,10))
    deoptim_fit_params = deoptim_fit_params %>% 
      rbind(data.frame(N0= deoptim_out$optim$bestmem[1], lag =  deoptim_out$optim$bestmem[2], week = i, type = type, bestval = deoptim_out$optim$bestval))
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
  ylim(c(0,1)) +
  ylab('fitted N0')


ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = lag, col = type)) + 
  geom_line(aes(y = lag, col = type)) +
  ylim(c(0,10)) +
  ylab('fitted lag')

ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = bestval, col = type)) + 
  geom_line(aes(y = bestval, col = type)) +
  ylim(c(0,10)) +
  ylab('sum of squared errors for the best fit')


# DeOPtim gives same results like multistart but is quicker and doesn't require setting arbitrary start points

# Now do sensitivity analysis try to fit params for slightly shifted lags
