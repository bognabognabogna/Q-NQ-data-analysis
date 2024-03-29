# set directory o the directory of this script!


# UWAGA: MY ZNAJDUJEMY N0 jako %poczatkowej biomassy widzianej w naszym OD.
# TO NIE JEST %glodzonych komorek, ktore wracaja do podzialu
# bo czesc komorek glodzonych mogla sie dawno rozpaść i nie jest widoczna w poczatkowym OD
# jeśli chcemy %Głodzonych komórek które wracają do podziału
# To musimy wziąć N0*initial biomass / biomass at the beginning of starvation (density)

suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(lubridate)))
suppressWarnings(suppressPackageStartupMessages(library(bbmle)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
source('~/Q-NQ-data-analysis/Q-NQ-helper-functions.R')               
deopticontrol = DEoptim.control(itermax = 1000, reltol = 10^(-8), trace = 100)

# molar mass of glucose: 180g / mol = 0.180 g per mmol
# 2% glucose i.e. 20g per 1L = 20[g/L] * 1/0.18 [mmol/g] = 20/0.18 [mmol/L]
H0 = 111
# we will rescale number of cells to miligrams so that we have a reasnable order of magnitude
# assume 1mg = 33*10^6 (MacLean et.al. 2010)
mg_count=33*10^6
max_time = 16
#freshCells = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/w0_fresh_cells/W0_fresh cells-summary.txt", header = TRUE, row.names = NULL)
freshCells = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/w0_fresh_cells/W0_fresh cells.txt", header = TRUE, row.names = NULL)

cols_to_group_by = c("time_hours")
freshCells = GetMeanValues(freshCells,cols_to_group_by)
freshCells = freshCells %>% 
  filter(time_hours <= max_time) %>%
  mutate(biomass = meanBiomass/mg_count) %>%
  select(time_hours, biomass)
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
#Parameters data, a,N0 and H0 are taken from the global environment
deoptim_out = DEoptim(fn = sumLeastSquaresFitGrowthToDeoptim, lower = c(0,0,0), upper=c(300,300,3),
                      control = deopticontrol)


Vh=deoptim_out$optim$bestmem[1] %>% as.numeric()
Kh=deoptim_out$optim$bestmem[2] %>% as.numeric()
freshLag = deoptim_out$optim$bestmem[3] %>% as.numeric()

saveRDS(a, "~/Q-NQ-data-analysis/a.Rds")
saveRDS(Kh, "~/Q-NQ-data-analysis/Kh.Rds")
saveRDS(Vh, "~/Q-NQ-data-analysis/Vh.Rds")



deoptim_fit_params = data.frame(N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0), bestval = numeric(0))
#Vh=106;Kh=106;a=0.0084 % fitted to freshcells when lag = 0
#Vh=165.63; Kh=127.81; a=0.0079; # fitted to fresh cells by Matlab when lag = 2.5
minLag = 0.5
types = c("Q", "NQ", "S")
weeks = 1:6


for (selected_type in types) {
osberved_lags = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/LAG_Bogna.txt"), header = TRUE)
#timeslag = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/LAG_Bogna.txt"), header = TRUE);
  for (selected_week in weeks) {
    max_time = 12
    #max_time = ifelse(selected_type == "NQ" & selected_week > 2, 24, 16)
    maxLag = osberved_lags %>% filter(week == selected_week) %>% pull(meanLAG)
    #data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/W", selected_week, "_",selected_type, " H2O_summary.txt"), header = TRUE);
    data_raw = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
    data_raw = data_raw %>% rename(biomass = meanBiomass)
    data_mean = GetMeanValues(data_raw,c("time_hours")) %>%
      mutate(biomass = meanBiomass/mg_count)
    
    data = data_mean %>%
      filter(time_hours <= max_time) 
    N0 = data$biomass[1]

    # Parameters a,Vh,Kh,N0,H0, data are taken from the global env
    deoptim_out = DEoptim(fn = sumLeastSquaresFitNprop3, lower = c(0,minLag), upper=c(1,maxLag),
                          control = deopticontrol)
    
    this_deoptim_fit_params = data.frame(N0= deoptim_out$optim$bestmem[1], 
                                          lag =  deoptim_out$optim$bestmem[2],  
                                          week = selected_week, 
                                          type = selected_type, 
                                          bestval = deoptim_out$optim$bestval)
    # check
    print(sumLeastSquaresFitNprop3(c(this_deoptim_fit_params$N0, this_deoptim_fit_params$lag)) == this_deoptim_fit_params$bestval)
    
    deoptim_fit_params = deoptim_fit_params %>% 
      rbind(this_deoptim_fit_params)
    
    
    data_new = data_mean
    # Here data is taken from the global environment
    data_mean$predicted = simulate_regrowth_with_lag(a, Vh, Kh,this_deoptim_fit_params$lag,N0, Nprop=this_deoptim_fit_params$N0, H0, data_mean)
    data_mean = data_mean %>%
      tidyr::gather(key="selected_type", value="biomass", biomass, predicted)
    p1=ggplot(data_mean) + 
      geom_point(aes(time_hours,biomass, color = selected_type)) +
      ggtitle(paste0("Pop: ", selected_type, " week: ", selected_week)) +
      theme_bw()
    print(p1)
    
    deoptim_out = NULL
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
  ylab('mean squared error for the best fit')

saveRDS(deoptim_fit_params , "~/Q-NQ-data-analysis/deoptim_fit_params_and_NQ2")


# DeOPtim gives same results like multistart but is quicker and doesn't require setting arbitrary start points
eps = 0.0001
#eps = 0.0004 # on average 0.02 difference per point (0.0024 is 0.01 OD)
# 5% from the minimal mse
mse_landscape = data.frame(N0=numeric(0), lag=numeric(0), week=numeric(0), type = numeric(0), mse = numeric(0))
for (selected_type in types) {
  osberved_lags = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/LAG_Bogna.txt"), header = TRUE)
for (selected_week in weeks) {
  maxLag = osberved_lags %>% filter(week == selected_week) %>% pull(meanLAG)
  minLag = 0.5
  #max_time = ifelse(selected_type == "NQ" & selected_week > 2, 24, 16)
  max_time = 12
  this_deoptim_fit_params = deoptim_fit_params %>% filter(week == selected_week & type == selected_type)
  data_raw = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
  data_raw = data_raw %>% rename(biomass = meanBiomass)
  data_mean = GetMeanValues(data_raw,c("time_hours")) %>%
    mutate(biomass = meanBiomass/mg_count)
  data = data_mean %>%
    filter(time_hours <= max_time) 
  N0 = data$biomass[1]
  # data, a, Vh, Kh, tlag, N0, H0 are taken from the global env
  print(paste0(" a: ",a, " Vh: ", Vh, " Kh: ",Kh, " N0: ", N0, " H0: ", H0))
  
  # Note a tiny change in tlag may make a big difference in the MSE. As tlag = 0.50001 -> we assume N(0.5) = N0, tlag = 0.499999 -> N(0.5) > N0!
  #opt_mse =  sumLeastSquaresFitNpropFirst10Obs(c(this_deoptim_fit_params$N0, this_deoptim_fit_params$lag))
  opt_mse =  sumLeastSquaresFitNprop3(c(this_deoptim_fit_params$N0, this_deoptim_fit_params$lag))
  lags = seq(minLag, maxLag, 0.1)
  N0s = seq(0,1, 0.05)
  for (no in N0s) {
    for (lag in lags) {
      mse =  sumLeastSquaresFitNprop3(c(no, lag))
      #mse = sumLeastSquaresFitNpropFirst10Obs(c(no, lag))
      this_mse = data.frame(N0=no, lag=lag, week=selected_week, type = selected_type, mse = mse)
      mse_landscape = rbind(mse_landscape, this_mse)
    }}
  min_mse = mse_landscape %>% filter(type == selected_type & week == selected_week) %>% pull(mse) %>% min()
  p=ggplot(mse_landscape %>% 
             filter(type == selected_type & week == selected_week) %>% 
             #mutate(mse = ifelse(mse < min_mse*(1 + eps), mse, NA),
             mutate(mse = ifelse(mse < opt_mse + eps, mse, NA),
                    week = paste0("week ", week))) + 
    geom_tile(aes(y=N0, x=lag, fill=mse)) +
    geom_text(data = data.frame(x= this_deoptim_fit_params$lag, y = this_deoptim_fit_params$N0),
               aes(x=x,y=y), label = "*", col = "red") +
    ggtitle("MSE for specific pairs possible of parameters") +
    scale_fill_gradient(high="gray90", low="darkblue", na.value = "white") +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle=30,size = 6),
          axis.text.y = element_text(size = 6),
          panel.background = element_rect(fill = "white")) +
    facet_grid(week ~ type)
  print(p)
}
}

Min = function(x) {if (length(x)==0){NA}else{min(x)}}
Max = function(x) {if (length(x)==0){NA}else{max(x)}}
dataWithConfidence = mse_landscape %>%
  inner_join(deoptim_fit_params %>% select(N0optim = N0, lagoptim = lag, bestval, week, type), by = c("week", "type")) %>%
  rowwise() %>%
  mutate(goodFit = mse - bestval <= 0.0004) %>%
  #mutate(goodFit = mse <= 1.05*bestval) %>%
  group_by(week, type, N0optim) %>%
  #filter(goodFit) %>%
  summarise(minN0 = Min(N0[goodFit]),
            maxN0 = Max(N0[goodFit])) %>%
  ungroup() %>%
  mutate(minN0 = ifelse(is.na(minN0), N0optim, minN0 ),
         maxN0 = ifelse(is.na(maxN0), N0optim, maxN0 ))

ggplot(data = dataWithConfidence) +
  geom_line(aes(x=week, y = N0optim, col = type)) +
  geom_point(aes(x=week, y = N0optim, col = type)) +
  geom_pointrange(aes(x = week, y = N0optim, ymin = minN0, ymax = maxN0, col = type)) +
  ylab("Fitted N0")


# For now do not run the code below
if (FALSE) {
# Now do sensitivity analysis try to fit params for slightly shifted lags
# First vary lag by some coefficient constant for each week and type
epsilons = c(0.8, 1.5)
deoptim_fit_params_deviated = data.frame(epsilon = numeric(0), N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0), bestval = numeric(0))
for (selected_epsilon in epsilons)
  for (selected_type in types) {
    osberved_lags = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/LAG_Bogna.txt"), header = TRUE)
    #timeslag = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/LAG_Bogna.txt"), header = TRUE);
    for (selected_week in weeks) {
      optim_params = deoptim_fit_params %>% filter(week == selected_week & type == selected_type)
      maxLag = osberved_lags %>% filter(week == selected_week) %>% pull(meanLAG)
      deviated_lag = min(max(minLag, selected_epsilon*optim_params$lag), maxLag)
      max_time = ifelse(selected_type == "NQ" & selected_week > 2, 24, 16)
      
      #data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/W", selected_week, "_",selected_type, " H2O_summary.txt"), header = TRUE);
      data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
      data = data %>% rename(biomass = meanBiomass)
      data = GetMeanValues(data,c("time_hours"))
      data = data %>%
        filter(time_hours <= max_time) %>%
        mutate(biomass = meanBiomass/mg_count)
      N0 = data$biomass[1]
      tlag = deviated_lag
      # data, a, Vh, Kh, tlag, N0, H0 are taken from the global env
      print(paste0(" a: ",a, " Vh: ", Vh, " Kh: ",Kh, " N0: ", N0, " H0: ", H0, " tlag: ", tlag))
      deoptim_out = DEoptim(fn = sumLeastSquaresFitNprop, lower = c(0), upper=c(1),
                            control = deopticontrol)
      this_deoptim_fit_params_deviated = data.frame(N0= deoptim_out$optim$bestmem[1], 
                                                    lag =  tlag,  
                                                    week = selected_week, 
                                                    type = selected_type, 
                                                    bestval = deoptim_out$optim$bestval, 
                                                    epsilon = selected_epsilon)
      # check
      print(sumLeastSquaresFitNprop3(c(this_deoptim_fit_params_deviated$N0, this_deoptim_fit_params_deviated$lag)) == this_deoptim_fit_params_deviated$bestval)
      
      deoptim_fit_params_deviated = deoptim_fit_params_deviated %>% 
        rbind(this_deoptim_fit_params_deviated)
      deoptim_out = NULL
      
    }
  }



rownames(deoptim_fit_params_deviated) = NULL
deoptim_fit_params = deoptim_fit_params %>%
  mutate(epsilon = 1.0) 
deoptim_fit_params = rbind(deoptim_fit_params, deoptim_fit_params_deviated) %>% 
  mutate(epsilon = as.factor(epsilon))

ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = bestval, col = epsilon)) + 
  geom_line(aes(y = bestval, col = epsilon, group = epsilon)) +
  facet_wrap("type") +
  ylab(paste0('mean squared error for the best fit: when tlag = optminal lag*epsilon '))

#@#ggplot(deoptim_fit_params, aes(x = week)) +
#  geom_point(aes(y = N0, col = epsilon)) + 
#  geom_line(aes(y = N0, col = epsilon)) +
#  facet_wrap("type") +
#  ylab(paste0('fitted N0: when tlag = optminal lag*epsilon'))

# Find points with relativelu close bestval

ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = lag, col = epsilon)) + 
  geom_line(aes(y = lag, col = epsilon)) +
  facet_wrap("type") +
  ylab(paste0('tlag: when tlag = optminal lag*epsilon'))

ggplot(deoptim_fit_params, aes(x = week)) +
  geom_point(aes(y = N0, col = type)) + 
  geom_line(aes(y = N0, col = type)) +
  ylim(c(0,1)) +
  facet_wrap("epsilon") +
  ylab(paste0('fitted N0: when tlag = optminal lag*epsilon'))



# Now perfrom sensitivity analyss by adding random noise
set.seed(1)
stdv = 1
stochastic_runs = 1:10


N = length(stochastic_runs)*length(types)*length(weeks)
noise = runif(N, -sqrt(3)*stdv, sqrt(3)*stdv) #rnorm(N, 0, stdv)
deoptim_fit_params_rand_deviated = data.frame(epsilon = numeric(0), N0 = numeric(0), lag = numeric(0), week = integer(0), type = character(0), bestval = numeric(0))

noise_index = 0
for (stochastic_run in stochastic_runs)
  for (selected_type in types) {
    osberved_lags = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/LAG_Bogna.txt"), header = TRUE)
    #timeslag = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", type, "_H20/LAG_Bogna.txt"), header = TRUE);
    for (selected_week in weeks) {
      noise_index = noise_index +1
      optim_params = deoptim_fit_params %>% filter(week == selected_week & type == selected_type & epsilon == 1)
      maxLag = osberved_lags %>% filter(week == selected_week) %>% pull(meanLAG)
      deviated_lag = min(max(minLag, optim_params$lag + noise[noise_index]), maxLag)
      
      max_time = ifelse(selected_type == "NQ" & selected_week > 2, 24, 16)
      #data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/W", selected_week, "_",selected_type, " H2O_summary.txt"), header = TRUE);
      data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
      data = data %>% rename(biomass = meanBiomass)
      data = GetMeanValues(data,c("time_hours"))
      data = data %>%
        filter(time_hours <= max_time) %>%
        mutate(biomass = meanBiomass/mg_count)
      N0 = data$biomass[1]
      tlag = deviated_lag
      # data, a, Vh, Kh, tlag, N0, H0 are taken from the global env
      print(paste0(" a: ",a, " Vh: ", Vh, " Kh: ",Kh, " N0: ", N0, " H0: ", H0, " tlag: ", tlag))
      deoptim_out = DEoptim(fn = sumLeastSquaresFitNprop, lower = c(0), upper=c(1),
                            control = deopticontrol)
      this_deoptim_fit_params_rand_deviated = data.frame(N0= deoptim_out$optim$bestmem[1], 
                                                         lag =  tlag,  
                                                         week = selected_week, 
                                                         type = selected_type, 
                                                         bestval = deoptim_out$optim$bestval, 
                                                         epsilon = paste0("run", stochastic_run))
      # check
      #print(sumLeastSquaresFitNprop3(c(this_deoptim_fit_params_rand_deviated$N0, this_deoptim_fit_params_rand_deviated$lag)) == this_deoptim_fit_params_rand_deviated$bestval)
      
      deoptim_fit_params_rand_deviated = deoptim_fit_params_rand_deviated %>% 
        rbind(this_deoptim_fit_params_rand_deviated)
      deoptim_out = NULL
      
    }
  }



rownames(deoptim_fit_params_rand_deviated) = NULL
deoptim_fit_params = rbind(deoptim_fit_params %>% filter(epsilon == 1)  %>% mutate(run_type = "main"), 
                           deoptim_fit_params_rand_deviated  %>% mutate(run_type = "with noise"))  %>%
  mutate(group = paste0(run_type, "_", epsilon, "_", type))



#@#ggplot(deoptim_fit_params_rand_deviated, aes(x = week)) +
#  geom_point(aes(y = N0, col = epsilon)) + 
#  geom_line(aes(y = N0, col = epsilon)) +
#  facet_wrap("type") +
#  ylab(paste0('fitted N0: when tlag = optminal lag*epsilon'))

# OUR biomass are between 0 and 1. So eg et 0.5 they may differ by 0.01 -> MSE = 0.0001
real_bestval = deoptim_fit_params %>% filter(epsilon == 1) %>% rename(realbestval = bestval) %>% select(week, type, realbestval)
deoptim_fit_params_ok = deoptim_fit_params %>%
  inner_join(real_bestval, by = c("week", "type")) %>%
  mutate(bestval_diff = bestval - realbestval) %>%
  filter(bestval_diff < diff_threshold)

ggplot(deoptim_fit_params_ok) +
  geom_point(aes(x = week, y = bestval, col = run_type), size = 0.5) + 
  facet_wrap("type") +
  ylab(paste0('MSE for the best fit')) +
  theme_bw()

ggplot(deoptim_fit_params_ok, aes(x = week)) +
  geom_point(aes(y = lag, col = run_type)) + 
  geom_line(aes(y = lag, col = run_type, group = group)) +
  facet_wrap("type") +
  ylim(c(0, maxLag)) +
  ylab(paste0('tlag: when tlag = optminal lag with noise')) +
  theme_bw()

ggplot(deoptim_fit_params_ok, aes(x = week)) +
  geom_point(aes(y = N0, col = type)) + 
  geom_line(aes(y = N0, col = type, group = group)) +
  ylim(c(0,1)) +
  facet_wrap("run_type") +
  ylab(paste0('Fitted N_prop')) +
  theme_bw()

saveRDS(deoptim_fit_params, "/Users/bognasmug/Q-NQ-data-analysis/deoptim_fit_params_1000iter2.Rds")

aa=deoptim_fit_params_ok %>% 
  group_by(type, week) %>%
  arrange(bestval_diff) %>%
  summarise(N0 =  N0[1],
            lag = lag[1]) %>%
  ungroup() %>%
  mutate(epsilon = 0,
         type = as.character(type))

bb = deoptim_fit_params_ok %>% 
  filter(epsilon == 1) %>%
  mutate(epsilon = as.numeric(epsilon)) %>%
  select(type, week, N0,lag, epsilon) %>%
  rbind(aa)

ggplot(bb, aes(x = week)) +
  geom_point(aes(y = N0, col = type)) + 
  geom_line(aes(y = N0, col = type)) +
  ylim(c(0,1)) +
  ylab(paste0('best of Fitted N_prop')) +
  facet_wrap('epsilon') +
  theme_bw()  
}