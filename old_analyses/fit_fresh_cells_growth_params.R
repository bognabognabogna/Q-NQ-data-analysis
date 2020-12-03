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

# Here data is taken from the global environment
data$predicted =  simulate_regrowth_with_lag2(a, Vh, Kh,freshLag, N0, Nprop=1, H0, data)
data_gathered = data %>%
  rename(real = biomass) %>%
  tidyr::gather(key="Biomass", value="biomass", real, predicted)
p1=ggplot(data_gathered, aes(time_hours,biomass, color = Biomass)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
    theme_bw() +
    ylab("biomass [mg/L]") +
    xlab("time [h]")
print(p1)
ggsave(filename = paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/manuscript/figs/", "Parameters_fit",  ".jpg"), plot = p1, width = 20, height = 10, units = "cm")