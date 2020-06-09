suppressWarnings(suppressPackageStartupMessages(library(deSolve)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(lubridate)))
suppressWarnings(suppressPackageStartupMessages(library(bbmle)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(DEoptim)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
source('~/Q-NQ-data-analysis/Q-NQ-helper-functions.R')    
deopticontrol = DEoptim.control(itermax = 1000, reltol = 10^(-8), trace = 100)

# Here we try to fit the initial proporion of Q in S.
# Assuming  NQ and Q have lags and initial proportions as fitted to the data

# molar mass of glucose: 180g / mol = 0.180 g per mmol
# 2% glucose i.e. 20g per 1L = 20[g/L] * 1/0.18 [mmol/g] = 20/0.18 [mmol/L]
H0 = 111
# we will rescale number of cells to miligrams so that we have a reasnable order of magnitude
# assume 1mg = 33*10^6 (MacLean et.al. 2010)
mg_count=33*10^6
max_time = 24
# get parsms fitted to fresh cells
a=readRDS("~/Q-NQ-data-analysis/a.Rds")
Kh=readRDS("~/Q-NQ-data-analysis/Kh.Rds")
Vh=readRDS("~/Q-NQ-data-analysis/Vh.Rds")
deoptim_fit_params = readRDS("~/Q-NQ-data-analysis/deoptim_fit_params_and_NQ2")
weeks = 1:6

selected_type = "S"
prop_real = 0.75

all_prop_fited = data.frame(propQ = numeric(0),
                        week = numeric(0), 
                        bestval = numeric(0))

for (selected_week in weeks) {
#data = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20/W", selected_week, "_",selected_type, " H2O_summary.txt"), header = TRUE);
data_raw = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
data_raw = data_raw %>% rename(biomass = meanBiomass)
data_mean = GetMeanValues(data_raw,c("time_hours")) %>%
  mutate(biomass = meanBiomass/mg_count)

data = data_mean %>%
  filter(time_hours <= max_time) 
N0 = data$biomass[1]


NpropQ =  deoptim_fit_params %>% filter(week == selected_week & type == "Q") %>% pull(N0)
tlagQ = deoptim_fit_params %>% filter(week == selected_week & type == "Q") %>% pull(lag)
NpropNQ =  deoptim_fit_params %>% filter(week == selected_week & type == "NQ") %>% pull(N0)
tlagNQ = deoptim_fit_params %>% filter(week == selected_week & type == "NQ") %>% pull(lag)


optim_prop_out = DEoptim(fn = sumLeastSquaresFitPropQAndNQ, 
                      lower = c(0), upper=c(1),
                      control = deopticontrol)

prop_fited = data.frame(propQ = optim_prop_out$optim$bestmem[1],
                        week = selected_week, 
                        bestval = optim_prop_out$optim$bestval)
                        
all_prop_fited = rbind(all_prop_fited,prop_fited)
data_new = data




# Here data is taken from the global environment
data_new$fitted = simulate_regrowth_two_srains_with_lag(a, Vh, Kh, tlagQ,tlagNQ, N0, NpropQ, NpropNQ, prop_fited$propQ, H0, data$time_hours)
data_new$predicted = simulate_regrowth_two_srains_with_lag(a, Vh, Kh, tlagQ,tlagNQ, N0, NpropQ, NpropNQ, prop_real, H0, data$time_hours)
data_new = data_new %>%
  tidyr::gather(key="selected_type", value="total_biomass", biomass, predicted, fitted)
p1=ggplot(data_new) + 
  geom_point(aes(time_hours,total_biomass, color = selected_type)) +
  geom_line(aes(time_hours,total_biomass, color = selected_type)) +
  ggtitle(paste0("Pop: ", selected_type, " week: ", selected_week)) +
  theme_bw()
print(p1)
}


