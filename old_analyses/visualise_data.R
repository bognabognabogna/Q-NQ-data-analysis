source('~/Q-NQ-data-analysis/Q-NQ-helper-functions.R')
library(ggplot2)
library(dplyr)

mg_count=33*10^6
data = data.frame(time_hours = numeric(0), biomass = numeric(0), type=character(0), week = integer(0))
for (selected_week in 1:6) {
  for (selected_type in c("Q", "NQ", "S")) {
    data_raw = read.table(paste0("/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/experiment2/input/", selected_type, "_H20_all_no_correction/W_", selected_week, " ",selected_type, " H2O - not corrected.txt"), header = TRUE);
    data_raw = data_raw %>% rename(biomass = meanBiomass)
    data_mean = GetMeanValues(data_raw,c("time_hours")) %>%
      mutate(biomass = meanBiomass/mg_count,
             type = selected_type,
             week =  selected_week) %>%
      select(time_hours, biomass, type, week)
    data = rbind(data, data_mean) 
  }
  }


ggplot(data %>% filter(time_hours <= 24)) +
  geom_line(aes(x = time_hours, color = type, y = biomass)) +
  facet_wrap('week')

data_long = data %>% mutate(total_time = week*7*24 + time_hours)
ggplot(data = data_long, aes(x= total_time, y = biomass, color = type)) +
  geom_line() + geom_point(alpha = 0.5) +
  geom_point(data = data_long %>% filter(total_time == 7*24+0.5), 
             aes(x= total_time, y = biomass, color = type), size = 5) +
  theme_bw()


data1=data %>% filter(time_hours %in% c(0.5, 2, 5, 8, 12, 18, 24))

ggplot(data1) +
  geom_line(aes(x = week, color = type, y = biomass)) +
  geom_point(aes(x = week, color = type, y = biomass)) +
  ylab("initial biomass") +
  theme_bw() +
  facet_grid(.~time_hours)


datInit = data %>% filter(time_hours == 1) %>% select(week, biomass,type)
dataHypotheticalSInit = datInit %>% 
  group_by(week) %>% 
  summarise(biomass = 0.75*biomass[type == "Q"] + 0.25*biomass[type == "NQ"]) %>%
  ungroup() %>%
  mutate(type = "75%Q + 25%NQ")
  datInit = datInit %>% rbind(dataHypotheticalSInit)

ggplot(datInit) +
    geom_line(aes(x = week, color = type, y = biomass)) +
    geom_point(aes(x = week, color = type, y = biomass)) +
    ylab("initial biomass after starvation") +
    xlab('weeks in starvation')
    theme_bw() 
  
  
    
ggplot(data1) +
  geom_line(aes(x = week, color = type, y = biomass)) +
  geom_point(aes(x = week, color = type, y = biomass)) +
  ylab("initial biomass") +
  theme_bw() +
  facet_grid(.~time_hours)


data2=data %>%
  group_by(time_hours, week) %>%
  mutate(QtoS = biomass[type == "Q"]/biomass[type=="S"],
         NtoS = biomass[type == "NQ"]/biomass[type=="S"]
  ) %>%
  tidyr::gather("key", "value", QtoS, NtoS)

ggplot(data2, aes(x = time_hours, y =value, col = key )) +
  geom_point() +
  geom_line() +
  facet_grid(week ~.) +
  geom_abline(aes(slope = 0, intercept = 1), col = "black")



initial_values = data %>% filter(week == 1 & time_hours == 0.5) %>%
  select(initial_biomass = biomass, type)
data3=data %>%
  left_join(initial_values) %>%
  group_by(time_hours, week) %>%
  mutate(QtoS = (biomass[type == "Q"]/initial_biomass[type == "Q"])/(biomass[type == "S"]/initial_biomass[type=="S"]),
         NtoS = (biomass[type == "NQ"]/initial_biomass[type == "NQ"])/(biomass[type == "S"]/initial_biomass[type=="S"])) %>%
  tidyr::gather("key", "value", QtoS, NtoS)

ggplot(data3, aes(x = time_hours, y =value, col = key )) +
  geom_point() +
  geom_line() +
  facet_grid(week ~.) +
  geom_abline(aes(slope = 0, intercept = 1), col = "black")
