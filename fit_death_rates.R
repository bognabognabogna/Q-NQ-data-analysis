library(dplyr)
library(ggplot2)

# from Allen
# N is CFU/CFU(0) we expect that CFU/CFU(0) = exp(-dt)
data_Q = data.frame(days  = c(7,14,21,28), N=c(100, 65, 65, 15)/100) %>% mutate(type = "Q", time = 24*days, logN=log(N))
data_NQ = data.frame(days  = c(7,14,21,28), N=c(35, 15, 10, 2)/100) %>% mutate(type = "NQ", time = 24*days, logN=log(N))
data = rbind(data_Q,data_NQ)

mod = lm(logN ~ 0+ time:type, data = data)
data$predicted = predict(mod)


ggplot(data = data %>%  
         rename(experimental.biomass = N) %>%
         mutate(predicted.biomass = exp(predicted)) %>% 
         select(days, experimental.biomass, predicted.biomass, type) %>%
         tidyr::gather(key = "data.type", value = "Biomass", predicted.biomass, experimental.biomass), 
       aes(x = days, y = Biomass, col = data.type)) +
  geom_point() +
  geom_line() +
  ylab('% CFU') +
  xlab('Time [days]') +
  facet_wrap("type")

summary(mod)

# so we have dQ = 0.002
# dNQ = 0.005


## do the same for our data
SURVIVING_DATA_PATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/surviving_data.txt"
surviving_data = read.table(SURVIVING_DATA_PATH, header = TRUE)
TIME_REGROWTH = 3
CULTURE = "H2O"
ALLEN_WEEK_1_PROP_Q = 1
ALLEN_WEEK_1_PROP_NQ = 0.35


our_surviving_data = surviving_data %>% 
  filter(time_hours == TIME_REGROWTH & starvation_medium == CULTURE) %>%
 #mutate(surviving = ifelse(cell_type == "Q", surviving*ALLEN_WEEK_1_PROP_Q, surviving*ALLEN_WEEK_1_PROP_NQ)) %>%
  mutate(N = surviving/100, logN=log(N), days = week*7, time = days*24) %>%
  select(logN, time, days, type = cell_type, N)

mod_our_data = lm(logN ~ 0+ time:type, data = our_surviving_data)
our_surviving_data$predicted = predict(mod_our_data)


ggplot(data = our_surviving_data %>%  
         rename(experimental.biomass = N) %>%
         mutate(predicted.biomass = exp(predicted)) %>% 
         select(days, experimental.biomass, predicted.biomass, type) %>%
         tidyr::gather(key = "data.type", value = "Biomass", predicted.biomass, experimental.biomass), 
       aes(x = days, y = Biomass, col = data.type)) +
  geom_point() +
  geom_line() +
  ylab('% CFU') +
  xlab('Time [days]') +
  facet_wrap("type")

summary(mod_our_data)
# this is not a problem that we start at week 1: the estimated death rates should be still ok
# we need to justfiy that we fit it to regrowth after 3 hours!

# dQ = 0.001;
# dNQ = 0.0007;



# Now do the same for Lee Leu data
LEELEU_DATA_PATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Q-NQ/data/LeeLeu - rebudding freq.txt"
leeleu_data = read.table(LEELEU_DATA_PATH, header = TRUE)


our_leeleu_data = leeleu_data %>% 
  mutate(N = rebudding/100, logN=log(N),  time = days*24) %>%
  select(logN, time, days, type = QNQ, N)

mod_leeleu_data = lm(logN ~ 0+ time:type, data = our_leeleu_data)
our_leeleu_data$predicted = predict(mod_leeleu_data)


ggplot(data = our_leeleu_data %>%  
         rename(experimental.biomass = N) %>%
         mutate(predicted.biomass = exp(predicted)) %>% 
         select(days, experimental.biomass, predicted.biomass, type) %>%
         tidyr::gather(key = "data.type", value = "Biomass", predicted.biomass, experimental.biomass), 
       aes(x = days, y = Biomass, col = data.type)) +
  geom_point() +
  geom_line() +
  ylab('% CFU') +
  xlab('Time [days]') +
  facet_wrap("type")

summary(mod_leeleu_data)

# dQ = 0.0001
# dNQ = 0.003

