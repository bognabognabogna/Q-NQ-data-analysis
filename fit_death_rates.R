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