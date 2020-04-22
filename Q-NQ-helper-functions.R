# # our model of growth (using differential euqutions)
simulate_regrowth_single_strain = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G) 
    Gdot   = -Jg*N 
    Ndot   = a*Jg*N
    return(list(c(Gdot, Ndot)))
  })
}


# FITTING PARAMETERS TO FRESH CELLS
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
# data is taken from the global environment
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
#Parameters data, a,N0 and H0 are taken from the global environment
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
# data, a, Vh, Kh, tlag, N0, H0 are taken from the global env
sumLeastSquaresFitNprop = function(param) {
  Nprop = param[1]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
  difference = data$biomass - simulatedN
  num_obs = length(difference)
  return(sum(difference^2)/num_obs)
}


# Same as above but a, Vh, Kh, tlag,N0 and H0 are taken from the global environment
# This is because function DEoptim can only handle a function of one input = param
# Parameters a,Vh,Kh,N0,H0, data are taken from the global env
sumLeastSquaresFitNprop3 = function(param) {
  Nprop = param[1]
  tlag = param[2]
  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
  difference = data$biomass - simulatedN
  num_obs = length(difference)
  return(sum(difference^2)/num_obs)
}



# This is a function used for MLE fitting
#nll = function(a, Vh, Kh, mu, sigma, tlag, N0, Nprop, H0) {
#  simulatedN = simulate_regrowth_with_lag(a, Vh, Kh,tlag, N0, Nprop, H0)
#  difference = data$biomass - simulatedN
#  logf = suppressWarnings(dnorm(difference[-1],mu,sigma,log = TRUE))
#  ll = sum(logf)
#  return(-ll)
#}