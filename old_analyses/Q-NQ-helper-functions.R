# # our model of growth (using differential euqutions)
simulate_regrowth_single_strain = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G) 
    Gdot   = -Jg*N 
    Ndot   = a*Jg*N
    return(list(c(Gdot, Ndot)))
  })
}

simulate_regrowth_single_strain_lag = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg = Vh*G/(Kh+G) 
    if (is.null(Time)) {
      JgN = 0
    } else{
      if(Time > tlag) {JgN =  Jg*N} else{JgN =  0}
    }
    Gdot   = -JgN 
    Ndot   = a*JgN
    return(list(c(Gdot, Ndot)))
  })
}

simulate_regrowth_two_strains_lag = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg = Vh*G/(Kh+G) 
    if (is.null(Time)) {
      JgN1 = 0
      JgN2 = 0
      } else{
      if(Time > tlag1) {JgN1 =  Jg*N1} else{JgN1 =  0}
      if(Time > tlag2) {JgN2 =  Jg*N2} else{JgN2 =  0}
      }
    Gdot  = -JgN1 - JgN2
    N1dot = a*JgN1
    N2dot = a*JgN2
    return(list(c(Gdot, N1dot, N2dot)))
  })
}

# FITTING PARAMETERS TO FRESH CELLS
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
# data is taken from the global environment
simulate_regrowth_with_lag = function(a, Vh, Kh, tlag, N0, Nprop, H0, data = data) {
  timesGrowth = data %>% rowwise() %>% mutate(time = max(0, time_hours - tlag)) %>% pull(time)
  aliveN0 = Nprop*N0
  deadN0 = (1-Nprop)*N0
  pars <- c(a = a, Vh = Vh, Kh = Kh)
  inits = c(G=H0, N=aliveN0)
  simulation <- as.data.frame(ode(inits, timesGrowth, simulate_regrowth_single_strain, pars)) %>% 
    mutate(N = N + deadN0)
  simulatedN <- simulation$N
  return(simulatedN)
}

simulate_regrowth_with_lag2 = function(a, Vh, Kh, tlag, N0, Nprop, H0, data = data) {
  aliveN0 = Nprop*N0
  deadN0 = (1-Nprop)*N0
  pars <- c(a = a, Vh = Vh, Kh = Kh, tlag = tlag)
  inits = c(G=H0, N=aliveN0)
  simulation <- as.data.frame(ode(inits, data$time_hours, simulate_regrowth_single_strain_lag, pars)) %>% 
    mutate(N = N + deadN0)
  simulatedN <- simulation$N
  return(simulatedN)
}

simulate_regrowth_two_srains_with_lag = function(a, Vh, Kh, tlag1,tlag2, N0, Nprop1, Nprop2, prop1, H0, timesGrowth) {
  N1_init_biomass = N0*prop1
  N2_init_biomass = N0*(1-prop1)
  
  aliveN01 = Nprop1*N1_init_biomass
  deadN01 = (1-Nprop1)*N1_init_biomass
  
  aliveN02 = Nprop2*N2_init_biomass
  deadN02 = (1-Nprop2)*N2_init_biomass
  
  pars <- c(a = a, Vh = Vh, Kh = Kh, tlag1=tlag1, tlag2=tlag2)
  inits = c(G=H0, N1=aliveN01, N2=aliveN02)
  simulation <- as.data.frame(ode(inits, timesGrowth, simulate_regrowth_two_strains_lag, pars)) %>% 
    mutate(N1 = N1 + deadN01,
           N2 = N2 + deadN02,
           N = N1 + N2)
  simulatedN <- simulation$N
  return(simulatedN)
}



# we will use that function to fit the metabolic parameters Vh and Kh which are assumed to be the same for fresh and starved cells. Here we assume tlag is known
# it returns sum of squared errors so we want to find Kh and Vh such that they minimise this function
# Remember: data comes from the global environment. Set this variable before calling this function
# Here data should be a data frame of times and biomass: one observation per each time point
sumLeastSquaresFitGrowth = function(param, a, tlag, N0, H0) {
  Vh = param[1]
  Kh = param[2]
  simulatedN = simulate_regrowth_with_lag2(a, Vh, Kh,tlag, N0, Nprop=1, H0, data)
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
  simulatedN = simulate_regrowth_with_lag2(a, Vh, Kh,tlag, N0, Nprop=1, H0, data)
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
  simulatedN = simulate_regrowth_with_lag2(a, Vh, Kh,tlag, N0, Nprop, H0, data)
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
  simulatedN = simulate_regrowth_with_lag2(a, Vh, Kh,tlag, N0, Nprop, H0, data)
  difference = data$biomass - simulatedN
  num_obs = length(difference)
  return(sum(difference^2)/num_obs)
}


sumLeastSquaresFitNpropFirst10Obs = function(param) {
  Nprop = param[1]
  tlag = param[2]
  simulatedN = simulate_regrowth_with_lag2(a, Vh, Kh,tlag, N0, Nprop, H0, data)
  difference = head(data$biomass - simulatedN,10)
  return(sum(difference^2)/10)
}


GetMeanValues = function(data, cols_to_group_by) {
  averagedData = data %>%
    group_by_at(cols_to_group_by) %>%
    summarise(meanBiomass = mean(biomass, na.rm = TRUE))
  return(averagedData)
}



sumLeastSquaresFitPropQAndNQ = function(param) {
  Qprop = param[1]
  simulated = simulate_regrowth_two_srains_with_lag(a, Vh, Kh, tlagQ,tlagNQ, N0, NpropQ, NpropNQ, Qprop, H0, data$time_hours)
  difference = data$biomass - simulated
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