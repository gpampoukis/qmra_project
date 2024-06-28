
## QMRA modelling ----------------------------------------------------------

library(tidyverse)
library(biorisk)

## Initial concentration: https://doi.org/10.1016/j.ijfoodmicro.2014.12.005 ; conventional farms

logN0 <- Normal$new("logN0")$
  map_input("mu", Constant$new("mu_logN0", 0.13))$
  map_input("sigma", Constant$new("sigma_logN0", 0.55))

prevalence <- 7/100

## Inactivation

reductions <- Normal$new("washing")$
  map_input("mu", Constant$new("mu_washing", -3.45))$
  map_input("sigma", Constant$new("sigma_washing", 0.56))

logN_stor <- ElementPlus$new("logN_stor")$
  map_input("a", logN0)$
  map_input("b", reductions)

## Growth model: https://doi.org/10.1016/S0168-1605(02)00252-0

stor_temp <- Triangular$new("refrig_temp")$
  map_input("a", Constant$new("min_temp", 3.9))$
  map_input("b", Constant$new("max_temp", 7.4))$
  map_input("c", Constant$new("most_likely_temp", 5.4))


stor_pH <- Uniform$new("pH")$
  map_input("min", Constant$new("pH_min", 4))$
  map_input("max", Constant$new("pH_max", 5))

mu_stor <- FullRatkowsky_pH_model$new("mu_stor")$
  map_input("temperature", stor_temp)$
  map_input("pH", stor_pH)$
  map_input("b",
            Normal$new("b")$
              map_input("mu", Constant$new("mu_b", 0.2345/log(10)))$
              map_input("sigma", Constant$new("sigma_b", 0.0083/log(10)))
  )$
  map_input("Tmin",
            Normal$new("Tmin")$
              map_input("mu", Constant$new("mu_Tmin", 4.14))$
              map_input("sigma", Constant$new("sigma_Tmin", 0.63))
  )$
  map_input("Tmax",
            Normal$new("Tmax")$
              map_input("mu", Constant$new("mu_Tmax", 49.55))$
              map_input("sigma", Constant$new("sigma_Tmax", 0.42))
  )$
  map_input("c",
            Normal$new("c")$
              map_input("mu", Constant$new("mu_c", 0.2636))$
              map_input("sigma", Constant$new("sigma_c", 0.038))
  )$
  map_input("pHmin",
            Normal$new("pHmin")$
              map_input("mu", Constant$new("mu_pHmin", 3.909))$
              map_input("sigma", Constant$new("sigma_pHmin", 0.031))
  )$
  map_input("pHmax",
            Normal$new("pHmax")$
              map_input("mu", Constant$new("mu_pHmax", 8.86))$
              map_input("sigma", Constant$new("sigma_pHmax", 0.19))
  )

 mu_stor$simulate(100)    
 mu_stor$density_plot()
 plot_model(mu_stor)

## Growth

stor_time <- Triangular$new("stor_time")$
  map_input("a", Constant$new("min_time", 24))$
  map_input("c", Constant$new("max_time", 24*4))$
  map_input("b", Constant$new("most_likely_time", 24*10))

exposure <- ExponentialGrowthNmax$new("exposure")$
  map_input("logN0", logN_stor)$
  map_input("t", stor_time)$
  map_input("mu", mu_stor)$
  map_input("logNmax", Constant$new("logNmax", 8))

 plot_model(exposure)
 
 exposure$simulate(1000)
 exposure$simulations
 exposure$density_plot()

## Dose-response: https://doi.org/10.1111/j.0272-4332.2004.00441.x

dose <- Concentration2Dose$new("dose")$
  map_input("logN", exposure)$
  map_input("size", 
            Uniform$new("serving_size")$
              map_input("min", Constant$new("min_serving", 10))$
              map_input("max", Constant$new("max_serving", 60))
  )

Pill_adults <- DoseResponse_BetaPoisson$new("Pill_adults")$
  map_input("dose", dose)$
  map_input("alpha",
            Constant$new("alpha_adults", 0.0496)
  )$
  map_input("beta",
            Constant$new("beta_adults", 1.001)
  )

Pill_kids <- DoseResponse_BetaPoisson$new("Pill_kids")$
  map_input("dose", dose)$
  map_input("alpha",
            Constant$new("alpha_kids", 0.0844)
  )$
  map_input("beta",
            Constant$new("beta_kids", 1.442)
  )

 Pill_adults$simulate(100)
 Pill_adults$simulations
 
 Pill_kids$simulate(100)
 Pill_kids$simulations
 
 Pill_adults$histogram() + scale_x_log10()
 Pill_kids$histogram() + scale_x_log10()

## Number of cases

cases_adults <- Pill2Cases_N$new("cases")$
  map_input("Pill", Pill_adults)$
  map_input("servings", Constant$new("N", 1e6))

cases_kids <- Pill2Cases_N$new("cases")$
  map_input("Pill", Pill_kids)$
  map_input("servings", Constant$new("N", 1e6))

## Simulations for adults

cases_adults$simulate(1e7, seed = 1242)

quantile_table(exposure, chosen = "exposure", probs = c(.1, .5, .90))

quantile_table(cases_adults, chosen = "cases", probs = c(.1, .5, .90)) %>%
  select(-node) %>%
  pivot_longer(everything()) %>%
  mutate(out = value/1e6*prevalence)

quantile_table(Pill_adults, chosen = "Pill_adults", probs = c(.1, .5, .90)) %>%
  select(-node) %>%
  pivot_longer(everything()) %>%
  mutate(out = value*prevalence)

## Simulations for kids

cases_kids$simulate(1e7, seed = 12142)

# cases_kids$density_plot()

quantile_table(cases_kids, chosen = "cases", probs = c(.1, .5, .9)) %>%
  select(-node) %>%
  pivot_longer(everything()) %>%
  mutate(out = value/1e6*prevalence)

quantile_table(Pill_kids, chosen = "Pill_kids", probs = c(.1, .5, .9)) %>%
  select(-node) %>%
  pivot_longer(everything()) %>%
  mutate(out = value*prevalence)

##

chlorine <- exposure$simulations

