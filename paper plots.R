## Plots for paper ####

setwd("C:/Users/dattas/OneDrive - NIWA/Projects/2025/OCES2501 Size-based modelling/FishMIP/impacts-of-future-climate/")
rm(list = ls())

# Libraries

library(therMizer) # remotes::install_github("sizespectrum/therMizer") # if needed to install
library(tidyverse)
library(readxl)
library(gridExtra)

# Load in functions

source('functions for paper.R')

## Do run if needed ####

do_run = T

if (do_run == T) {

  ## Read in models and data

  # Mizer models

  # TBGB

  #Read in model
  tbgb_model = readRDS('inputs/tbgb_model_modified.rds')
  species_names_tbgb <- as.character(tbgb_model@species_params$species) # getting species names
  sizes_tbgb <- tbgb_model@w # getting w vector
  full_x_tbgb <- log10(tbgb_model@w_full) # 211 long

  temp = tbgb_model@species_params %>%
    as_tibble() %>%
    select(species, 'Group name', 'Main species') %>%
    arrange(species)

  # CR

  #Read in model
  cr_model <- readParams("inputs/cr_params_final_fit_updated.RDS") # pulling in parameters
  species_names_cr <- as.character(cr_model@species_params$species) # getting species names
  sizes_cr <- cr_model@w # getting w vector
  full_x_cr <- log10(cr_model@w_full) # 211 long

  temp = cr_model@species_params %>%
    as_tibble() %>%
    select(species, 'Group name', 'Main species') %>%
    arrange(species)

  # Thermal tolerances from Fishbase

  all_thermal_tolerances = readRDS('inputs/thermal_tolerances.rds')

  tbgb_tolerances = all_thermal_tolerances %>%
    filter(area == 'Tasman and Golden Bay') %>%
    select(species, temp_min, temp_max)

  tbgb_model@species_params = tbgb_model@species_params %>%
    left_join(tbgb_tolerances, by = 'species')

  cr_tolerances = all_thermal_tolerances %>%
    filter(area == 'Chatham Rise') %>%
    select(species, temp_min, temp_max)

  cr_model@species_params = cr_model@species_params %>%
    left_join(cr_tolerances, by = 'species')

  same_species = tbgb_model@species_params %>%
    inner_join(cr_model@species_params, by = 'species') %>%
    pull(species) # find overlapping species

  # Ocean temperatures by depth

  all_temps_tbgb = readRDS(file = "inputs/tbgb_all_temperatures.rds") %>%
    relocate(BOT, .after = 7)
  all_temps_cr = readRDS(file = "inputs/cr_all_temperatures.rds") %>%
    relocate(BOT, .after = 7)

  # Add the species thermal tolerances, space throughout time period

  to_add = tbgb_model@species_params %>%
    select(species, temp_min, temp_max) %>%
    arrange(temp_min, temp_max) %>%
    mutate(date = round_date(seq(min(all_temps_tbgb$date), max(all_temps_tbgb$date), length.out = nrow(tbgb_model@species_params)), '1 month'))
  all_temps_tbgb_with_tols = all_temps_tbgb %>% left_join(to_add, by = 'date') %>%
    pivot_longer(-c('date', 'date_adjust', 'species', 'temp_min', 'temp_max')) %>%
    mutate(name = factor(name, levels = c('SST', '10m', '20m', '35m', 'BOT')))

  to_add = cr_model@species_params %>%
    select(species, temp_min, temp_max) %>%
    arrange(temp_min, temp_max) %>%
    mutate(date = round_date(seq(min(all_temps_cr$date), max(all_temps_cr$date), length.out = nrow(cr_model@species_params)), '1 month'))
  all_temps_cr_with_tols = all_temps_cr %>% left_join(to_add, by = 'date') %>%
    pivot_longer(-c('date', 'date_adjust', 'species', 'temp_min', 'temp_max')) %>%
    mutate(name = factor(name, levels = c('SST', '100m', '400m', '600m', 'BOT')))

  # Feed temperatures into the ocean temp array

  # TBGB

  time_steps = as.character(all_temps_tbgb$date_adjust) # vector of time steps
  realm_names_tbgb = c('SST', '10m', '20m', '35m', 'BOT') # use realm names from all_temps

  new_ocean_temp_array <- array(NA, dim = c(nrow(all_temps_tbgb),
                                            length(realm_names_tbgb)),
                                dimnames = list(time = time_steps, realm = realm_names_tbgb)) # need date as character to work

  new_ocean_temp_array[, 1] <- all_temps_tbgb$SST
  new_ocean_temp_array[, 2] <- all_temps_tbgb$`10m`
  new_ocean_temp_array[, 3] <- all_temps_tbgb$`20m`
  new_ocean_temp_array[, 4] <- all_temps_tbgb$`35m`
  new_ocean_temp_array[, 5] <- all_temps_tbgb$BOT

  new_ocean_temp_array_tbgb = new_ocean_temp_array # copy to TBGB

  # CR

  realm_names_cr = c('SST', '100m', '400m', '600m', 'BOT')
  # use realm names from all_temps

  new_ocean_temp_array <- array(NA, dim = c(nrow(all_temps_cr),
                                            length(realm_names_cr)),
                                dimnames = list(time = time_steps, realm = realm_names_cr)) # need date as character to work

  new_ocean_temp_array[, 1] <- all_temps_cr$SST
  new_ocean_temp_array[, 2] <- all_temps_cr$`100m`
  new_ocean_temp_array[, 3] <- all_temps_cr$`400m`
  new_ocean_temp_array[, 4] <- all_temps_cr$`600m`
  new_ocean_temp_array[, 5] <- all_temps_cr$BOT

  new_ocean_temp_array_cr = new_ocean_temp_array

  # set up the vertical_migration array

  # TBGB

  # Create the vertical migration array and fill it
  vertical_migration_array <- array(0, dim = (c(length(realm_names_tbgb),
                                                length(species_names_tbgb),
                                                length(sizes_tbgb))),
                                    dimnames = list(realm = realm_names_tbgb, sp = species_names_tbgb, w = signif(sizes_tbgb, 4))) #realm x species x size

  # The vertical migration array we have created is currently filled with 0s

  # Use data from Vidette

  for (j in 1:nrow(tbgb_model@species_params)) {
    vert_mixing = read_excel('inputs/tbgb_Initial_biomass_distributions.xlsx',
                             sheet = tbgb_model@species_params$species[j],
                             range = "A5:H30") %>%
      select(contains('Layer')) %>%
      colMeans()
    vertical_migration_array[, j, ] = vert_mixing/sum(vert_mixing)
  }

  vertical_migration_array_tbgb = vertical_migration_array

  # CR

  # Create the vertical migration array and fill it
  vertical_migration_array <- array(0, dim = (c(length(realm_names_cr),
                                                length(species_names_cr),
                                                length(sizes_cr))),
                                    dimnames = list(realm = realm_names_cr, sp = species_names_cr, w = signif(sizes_cr, 4))) #realm x species x size

  # vertical_migration_array [1, , ] <- 1 # all in top layer
  # vertical_migration_array [, , ] <- 1/dim(vertical_migration_array)[1] # equally mixed

  # Use data from Vidette
  cr_sheets <- excel_sheets('inputs/cr_Initial_biomass_distributions.xlsx') # list sheet names

  for (j in 1:nrow(cr_model@species_params)) {
    vert_mixing = read_excel('inputs/cr_Initial_biomass_distributions.xlsx',
                             sheet = which(str_detect(sheets, cr_model@species_params$species[j])), # to find correct sheet
                             range = "A5:H29") %>%
      select(contains('Layer')) %>%
      colMeans()
    vertical_migration_array[, j, ] = vert_mixing/sum(vert_mixing)
  }

  vertical_migration_array_cr = vertical_migration_array

  # exposure array is all 1s

  # TBGB

  exposure_array <- array(0, dim = (c(length(realm_names_tbgb),
                                      length(species_names_tbgb))),
                          dimnames = list (realm = realm_names_tbgb, sp = species_names_tbgb)) #realm x species

  exposure_array[, ] = 1 # present in all layers

  exposure_array_tbgb = exposure_array

  # CR

  exposure_array <- array(0, dim = (c(length(realm_names_cr),
                                      length(species_names_cr))),
                          dimnames = list (realm = realm_names_cr, sp = species_names_cr)) #realm x species

  exposure_array[, ] = 1 # present in all layers

  exposure_array_cr = exposure_array

  # Resource spectrum


  # Use Kieran's method to rescale plankton up to original mizer abundances.

  # Below we change the array that will be used to force the plankton dynamics through time to match the scale of the original mizer model resource.

  # Load saved ISIMIP spectra

  # TBGB
  out_isimip = read.table(file = "inputs/tbgb_resource_spectra_fitted.dat")
  out_isimip = as(out_isimip, "matrix")
  rownames(out_isimip) <- time_steps
  colnames(out_isimip) <- signif(full_x_tbgb, 3)

  n_pp_array_rescaled_tbgb <- array(NA, dim = c(nrow(out_isimip), ncol(out_isimip)), dimnames = list(time = time_steps, w = signif(full_x_tbgb, 3))) # setup an empty array for the corrected plankton

  stable_period_rows <- which(all_temps_tbgb$date >= '1960-01-01' & all_temps_tbgb$date < '1970-01-01') # stable period rows

  for (j in seq(1,length(time_steps), 1)) {
    n_pp_array_rescaled_tbgb[j, ] <- (out_isimip[j, ] - colMeans(out_isimip[stable_period_rows, ])) + log10(tbgb_model@initial_n_pp*tbgb_model@dw_full) #subtracting anomaly and adding back the mean
  }
  n_pp_array_rescaled_tbgb[, w_full(tbgb_model) >= resource_params(tbgb_model)$w_pp_cutoff] <- 0 # cut off above w_pp_cutoff

  # for 1961-1970

  ss_n_pp_by_layer = n_pp_array_rescaled_tbgb[1:120, ] %>% # limit to 1961-1970
    colMeans() # take mean of each weight bin
  ss_n_pp_array_rescaled_tbgb = n_pp_array_rescaled_tbgb # take first two rows
  for (j in 1:nrow(ss_n_pp_array_rescaled_tbgb)) ss_n_pp_array_rescaled_tbgb[j, ] = ss_n_pp_by_layer # assign values


  # CR
  out_isimip = read.table(file = "inputs/cr_resource_spectra_fitted.dat")
  out_isimip = as(out_isimip, "matrix")
  rownames(out_isimip) <- time_steps
  colnames(out_isimip) <- signif(full_x_cr, 3)

  n_pp_array_rescaled_cr <- array(NA, dim = c(nrow(out_isimip), ncol(out_isimip)), dimnames = list(time = time_steps, w = signif(full_x_cr, 3))) # setup an empty array for the corrected plankton

  for (j in seq(1,length(time_steps), 1)) {
    n_pp_array_rescaled_cr[j, ] <- (out_isimip[j, ] - colMeans(out_isimip[stable_period_rows, ])) + log10(cr_model@initial_n_pp*cr_model@dw_full) #subtracting anomaly and adding back the mean
  }
  n_pp_array_rescaled_cr[, w_full(cr_model) >= resource_params(cr_model)$w_pp_cutoff] <- 0 # cut off above w_pp_cutoff

  # for 1961-1970

  ss_n_pp_by_layer = n_pp_array_rescaled_cr[1:120, ] %>% # limit to 1961-1970
    colMeans() # take mean of each weight bin
  ss_n_pp_array_rescaled_cr = n_pp_array_rescaled_cr # take first two rows
  for (j in 1:nrow(ss_n_pp_array_rescaled_cr)) ss_n_pp_array_rescaled_cr[j, ] = ss_n_pp_by_layer # assign values

  # Adjusting tolerances to match depth presence

  ## TBGB

  # temperature
  ss_temp_by_layer = new_ocean_temp_array_tbgb[1:120, ] %>% # limit to 1961-1970
    colMeans() # take mean of each layer
  ss_ocean_temp_array_tbgb = new_ocean_temp_array_tbgb # take first two rows (can't just use one for some reason)
  for (j in 1:nrow(ss_ocean_temp_array_tbgb)) ss_ocean_temp_array_tbgb[j, ] = ss_temp_by_layer # assign values

  ss_temp_tibble_tbgb = t(ss_ocean_temp_array_tbgb[1, ]) %>%
    as_tibble() %>%
    add_column(time = 0) %>%
    pivot_longer(-time, values_to = 'mean') %>%
    mutate(name = factor(name, levels = c('SST', '10m', '20m', '35m', 'BOT'))) %>%
    select(-time)

  # Adjust thermal tolerances so they always overlap with all five layers

  tolerances_adjusted_tbgb = tibble(
    species = tbgb_model@species_params$species,
    temp_min = tbgb_model@species_params$temp_min,
    temp_max = tbgb_model@species_params$temp_max) %>%
    mutate(temp_min_adjust = case_when(temp_min > min(ss_temp_tibble_tbgb$mean) ~ min(ss_temp_tibble_tbgb$mean), T ~ temp_min),
           temp_max_adjust = case_when(temp_max < max(ss_temp_tibble_tbgb$mean) ~ max(ss_temp_tibble_tbgb$mean), T ~ temp_max)) %>%
    relocate(temp_min_adjust, .after = temp_min) %>%
    relocate(temp_max_adjust, .after = temp_max) %>%
    arrange(temp_min_adjust, temp_max_adjust) %>%
    mutate(date = round_date(seq(min(all_temps_tbgb$date), max(all_temps_tbgb$date), length.out = nrow(tbgb_model@species_params)), '1 month'))


  ## Add species tolerances, space throughout time period

  all_temps_tbgb_with_tols = all_temps_tbgb %>% left_join(tolerances_adjusted_tbgb, by = 'date') %>%
    pivot_longer(-c('date', 'date_adjust', 'species', 'temp_min', 'temp_max', 'temp_min_adjust', 'temp_max_adjust')) %>%
    mutate(name = factor(name, levels = c('SST', '10m', '20m', '35m', 'BOT'))) %>%
    left_join(ss_temp_tibble_tbgb, by = 'name')

  ## CR

  # temperature
  ss_temp_by_layer = new_ocean_temp_array_cr[1:120, ] %>% # limit to 1961-1970
    colMeans() # take mean of each layer
  ss_ocean_temp_array_cr = new_ocean_temp_array_cr # take first two rows (can't just use one for some reason)
  for (j in 1:nrow(ss_ocean_temp_array_cr)) ss_ocean_temp_array_cr[j, ] = ss_temp_by_layer # assign values

  ss_temp_tibble_cr = t(ss_ocean_temp_array_cr[1, ]) %>%
    as_tibble() %>%
    add_column(time = 0) %>%
    pivot_longer(-time, values_to = 'mean') %>%
    mutate(name = factor(name, levels = c('SST', '100m', '400m', '600m', 'BOT'))) %>%
    select(-time)

  # Adjust thermal tolerances so they always overlap with all five layers
  tolerances_adjusted_cr = tibble(
    species = cr_model@species_params$species,
    temp_min = cr_model@species_params$temp_min,
    temp_max = cr_model@species_params$temp_max) %>%
    mutate(time_at_surface = vertical_migration_array_cr[1, , 1],
           time_at_bottom = vertical_migration_array_cr[5, , 1]) # add time at surface and bottom

  # Go through species by species, if it spends >=5% of time at surface then include at temp_max,
  # if it spends >=5% of time at bottom then include at temp_min

  tolerances_adjusted_cr = tolerances_adjusted_cr %>%
    mutate(temp_min_adjust = case_when(species == 'PFS' ~ 10, # manually set
                                       time_at_bottom >= 0.05 & temp_min > min(ss_temp_tibble_cr$mean) ~ min(ss_temp_tibble_cr$mean), T ~ temp_min),
           temp_max_adjust = case_when(species == 'PFL' ~ 16, # manually set
                                       time_at_surface >= 0.05 & temp_max < max(ss_temp_tibble_cr$mean)  ~ max(ss_temp_tibble_cr$mean), T ~ temp_max)) %>%
    relocate(temp_min_adjust, .after = temp_min) %>%
    relocate(temp_max_adjust, .after = temp_max) %>%
    arrange(temp_min_adjust, temp_max_adjust) %>%
    mutate(date = round_date(seq(min(all_temps_cr$date), max(all_temps_cr$date), length.out = nrow(cr_model@species_params)), '1 month'))

  all_temps_cr_with_tols = all_temps_cr %>% left_join(tolerances_adjusted_cr, by = 'date') %>%
    pivot_longer(-c('date', 'date_adjust', 'species', 'temp_min', 'temp_max', 'temp_min_adjust', 'temp_max_adjust', 'time_at_surface', 'time_at_bottom')) %>%
    mutate(name = factor(name, levels = c('SST', '100m', '400m', '600m', 'BOT'))) %>%
    left_join(ss_temp_tibble_cr, by = 'name')

  # Add adjusted thermal tolerances to species params objects

  tbgb_model@species_params = tbgb_model@species_params %>%
    left_join(tolerances_adjusted_tbgb %>% select(species, temp_min_adjust, temp_max_adjust), by = 'species')

  cr_model@species_params = cr_model@species_params %>%
    left_join(tolerances_adjusted_cr %>% select(species, temp_min_adjust, temp_max_adjust), by = 'species')

  # Plot thermal tolerances

  plot_temp_plus_tolerance(all_temps_tbgb_with_tols, same_species, title_pick = 'TBGB mean temperature', plot_mean = F, include_tolerance = F) # test without tolerances

  # before adjustment

  g1 = plot_temp_plus_tolerance(all_temps_tbgb_with_tols, same_species, title_pick = 'TBGB temperatures and thermal tolerances')

  g2 = plot_temp_plus_tolerance(all_temps_cr_with_tols, same_species, title = 'CR temperatures and thermal tolerances')

  both = grid.arrange(g1, g2, nrow = 2)
  ggsave(filename = "paper/PNGs/Tolerances unadjusted.png", both, width = 9, height = 9)

  # after adjustment

  g3 = plot_temp_plus_tolerance(all_temps_tbgb_with_tols, same_species, adjust_tol = T,
                                title_pick = 'TBGB temperatures and adjusted thermal tolerances',
                                dodge_labels = T)


  g4 = plot_temp_plus_tolerance(all_temps_cr_with_tols, same_species, adjust_tol = T,
                                title = 'CR temperatures and adjusted thermal tolerances',
                                dodge_labels = T)

  both = grid.arrange(g3, g4, nrow = 2)
  ggsave(filename = "paper/PNGs/Tolerances adjusted.png", both, width = 9, height = 9)

  # before and after with mean temperatures

  g1 = plot_temp_plus_tolerance(all_temps_tbgb_with_tols, same_species, title_pick = 'TBGB mean temperatures and thermal tolerances', plot_mean = T)
  g2 = plot_temp_plus_tolerance(all_temps_cr_with_tols, same_species, title = 'CR mean temperatures and thermal tolerances', plot_mean = T)
  g3 = plot_temp_plus_tolerance(all_temps_tbgb_with_tols, same_species, adjust_tol = T,
                                title_pick = 'TBGB adjusted thermal tolerances', plot_mean = T)
  g4 = plot_temp_plus_tolerance(all_temps_cr_with_tols, same_species, adjust_tol = T,
                                title = 'CR adjusted thermal tolerances', plot_mean = T)

  both = grid.arrange(g1, g2, g3, g4, nrow = 2)
  ggsave(filename = "paper/PNGs/Both adjusting tolerances.png", both, width = 16, height = 12)


  ## Running simulations

  # Mizer steady state

  # TBGB

  sim_orig = project(tbgb_model) # run simulation

  temp_model = tbgb_model
  initialN(temp_model) = sim_orig@n[101, , ]
  initialNResource(temp_model) = sim_orig@n_pp[101, ]

  sim_repeat_tbgb = project(temp_model) # run simulation

  # CR

  sim_orig = project(cr_model) # run simulation

  temp_model = cr_model
  initialN(temp_model) = sim_orig@n[101, , ]
  initialNResource(temp_model) = sim_orig@n_pp[101, ]

  sim_repeat_cr = project(temp_model) # run simulation

  # Plotting

  g1 = plotSpectra(sim_repeat_tbgb, total = T) +
    labs(title = 'Tasman and Golden Bay spectra') +
    theme_classic()

  g2 = plotBiomassObservedVsModel(sim_repeat_tbgb, ratio = T) +
    labs(title = 'Tasman and Golden Bay biomass ratios') + # basically the same
    theme_classic()

  g3 = plotSpectra(sim_repeat_cr, total = T) +
    labs(title = 'Chatham Rise spectra') +
    theme_classic()

  g4 = plotBiomassObservedVsModel(sim_repeat_cr, ratio = T) +
    labs(title = 'Chatham Rise biomass ratios') + # basically the same
    theme_classic()

  both = grid.arrange(g1, g2, g3, g4, nrow = 2, layout_matrix= rbind(c(1, 3), c(2, 4)))
  ggsave(filename = "paper/PNGs/Both steady states.png", both, width = 14, height = 10)


  ## TherMizer steady state using average of 1961-1970 temperature and resource spectrum

  # Tasman Bay and Golden Bay

  ss_tbgb_therm = upgradeTherParams(tbgb_model,
                                    temp_min = tbgb_model@species_params$temp_min_adjust,
                                    temp_max = tbgb_model@species_params$temp_max_adjust,
                                    ocean_temp_array = ss_ocean_temp_array_tbgb,
                                    n_pp_array = ss_n_pp_array_rescaled_tbgb,
                                    vertical_migration_array = vertical_migration_array_tbgb,
                                    exposure_array = exposure_array_tbgb,
                                    aerobic_effect = TRUE,
                                    metabolism_effect = TRUE)

  initialNResource(ss_tbgb_therm) = 10^(other_params(ss_tbgb_therm)$n_pp_array[1, ])/ss_tbgb_therm@dw_full # set initial resource (necessary?)

  # Check steady state is possible
  # sim_test = projectToSteady(ss_tbgb_therm, t_max = 50, 1e-8, return_sim = T)
  # sim_test = project(ss_tbgb_therm, t_max = 50) # if not, see what is dying off
  # plotBiomass(sim_test)

  # Now try running to steady state, simulate and take final population - rinse and repeat
  ss_tbgb_therm = ss_tbgb_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_test = project(ss_tbgb_therm, t_max = 100)
  ss_tbgb_therm@initial_n = sim_test@n[dim(sim_test@n)[1],,] # use final population as starting point

  ss_tbgb_therm = ss_tbgb_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_test = project(ss_tbgb_therm, t_max = 100)
  ss_tbgb_therm@initial_n = sim_test@n[dim(sim_test@n)[1],,] # use final population as starting point

  ss_tbgb_therm = ss_tbgb_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_after_tune = project(ss_tbgb_therm, t_max = 200)

  # Chatham Rise


  ss_cr_therm = upgradeTherParams(cr_model,
                                  temp_min = cr_model@species_params$temp_min_adjust,
                                  temp_max = cr_model@species_params$temp_max_adjust,
                                  ocean_temp_array = ss_ocean_temp_array_cr,
                                  n_pp_array = ss_n_pp_array_rescaled_cr,
                                  vertical_migration_array = vertical_migration_array_cr,
                                  exposure_array = exposure_array_cr,
                                  aerobic_effect = TRUE,
                                  metabolism_effect = TRUE)
  initialNResource(ss_cr_therm) = 10^(other_params(ss_cr_therm)$n_pp_array[1, ])/ss_cr_therm@dw_full # set initial resource (necessary?)

  # Check steady state is possible
  # sim_test = projectToSteady(ss_cr_therm, t_max = 200, 1e-8, return_sim = T)
  # sim_test = project(ss_cr_therm, t_max = 50) # if not, see what is dying off
  # plotBiomass(sim_test)

  # Now try running to steady state, simulate and take final population - rinse and repeat
  ss_cr_therm = ss_cr_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_test = project(ss_cr_therm, t_max = 100)
  ss_cr_therm@initial_n = sim_test@n[dim(sim_test@n)[1],,] # use final population as starting point

  ss_cr_therm = ss_cr_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_test = project(ss_cr_therm, t_max = 100)
  ss_cr_therm@initial_n = sim_test@n[dim(sim_test@n)[1],,] # use final population as starting point

  ss_cr_therm = ss_cr_therm |>calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> steady() # crashes, not enough reproduction

  sim_after_tune2 = project(ss_cr_therm, t_max = 200)

  # Plot steady states

  # TBGB

  g5 = plotSpectra(sim_after_tune) +
    labs(title = 'Tasman Bay and Golden Bay spectra') +
    theme_classic()
  plotBiomass(sim_after_tune, total = T) +
    labs(title = 'Tasman Bay and Golden Bay') +
    theme_classic()
  g6 = plotBiomassObservedVsModel(sim_after_tune, ratio = T) +
    labs(title = 'Tasman Bay and Golden Bay biomass ratios') +
    theme_classic()

  # CR

  g7 = plotSpectra(sim_after_tune2) +
    labs(title = 'Chatham Rise spectra') +
    theme_classic()
  plotBiomass(sim_after_tune2, total = T) +
    labs(title = 'Chatham Rise') +
    theme_classic()
  g8 = plotBiomassObservedVsModel(sim_after_tune2, ratio = T) +
    labs(title = 'Chatham Rise biomass ratios') +
    theme_classic()

  both = grid.arrange(g5, g6, g7, g8, nrow = 2, layout_matrix= rbind(c(1, 3), c(2, 4)))
  ggsave(filename = "paper/Both steady states therMizer.png", both, width = 12, height = 9)


  ## Running therMizer simulations

  tbgb_therm = upgradeTherParams(ss_tbgb_therm,
                                 temp_min = ss_tbgb_therm@species_params$temp_min_adjust,
                                 temp_max = ss_tbgb_therm@species_params$temp_max_adjust,
                                 ocean_temp_array = new_ocean_temp_array_tbgb,
                                 n_pp_array = n_pp_array_rescaled_tbgb,
                                 vertical_migration_array = vertical_migration_array_tbgb,
                                 exposure_array = exposure_array_tbgb,
                                 aerobic_effect = TRUE,
                                 metabolism_effect = TRUE)
  initialN(tbgb_therm) = sim_after_tune@n[201, , ] # set initial distribution
  initialNResource(tbgb_therm) = 10^(other_params(tbgb_therm)$n_pp_array[1, ])/tbgb_therm@dw_full # set initial resource

  cr_therm = upgradeTherParams(ss_cr_therm,
                               temp_min = ss_cr_therm@species_params$temp_min_adjust,
                               temp_max = ss_cr_therm@species_params$temp_max_adjust,
                               ocean_temp_array = new_ocean_temp_array_cr,
                               n_pp_array = n_pp_array_rescaled_cr,
                               vertical_migration_array = vertical_migration_array_cr,
                               exposure_array = exposure_array_cr,
                               aerobic_effect = TRUE,
                               metabolism_effect = TRUE)
  initialN(cr_therm) = sim_after_tune2@n[201, , ] # set initial distribution
  initialNResource(cr_therm) = 10^(other_params(cr_therm)$n_pp_array[1, ])/cr_therm@dw_full # set initial resource (necessary?)

  # Now try running simulation

  # TBGB

  sim_tbgb = project(tbgb_therm, t_max = 600)
  temp = pull_out_info(list(sim_tbgb), real_dates = all_temps_tbgb$date)
  g1 = plot_species_biomass(temp$by_species, title_pick = 'Tasman Bay and Golden Bay')

  # CR

  sim_cr = project(cr_therm, t_max = 600)
  temp2 = pull_out_info(list(sim_cr), real_dates = all_temps_tbgb$date)
  g2 = plot_species_biomass(temp2$by_species, title_pick = 'Chatham Rise')

  both = grid.arrange(g1, g2, nrow = 2)
  ggsave(filename = "paper/Both therMizer historical.png", both, width = 9, height = 9)

  # Let's look at the thermal tolerances of all the species in the model

  plotTherPerformance(tbgb_therm)
  plotTherPerformance(cr_therm)

  # Variation in each species in 2005-2010

  g1 = plot_biomass_variation(sim_tbgb, dates_vec = all_temps_tbgb$date, min_date = '2006-01-01',
                              title_pick = 'TBGB range of biomass 2006 - 2010')
  g2 = plot_biomass_variation(sim_cr, dates_vec = all_temps_cr$date, min_date = '2006-01-01',
                              title_pick = 'CR range of biomass 2006 - 2010')
  both = grid.arrange(g1, g2, nrow = 2)
  ggsave(filename = "paper/Both historical variation.png", both, width = 12, height = 12)


  ## Future climate forcing

  # Stretch out time to 2060

  # climate change vector

  temperature_increase_per_year <- 3/80 # rule of thumb - 3°C from 2021-2100 (i.e., around 80 years)
  temp_increase_vector = temperature_increase_per_year*((1:600)/12) # 50 years, from 2010 to 2060

  # TBGB

  chunk = all_temps_tbgb %>% filter(date >= '2006-01-01') %>% # take last five years
    slice(rep(1:n(), 10)) %>% # duplicate ten times
    mutate(date_adjust = 601:1200,
           date = seq(ymd('2011-01-01'), ymd('2060-12-01'), by = '1 month')) # shift dates forward

  projection_nocc_tbgb = all_temps_tbgb %>%
    add_row(chunk) %>%
    mutate(scenario = 'Baseline')

  dates_vec = projection_nocc_tbgb$date # date vector of entire time range

  chunk_cc = chunk %>%
    mutate(across(SST:BOT, ~ .x + temp_increase_vector))

  projection_cc_tbgb = all_temps_tbgb %>%
    add_row(chunk_cc) %>%
    mutate(scenario = 'With climate change')

  t1 = projection_nocc_tbgb %>%
    add_row(projection_cc_tbgb) %>%
    filter(date >= '2010-12-01') %>%
    pivot_longer(-c(scenario, date, date_adjust)) %>%
    mutate(name = factor(name, levels = c('SST', '10m', '20m', '35m', 'BOT')))

  # CR

  chunk = all_temps_cr %>% filter(date >= '2006-01-01') %>% # take last five years
    slice(rep(1:n(), 10)) %>% # duplicate ten times
    mutate(date_adjust = 601:1200,
           date = seq(ymd('2011-01-01'), ymd('2060-12-01'), by = '1 month')) # shift dates forward

  projection_nocc_cr = all_temps_cr %>%
    add_row(chunk) %>%
    mutate(scenario = 'Baseline')

  chunk_cc = chunk %>%
    mutate(across(SST:BOT, ~ .x + temp_increase_vector))

  projection_cc_cr = all_temps_cr %>%
    add_row(chunk_cc) %>%
    mutate(scenario = 'With climate change')

  t2 = projection_nocc_cr %>%
    add_row(projection_cc_cr) %>%
    filter(date >= '2010-12-01') %>%
    pivot_longer(-c(scenario, date, date_adjust)) %>%
    mutate(name = factor(name, levels = c('SST', '100m', '400m', '600m', 'BOT')))

  # plots

  run_colours_cc <- setNames(RColorBrewer::brewer.pal(4, "Set1")[c(2, 1)],
                             c('Baseline', 'With climate change')) # assign colours to each layer


  g1 = ggplot(t1, aes(x = date, y = value, colour = scenario, linetype = scenario)) +
    geom_line() +
    labs(x = "Year", y = "Temperature (°C)", title = 'TBGB projected temperature', colour = 'Scenario', linetype = 'Scenario') +
    theme_classic() +
    facet_wrap(~name) +
    scale_x_date(breaks = seq(min(t1$date), max(t1$date), by = "10 years"), date_labels = "%Y") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = run_colours_cc)

  g2 = ggplot(t2, aes(x = date, y = value, colour = scenario, linetype = scenario)) +
    geom_line() +
    labs(x = "Year", y = "Temperature (°C)", title = 'CR projected temperature', colour = 'Scenario', linetype = 'Scenario') +
    theme_classic() +
    facet_wrap(~name) +
    scale_x_date(breaks = seq(min(t2$date), max(t2$date), by = "10 years"), date_labels = "%Y") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = run_colours_cc)

  both = grid.arrange(g1, g2, nrow = 2)
  ggsave(filename = "paper/Both climate projections.png", both, width = 9, height = 9)


  # Build the ocean_temp_array and n_pp_array arrays for the simulations

  # Ocean temp arrays

  ## TBGB

  ### no climate change

  ocean_temp_array_tbgb_nocc <- array(NA, dim = c(nrow(projection_nocc_tbgb),
                                                  length(realm_names_tbgb)),
                                      dimnames = list(time = as.character(projection_nocc_tbgb$date_adjust), realm = realm_names_tbgb)) # need date as character to work

  ocean_temp_array_tbgb_nocc[, 1] <- projection_nocc_tbgb$SST
  ocean_temp_array_tbgb_nocc[, 2] <- projection_nocc_tbgb$`10m`
  ocean_temp_array_tbgb_nocc[, 3] <- projection_nocc_tbgb$`20m`
  ocean_temp_array_tbgb_nocc[, 4] <- projection_nocc_tbgb$`35m`
  ocean_temp_array_tbgb_nocc[, 5] <- projection_nocc_tbgb$BOT

  ### with climate change

  ocean_temp_array_tbgb_cc <- array(NA, dim = c(nrow(projection_cc_tbgb),
                                                length(realm_names_tbgb)),
                                    dimnames = list(time = as.character(projection_cc_tbgb$date_adjust), realm = realm_names_tbgb)) # need date as character to work

  ocean_temp_array_tbgb_cc[, 1] <- projection_cc_tbgb$SST
  ocean_temp_array_tbgb_cc[, 2] <- projection_cc_tbgb$`10m`
  ocean_temp_array_tbgb_cc[, 3] <- projection_cc_tbgb$`20m`
  ocean_temp_array_tbgb_cc[, 4] <- projection_cc_tbgb$`35m`
  ocean_temp_array_tbgb_cc[, 5] <- projection_cc_tbgb$BOT

  ## CR

  ### no climate change

  ocean_temp_array_cr_nocc <- array(NA, dim = c(nrow(projection_nocc_cr),
                                                length(realm_names_cr)),
                                    dimnames = list(time = as.character(projection_nocc_cr$date_adjust), realm = realm_names_cr)) # need date as character to work

  ocean_temp_array_cr_nocc[, 1] <- projection_nocc_cr$SST
  ocean_temp_array_cr_nocc[, 2] <- projection_nocc_cr$`100m`
  ocean_temp_array_cr_nocc[, 3] <- projection_nocc_cr$`400m`
  ocean_temp_array_cr_nocc[, 4] <- projection_nocc_cr$`600m`
  ocean_temp_array_cr_nocc[, 5] <- projection_nocc_cr$BOT

  ### with climate change

  ocean_temp_array_cr_cc <- array(NA, dim = c(nrow(projection_cc_cr),
                                              length(realm_names_cr)),
                                  dimnames = list(time = as.character(projection_cc_cr$date_adjust), realm = realm_names_cr)) # need date as character to work

  ocean_temp_array_cr_cc[, 1] <- projection_cc_cr$SST
  ocean_temp_array_cr_cc[, 2] <- projection_cc_cr$`100m`
  ocean_temp_array_cr_cc[, 3] <- projection_cc_cr$`400m`
  ocean_temp_array_cr_cc[, 4] <- projection_cc_cr$`600m`
  ocean_temp_array_cr_cc[, 5] <- projection_cc_cr$BOT

  # n_pp_arrays - keep these the same

  row_pick = 1:nrow(projection_nocc_tbgb) # row to pick
  row_pick[row_pick > 601] = 542 + (row_pick[row_pick > 601]-2) %% 60 # take modulo for repeating last five years
  check = projection_nocc_tbgb %>% select(date, date_adjust) %>%
    mutate(row = row_number(), date_loop = row_pick)
  check %>% filter(date >= '2006-01-01' & date <= '2020-12-01') %>% print(n = Inf) # check 2006 - 2020

  ## TBGB

  proj_n_pp_array_rescaled_tbgb <- array(NA, dim = c(nrow(projection_nocc_tbgb), ncol(n_pp_array_rescaled_tbgb)), dimnames = list(time = projection_nocc_tbgb$date_adjust, w = signif(full_x_tbgb, 3))) # setup an empty array for the corrected plankton

  for (j in 1:nrow(projection_nocc_tbgb)) {
    proj_n_pp_array_rescaled_tbgb[j, ] <- n_pp_array_rescaled_tbgb[row_pick[j], ] # looping 2006-2010 for future
  }

  ## CR

  proj_n_pp_array_rescaled_cr <- array(NA, dim = c(nrow(projection_nocc_cr), ncol(n_pp_array_rescaled_cr)), dimnames = list(time = projection_nocc_cr$date_adjust, w = signif(full_x_cr, 3))) # setup an empty array for the corrected plankton

  for (j in 1:nrow(projection_nocc_cr)) {
    proj_n_pp_array_rescaled_cr[j, ] <- n_pp_array_rescaled_cr[row_pick[j], ] # looping 2006-2010 for future
  }

  # Let us build fishing effort arrays - constant and using historical data

  # Set up fishing arrays - fixed at 0.2

  tbgb_effort_array <- array(0, dim = c(1201, nrow(tbgb_therm@gear_params)),
                             dimnames = list(time = 0:1200, gear = tbgb_therm@gear_params$gear)) # set up size
  tbgb_effort_array[602:1201, ] = 0.2 # whack in effort after 2010

  cr_effort_array <- array(0, dim = c(1201, nrow(cr_therm@gear_params)),
                           dimnames = list(time = 0:1200, gear = cr_therm@gear_params$gear)) # set up size
  cr_effort_array[602:1201, ] = 0.2 # whack in effort after 2010

  # Set up fishing arrays - using historical data

  # TBGB

  # Try and tailor it to monthly timesteps, 1961-2010

  tbgb_historical_effort = as(read.csv("../Alice/data/f_history_with_buffer.csv", row.names = 1), "matrix")
  gear_names = colnames(tbgb_historical_effort) # record original names
  gear_names_tib = tibble(gear_names) %>%
    mutate(gear_names_mod = case_when(gear_names == 'DEM' ~ 'DPI',
                                      gear_names == 'BAR' ~ 'PFM',
                                      gear_names == 'CEP' ~ 'ASQ',
                                      gear_names == 'ELI' ~ 'RSK',
                                      gear_names == 'ELP' ~ 'TSH',
                                      gear_names == 'PFL' ~ 'TRV',
                                      gear_names == 'PFS' ~ 'PIL',
                                      gear_names == 'RFI' ~ 'BTF',
                                      T ~ gear_names))

  fishing_years = as.numeric(rownames(tbgb_historical_effort))
  tbgb_historical_effort = tbgb_historical_effort[fishing_years > 1960 & fishing_years <= 2010, ]
  tbgb_historical_effort = tbgb_historical_effort[c(1, rep(1:nrow(tbgb_historical_effort), each = 12)), ] # make it monthly, add extra row at top

  # Now extend to 2060, looping 2006-10 fishing effort
  tbgb_historical_effort = tbgb_historical_effort[c(1:601, rep(542:601, times = 10)), ] # loop 2006-10 ten times
  dimnames(tbgb_historical_effort) = list(time = 1:nrow(tbgb_historical_effort)-1, gear = gear_names_tib$gear_names_mod) # make time 0-1200, fix species names

  # CR

  # Try and tailor it to monthly timesteps, 1961-2010

  cr_historical_effort = as(read.csv("../Samik/Chatham Rise/data/CR_effort.csv", row.names = 1), "matrix")
  fishing_years = as.numeric(rownames(cr_historical_effort))
  cr_historical_effort = cr_historical_effort[fishing_years > 1960 & fishing_years <= 2010, ]
  cr_historical_effort = cr_historical_effort[c(1, rep(1:nrow(cr_historical_effort), each = 12)), ] # make it monthly, add extra row at top

  # Now extend to 2060, looping 2006-10 fishing effort
  cr_historical_effort = cr_historical_effort[c(1:601, rep(542:601, times = 10)), ] # loop 2006-10 ten times
  dimnames(cr_historical_effort) = list(time = 1:nrow(cr_historical_effort)-1, gear = colnames(cr_historical_effort)) # make time 0-1200, species names already OK


  # Now let us run the eight simulations

  tbgb_proj_nocc = upgradeTherParams(ss_tbgb_therm,
                                     temp_min = ss_tbgb_therm@species_params$temp_min_adjust,
                                     temp_max = ss_tbgb_therm@species_params$temp_max_adjust,
                                     ocean_temp_array = ocean_temp_array_tbgb_nocc,
                                     n_pp_array = proj_n_pp_array_rescaled_tbgb,
                                     vertical_migration_array = vertical_migration_array_tbgb,
                                     exposure_array = exposure_array_tbgb,
                                     aerobic_effect = TRUE,
                                     metabolism_effect = TRUE)
  initialNResource(tbgb_proj_nocc) = 10^(other_params(tbgb_proj_nocc)$n_pp_array[1, ])/tbgb_proj_nocc@dw_full # set initial resource (necessary?)

  tbgb_proj_cc = upgradeTherParams(ss_tbgb_therm,
                                   temp_min = ss_tbgb_therm@species_params$temp_min_adjust,
                                   temp_max = ss_tbgb_therm@species_params$temp_max_adjust,
                                   ocean_temp_array = ocean_temp_array_tbgb_cc,
                                   n_pp_array = proj_n_pp_array_rescaled_tbgb,
                                   vertical_migration_array = vertical_migration_array_tbgb,
                                   exposure_array = exposure_array_tbgb,
                                   aerobic_effect = TRUE,
                                   metabolism_effect = TRUE)
  initialNResource(tbgb_proj_cc) = 10^(other_params(tbgb_proj_cc)$n_pp_array[1, ])/tbgb_proj_cc@dw_full # set initial resource (necessary?)

  cr_proj_nocc = upgradeTherParams(ss_cr_therm,
                                   temp_min = ss_cr_therm@species_params$temp_min_adjust,
                                   temp_max = ss_cr_therm@species_params$temp_max_adjust,
                                   ocean_temp_array = ocean_temp_array_cr_nocc,
                                   n_pp_array = proj_n_pp_array_rescaled_cr,
                                   vertical_migration_array = vertical_migration_array_cr,
                                   exposure_array = exposure_array_cr,
                                   aerobic_effect = TRUE,
                                   metabolism_effect = TRUE)
  initialNResource(cr_proj_nocc) = 10^(other_params(cr_proj_nocc)$n_pp_array[1, ])/cr_proj_nocc@dw_full # set initial resource (necessary?)

  cr_proj_cc = upgradeTherParams(ss_cr_therm,
                                 temp_min = ss_cr_therm@species_params$temp_min_adjust,
                                 temp_max = ss_cr_therm@species_params$temp_max_adjust,
                                 ocean_temp_array = ocean_temp_array_cr_cc,
                                 n_pp_array = proj_n_pp_array_rescaled_cr,
                                 vertical_migration_array = vertical_migration_array_cr,
                                 exposure_array = exposure_array_cr,
                                 aerobic_effect = TRUE,
                                 metabolism_effect = TRUE)
  initialNResource(cr_proj_cc) = 10^(other_params(cr_proj_cc)$n_pp_array[1, ])/cr_proj_cc@dw_full # set initial resource (necessary?)


  # TBGB

  sim_tbgb_baseline = project(tbgb_proj_nocc, t_max = 1200)
  sim_tbgb_fishing = project(tbgb_proj_nocc, t_max = 1200, effort = tbgb_historical_effort)
  sim_tbgb_cc = project(tbgb_proj_cc, t_max = 1200)
  sim_tbgb_cc_fishing = project(tbgb_proj_cc, t_max = 1200, effort = tbgb_historical_effort)

  # CR

  sim_cr_baseline = project(cr_proj_nocc, t_max = 1200)
  sim_cr_fishing = project(cr_proj_nocc, t_max = 1200, effort = cr_historical_effort)
  sim_cr_cc = project(cr_proj_cc, t_max = 1200)
  sim_cr_cc_fishing = project(cr_proj_cc, t_max = 1200, effort = cr_historical_effort)


  # Saving steady states and projection runs

  save(tbgb_model, cr_model,
       tbgb_proj_nocc, tbgb_proj_cc, cr_proj_nocc, cr_proj_cc,
       ss_tbgb_therm, ss_cr_therm, dates_vec,
       sim_tbgb_baseline, sim_tbgb_fishing, sim_tbgb_cc, sim_tbgb_cc_fishing,
       sim_cr_baseline, sim_cr_fishing, sim_cr_cc, sim_cr_cc_fishing,
       file = "files_for_paper.RData")

}


# Load in files saved from previous run ####

load('files_for_paper.RData') # comment if you don't need to load files
all_sims_tbgb = list(sim_tbgb_baseline, sim_tbgb_cc, sim_tbgb_fishing, sim_tbgb_cc_fishing)
all_sims_cr = list(sim_cr_baseline, sim_cr_cc, sim_cr_fishing, sim_cr_cc_fishing)
tbgb_biomass = pull_out_info(all_sims_tbgb, ecosystem_name = 'Tasman and Golden Bay')
cr_biomass = pull_out_info(all_sims_cr, ecosystem_name = 'Chatham Rise')

# Finding same species in both systems

same_species = tbgb_proj_cc@species_params %>%
  inner_join(cr_proj_cc@species_params, by = 'species') %>%
  pull(species) # find overlapping species

## Plots ####

## (1) Community biomass relative to baseline

# Historical

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species %>% filter(run %in% c('Baseline', 'With fishing')),
                            title_pick = 'Tasman and Golden Bay community biomass',
                            max_date = '2011-01-01', smooth_lines = T)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species %>% filter(run %in% c('Baseline', 'With fishing')),
                            title_pick = 'Chatham Rise community biomass',
                            max_date = '2011-01-01', smooth_lines = T)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Community biomass historical.png", both, width = 12, height = 6)

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species %>% filter(run %in% c('Baseline', 'With fishing')),
                            title_pick = 'Tasman and Golden Bay community biomass relative to baseline',
                            max_date = '2010-12-01', proportion = T, smooth_lines = T,
                            min_y = 0.75, max_y = 1.05)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species %>% filter(run %in% c('Baseline', 'With fishing')),
                            title_pick = 'Chatham Rise community biomass relative to baseline',
                            max_date = '2011-01-01', proportion = T, smooth_lines = T,
                            min_y = 0.7, max_y = 1.05)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Community biomass historical relative.png", both, width = 12, height = 6)

# Future

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species,
                            title_pick = 'Tasman and Golden Bay community biomass',
                            min_date = '2000-01-01', smooth_lines = T)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species,
                            title_pick = 'Chatham Rise community biomass',
                            min_date = '2000-01-01', smooth_lines = T)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Community biomass future.png", both, width = 12, height = 6)

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species,
                            title_pick = 'Tasman and Golden Bay community biomass relative to baseline',
                            min_date = '2000-01-01', proportion = T, smooth_lines = T,
                            min_y = 0.7, max_y = 1.25)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species,
                            title_pick = 'Chatham Rise community biomass relative to baseline',
                            min_date = '2000-01-01', proportion = T, smooth_lines = T,
                            min_y = 0.7, max_y = 1.25)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Community biomass future relative.png", both, width = 12, height = 6)

## (2) Species biomass relative to baseline

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species, select_group = 'species', min_date = '2000-01-01',
                            title_pick = 'Tasman and Golden Bay species biomass relative to baseline', proportion = T, smooth_lines = T) +
  scale_x_date(breaks = seq(as.Date("2000-01-01"), as.Date("2060-01-01"), by = "20 years"), date_labels = "%Y")
ggsave(filename = "paper/PNGs/TBGB species biomass.png", g1, width = 16, height = 12)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species, select_group = 'species', min_date = '2000-01-01',
                            title_pick = 'Chatham Rise species biomass relative to baseline', proportion = T, smooth_lines = T) +
  scale_x_date(breaks = seq(as.Date("2000-01-01"), as.Date("2060-01-01"), by = "20 years"), date_labels = "%Y")
ggsave(filename = "paper/PNGs/CR species biomass.png", g2, width = 16, height = 14)


## (3) Weight group biomass relative to baseline

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_weight_group, select_group = 'weight_group', min_date = '2006-01-01', title_pick = 'Tasman and Golden Bay weight group biomass relative to baseline', proportion = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/TBGB weight group biomass.png", g1, width = 12, height = 16)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_weight_group, select_group = 'weight_group', min_date = '2006-01-01', title_pick = 'Chatham Rise species biomass relative to baseline', proportion = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/CR weight group biomass.png", g2, width = 12, height = 16)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both weight group biomass.png", both, width = 16, height = 16)


## (4) Biomass range plots by species

g1 = plot_biomass_range(sim_biomass = tbgb_biomass$by_species, min_date = '2056-01-01', title_pick = 'Tasman and Golden Bay species biomass')
g2 = plot_biomass_range(sim_biomass = cr_biomass$by_species, min_date = '2056-01-01', title_pick = 'Chatham Rise species biomass')
both = grid.arrange(g1, g2, nrow = 2)
ggsave(filename = "paper/PNGs/Both species biomass range.png", both, width = 12, height = 10)

g1 = plot_biomass_range(sim_biomass = tbgb_biomass$by_species, min_date = '2056-01-01',
                        title_pick = 'Tasman and Golden Bay species biomass relative to baseline',
                        proportion = T)
g2 = plot_biomass_range(sim_biomass = cr_biomass$by_species, min_date = '2056-01-01',
                        title_pick = 'Chatham Rise species biomass relative to baseline',
                        proportion = T)
both = grid.arrange(g1, g2, nrow = 2)
ggsave(filename = "paper/PNGs/Both species biomass range relative.png", both, width = 12, height = 10)

## (6) Biomass range plots by weight group

g1 = plot_biomass_range(sim_biomass = tbgb_biomass$by_weight_group, min_date = '2056-01-01', title_pick = 'Tasman and Golden Bay weight group biomass')
g2 = plot_biomass_range(sim_biomass = cr_biomass$by_weight_group, min_date = '2056-01-01', title_pick = 'Chatham Rise species biomass')
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both weight group biomass range.png", both, width = 12, height = 6)

g1 = plot_biomass_range(sim_biomass = tbgb_biomass$by_weight_group, min_date = '2056-01-01', title_pick = 'Tasman and Golden Bay weight group biomass', proportion = T, min_y = 0.5, max_y = 6)
g2 = plot_biomass_range(sim_biomass = cr_biomass$by_weight_group, min_date = '2056-01-01', title_pick = 'Chatham Rise species biomass', proportion = T, min_y = 0.5, max_y = 6)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both weight group biomass range relative.png", both, width = 12, height = 6)


# (7) Biomass range plots by species asymptotic body size

g1 = plot_biomass_range_species_winf(sim_biomass = tbgb_biomass$by_species, min_date = '2056-01-01', title_pick = 'Tasman and Golden Bay species biomass')
g2 = plot_biomass_range_species_winf(sim_biomass = cr_biomass$by_species, min_date = '2056-01-01', title_pick = 'Chatham Rise species biomass')
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both asymptotic weight biomass range.png", both, width = 12, height = 6)

g1 = plot_biomass_range_species_winf(sim_biomass = tbgb_biomass$by_species, min_date = '2056-01-01', title_pick = 'Tasman and Golden Bay species biomass relative to baseline', proportion = T)
g2 = plot_biomass_range_species_winf(sim_biomass = cr_biomass$by_species, min_date = '2056-01-01', title_pick = 'Chatham Rise species biomass relative to baseline', proportion = T)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both asymptotic weight biomass range relative.png", both, width = 12, height = 6)

# (8) # Biomass for same species from both systems

g1 = plot_compare_same_species(sim_biomass = tbgb_biomass$by_species, title_pick = 'Tasman and Golden Bay same species biomass', min_date = '2000-01-01', proportion = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/TBGB same species biomass relative.png", g1, width = 12, height = 8)
g2 = plot_compare_same_species(sim_biomass = cr_biomass$by_species, title_pick = 'Chatham Rise same species biomass', min_date = '2000-01-01', proportion = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/CR same species biomass relative.png", g2, width = 12, height = 8)

g = plot_compare_same_species2(sim_biomass = bind_rows(tbgb_biomass$by_species, cr_biomass$by_species), title_pick = 'Comparing species present in both ecosystems', min_date = '2000-01-01', proportion = T, smooth_lines = T)
p = grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], ncol = 1)
ggsave(filename = "paper/PNGs/Both same species biomass relative 2.png", p, width = 12, height = 18)

# (9) Community yields

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species %>% filter(run %in% c('With fishing', 'With climate change and fishing')),
                            title_pick = 'Tasman and Golden Bay community yield',
                            min_date = '2000-01-01', plot_yield = T, smooth_lines = T)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species %>% filter(run %in% c('With fishing', 'With climate change and fishing')),
                            title_pick = 'Chatham Rise community yield',
                            min_date = '2000-01-01', plot_yield = T, smooth_lines = T)
both = grid.arrange(g1, g2, nrow = 1)
ggsave(filename = "paper/PNGs/Both community yield.png", both, width = 12, height = 6)

# (10) Species yields - confine to highest catches from 2006-10

load('../Alice/data/storeCatchByFleet')
tbgb_catches = t(storeCatchByStock) %>% as_tibble() %>%
  setNames(fishedCodes) %>%
  mutate(year = 1900:2012) %>% # taken from 'f_history.txt'
  relocate(year) %>%
  rename('DPI' = 'DEM', 'PFM' = 'BAR', 'ASQ' = 'CEP', 'RSK' = 'ELI',
         'TSH' = 'ELP', 'TRV' = 'PFL', 'PIL' = 'PFS', 'BTF' = 'RFI')
tbgb_most_caught = tbgb_catches %>% filter(year >= 2006, year <= 2010) %>%
  select(-year) %>%
  colSums() %>%
  sort(decreasing = T) # HOK, ORH, SSO
tbgb_most_caught # sorted by catch

cr_catches = as(read.csv("../Samik/Chatham Rise/data/CR_catches.csv", row.names = 1), "matrix")
cr_catches = cbind(cr_catches, year = rownames(cr_catches))
cr_catches = cr_catches %>% as_tibble() %>%
  mutate_if(is.character, as.numeric) %>%
  relocate(year)
cr_most_caught = cr_catches %>% filter(year >= 2006, year <= 2010) %>%
  select(-year) %>%
  colSums() %>%
  sort(decreasing = T) # HOK, ORH, SSO
cr_most_caught # sorted by catch

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species %>%
                              mutate(species = factor(species, levels = names(tbgb_most_caught))) %>%
                              filter(species %in% names(tbgb_most_caught)[1:6], !(run %in% c('Baseline', 'With climate change'))),
                            min_date = '2000-01-01', title_pick = 'Tasman and Golden Bay species yield', select_group = 'species',
                            plot_yield = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/TBGB most caught yield.png", g1, width = 12, height = 6)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species %>%
                              mutate(species = factor(species, levels = names(cr_most_caught))) %>%
                              filter(species %in% names(cr_most_caught)[1:6], !(run %in% c('Baseline', 'With climate change'))),
                            min_date = '2000-01-01', title_pick = 'Chatham Rise species yield', select_group = 'species',
                            plot_yield = T, smooth_lines = T)
ggsave(filename = "paper/PNGs/CR most caught yield.png", g2, width = 12, height = 6)

# Look at


# (12) Correlation between maximum thermal tolerance and change in biomass

g1 = thermal_tolerance_correlation(data = tbgb_biomass$by_species,
                                    min_date = '2056-01-01',
                                    col_pick = RColorBrewer::brewer.pal(6, "Set1")[1],
                                    title_pick = 'TBGB comparing maximum thermal tolerance with biomass change')
g2 = thermal_tolerance_correlation(data = cr_biomass$by_species,
                                    min_date = '2056-01-01',
                                    col_pick = RColorBrewer::brewer.pal(6, "Set1")[2],
                                    title_pick = 'CR comparing maximum thermal tolerance with biomass change')
both = grid.arrange(g1, g2, nrow = 2)
ggsave(filename = "paper/Both biomass change fit.png", both, width = 9, height = 9)


## Supplementary figures ####

# (1) Plot TherMizer steady states

sim_after_tune = project(ss_tbgb_therm, t_max = 200)
sim_after_tune2 = project(ss_cr_therm, t_max = 200)

# TBGB

g5 = plotSpectra(sim_after_tune, total = T) +
  labs(title = 'Tasman and Golden Bay spectra') +
  theme_classic()
g6 = plotBiomassObservedVsModel(sim_after_tune, ratio = T) +
  labs(title = 'Tasman and Golden Bay biomass ratios') +
  theme_classic()

# CR

g7 = plotSpectra(sim_after_tune2, total = T) +
  labs(title = 'Chatham Rise spectra') +
  theme_classic()
g8 = plotBiomassObservedVsModel(sim_after_tune2, ratio = T) +
  labs(title = 'Chatham Rise biomass ratios') +
  theme_classic()

both = grid.arrange(g5, g6, g7, g8, nrow = 2, layout_matrix= rbind(c(1, 3), c(2, 4)))
ggsave(filename = "paper/Both steady states therMizer.png", both, width = 14, height = 10)

# Pulling out biomasses compared to mizer model (average of last fifty iterations)

ss_bio_tbgb_mizer = getBiomass(sim_repeat_tbgb) %>%
  as_tibble() %>%
  slice_tail(n = 50) %>%
  colMeans()
ss_bio_cr_mizer = getBiomass(sim_repeat_cr) %>%
  as_tibble() %>%
  slice_tail(n = 50) %>%
  colMeans()
ss_bio_tbgb_therm = getBiomass(sim_after_tune) %>%
  as_tibble() %>%
  slice_tail(n = 50) %>%
  colMeans()
ss_bio_cr_therm = getBiomass(sim_after_tune2) %>%
  as_tibble() %>%
  slice_tail(n = 50) %>%
  colMeans()

compare_tbgb_models = tibble(species = names(ss_bio_tbgb_mizer),
                             ratio = ss_bio_tbgb_therm/ss_bio_tbgb_mizer) %>%
  arrange(desc(ratio))
compare_cr_models = tibble(species = names(ss_bio_cr_mizer),
                             ratio = ss_bio_cr_therm/ss_bio_cr_mizer) %>%
  arrange(desc(ratio))

tbgb_biomass$by_species %>% filter(run == 'Baseline', species == 'SPD', year(date) == 2010) %>%
  pull(biomass)/ss_bio_tbgb_therm[names(ss_bio_tbgb_therm) == 'SPD']


## (2) Plots of species biomass through time (historical)

# TBGB

g1 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'Baseline'), title_pick = 'Tasman and Golden Bay (baseline)', max_date = '2011-01-01')
g2 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'With fishing'), title_pick = 'Tasman and Golden Bay (with fishing)', max_date = '2011-01-01')

# CR

g3 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'Baseline'), title_pick = 'Chatham Rise (baseline)', max_date = '2011-01-01')
g4 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'With fishing'), title_pick = 'Chatham Rise (with fishing)', max_date = '2011-01-01')

both = grid.arrange(g1, g3, g2, g4, nrow = 2)
ggsave(filename = "paper/PNGs/Both therMizer historical.png", both, width = 14, height = 9)


## (3) Plots of species biomass through time (projections)

# TBGB

g1 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'Baseline'), title_pick = 'Tasman and Golden Bay (baseline)', min_date = '2000-01-01')
g2 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'With climate change'), title_pick = 'Tasman and Golden Bay (with climate change)', min_date = '2000-01-01')
g3 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'With fishing'), title_pick = 'Tasman and Golden Bay (with fishing)', min_date = '2000-01-01')
g4 = plot_species_biomass(sim_biomass = tbgb_biomass$by_species %>% filter(run == 'With climate change and fishing'), title_pick = 'Tasman and Golden Bay (with climate change and fishing)', min_date = '2000-01-01')

# CR

g5 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'Baseline'), title_pick = 'Chatham Rise (baseline)', min_date = '2000-01-01')
g6 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'With climate change'), title_pick = 'Chatham Rise (with climate change)', min_date = '2000-01-01')
g7 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'With fishing'), title_pick = 'Chatham Rise (with fishing)', min_date = '2000-01-01')
g8 = plot_species_biomass(sim_biomass = cr_biomass$by_species %>% filter(run == 'With climate change and fishing'), title_pick = 'Chatham Rise (with climate change and fishing)', min_date = '2000-01-01')

both = grid.arrange(g1, g5, g2, g6, g3, g7, g4, g8, nrow = 4)
ggsave(filename = "paper/PNGs/Both therMizer projections.png", both, width = 14, height = 18)

# (4) Species yields

g1 = plot_biomass_community(sim_biomass = tbgb_biomass$by_species %>% filter(run %in% c('With fishing', 'With climate change and fishing')),
                            min_date = '2000-01-01', select_group = 'species', plot_yield = T, smooth_lines = T,
                            title_pick = 'Tasman and Golden Bay species yield') +
  scale_x_date(breaks = seq(as.Date("2000-01-01"), as.Date("2060-01-01"), by = "20 years"), date_labels = "%Y")
ggsave(filename = "paper/PNGs/TBGB species yield.png", g1, width = 16, height = 12)
g2 = plot_biomass_community(sim_biomass = cr_biomass$by_species %>% filter(run %in% c('With fishing', 'With climate change and fishing')),
                            min_date = '2000-01-01', select_group = 'species', plot_yield = T, smooth_lines = T,
                            title_pick = 'Chatham Rise species yield') +
  scale_x_date(breaks = seq(as.Date("2000-01-01"), as.Date("2060-01-01"), by = "20 years"), date_labels = "%Y")
ggsave(filename = "paper/PNGs/CR species yield.png", g2, width = 16, height = 14)

# (5) Heat map of fishing effort

tbgb_fishing = tbgb_historical_effort %>%
  as_tibble() %>%
  mutate(date = dates_vec) %>%
  relocate(date) %>%
  pivot_longer(-date, names_to = "species", values_to = "effort")

cr_fishing = cr_historical_effort %>%
  as_tibble() %>%
  mutate(date = dates_vec) %>%
  relocate(date) %>%
  pivot_longer(-date, names_to = "species", values_to = "effort")

my_breaks = c(-10, -0.1, 1, 10)

g1 = plot_fishing_effort_heat_map(tbgb_fishing, title_pick = 'Tasman and Golden Bay fishing effort')
g2 = plot_fishing_effort_heat_map(cr_fishing, title_pick = 'Chatham Rise fishing effort')
both = grid.arrange(g1, g2, ncol = 1)
ggsave(filename = "paper/PNGs/Both fishing effort.png", both, width = 9, height = 9)


