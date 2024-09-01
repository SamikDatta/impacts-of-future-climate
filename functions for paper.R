## Libraries ####

library(therMizer) # automatically loads mizer too
library(tidyverse)
library(readxl)
library(gridExtra)
library(ggrepel)


## Functions for paper ####

# Colours
run_colours <- setNames(RColorBrewer::brewer.pal(6, "Dark2")[1:4],
                        c('Baseline', 'With climate change', 'With fishing', 'With climate change and fishing')) # assign colours to each layer
all_run_names = names(run_colours)
species_palette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1")) # palette for species


# Functions

pull_out_info = function(sim_list, ecosystem_name = 'Ecosystem', run_names = all_run_names, real_dates = dates_vec, weight_breaks = c(30, 1e3)) {

  by_species = NULL
  by_weight_group = NULL
  weight_breaks_mod= c(0, weight_breaks, Inf)

  weight_breaks_names = NULL
  for (k in 1:length(weight_breaks_mod)) {
    if (k < length(weight_breaks_mod)) {
      weight_breaks_names = c(weight_breaks_names, paste0('(', weight_breaks_mod[k], ',', weight_breaks_mod[k+1], ']'))
    } else {
      weight_breaks_names = c(weight_breaks_names, 'total')
    }
  }

  # By species

  for (j in 1:length(sim_list)) {
    dummy = getBiomass(sim_list[[j]]) %>%
      as_tibble() %>%
      mutate(date = real_dates) %>%
      pivot_longer(-date, names_to = 'species', values_to = 'biomass') %>%
      group_by(date) %>%
      mutate(total_biomass = sum(biomass)) %>%
      mutate(run = run_names[j])

    temp = getYield(sim_list[[j]]) %>%
      as_tibble() %>%
      mutate(date = real_dates) %>%
      pivot_longer(-date, names_to = 'species', values_to = 'yield') %>%
      group_by(date) %>%
      mutate(total_yield = sum(yield)) %>%
      mutate(run = run_names[j])

    dummy = dummy %>% left_join(temp, by = join_by(date, species, run))

    by_species = bind_rows(by_species, dummy)

    # By weight group

    for (k in 1:(length(weight_breaks_mod))) {

      if (k < length(weight_breaks_mod)) { # do by weight group

        dummy = getBiomass(sim_list[[j]], min_w = weight_breaks_mod[k], max_w = weight_breaks_mod[k+1])

      } else { # add total

        dummy = getBiomass(sim_list[[j]])

      }

      dummy = dummy %>%
        as_tibble() %>%
        mutate(date = real_dates) %>%
        pivot_longer(-date, names_to = 'species') %>%
        group_by(date) %>%
        mutate(biomass = sum(value)) %>%
        select(-c(species, value)) %>%
        distinct() %>%
        mutate(run = run_names[j], weight_group = weight_breaks_names[k])

      by_weight_group = bind_rows(by_weight_group, dummy)

    }
  }


  # Get scaled values compared to baseline (by species and overall)

  by_species = by_species %>% mutate(run = factor(run, levels = run_names)) %>%
    group_by(date) %>%
    mutate(relative_total_biomass = total_biomass / first(total_biomass)) %>%
    group_by(date, species) %>%
    mutate(relative_biomass = biomass / first(biomass)) %>%
    ungroup() %>%
    arrange(date)

  # Get scaled values compared to baseline (by weight_group)

  by_weight_group = by_weight_group %>% mutate(run = factor(run, levels = run_names), weight_group = factor(weight_group, levels = weight_breaks_names)) %>%
    group_by(date, weight_group) %>%
    mutate(relative_biomass = biomass / first(biomass)) %>%
    ungroup() %>%
    arrange(date) %>%
    mutate(ecosystem = ecosystem_name)

  # Add in specific species parameters

  by_species = by_species %>% left_join(sim_list[[1]]@params@species_params %>% select(species, w_inf, temp_min_adjust, temp_max_adjust), by = join_by(species)) %>%
    mutate(ecosystem = ecosystem_name)

  return(list(by_species = by_species,
              by_weight_group = by_weight_group))

}

plot_temp_plus_tolerance = function(dummy, same_species, title_pick = 'Temperatures plus tolerances', include_tolerance = T, plot_mean = F, adjust_tol = F, dodge_labels = F) {

  g = ggplot(dummy, aes(x = date, y = value, colour = name))
  if (plot_mean == T) {
    g = g + geom_hline(aes(yintercept = mean, colour = name))
  } else {
    g = g + geom_line()
  }

  g = g + labs(x = "Year", y = "Temperature (°C)", title = title_pick, colour = 'Depth') +
    theme_classic() +
    scale_y_continuous(limits = c(0, 30), breaks = scales::pretty_breaks(n = 10)) +
    scale_x_date(breaks = seq(min(dummy$date), max(dummy$date), by = "5 years"), date_labels = "%Y") +
    scale_color_brewer(palette = 'Set3')# +
    # theme(axis.title.x=element_blank(),
    #       axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank())

  if (include_tolerance == T) {

    if (adjust_tol == F) {

      g = g + geom_errorbar(aes(ymin = temp_min, ymax = temp_max), width = 0.2, color = "purple4") +
        geom_errorbar(data = dummy %>% filter(species %in% same_species), aes(ymin = temp_min, ymax = temp_max), width = 0.2, color = "tomato3")

      if (dodge_labels == F) {
        g = g + geom_label(aes(x = date, y = temp_max, label = species), fill = 'white', colour = "purple4") +
          geom_label(data = dummy %>% filter(species %in% same_species), aes(x = date, y = temp_max, label = species), fill = 'white', colour = "tomato3")
      } else {
        g = g + geom_label_repel(data = dummy %>% filter(!is.na(species), !species %in% same_species, name == 'SST'), aes(x = date, y = temp_max, label = species), fill = 'white', colour = "purple4", box.padding = 0.25, max.overlaps = 30) +
          geom_label_repel(data = dummy %>% filter(species %in% same_species, name == 'SST'), aes(x = date, y = temp_max, label = species), fill = 'white', colour = "tomato3", box.padding = 0.25, max.overlaps = 30)
      }

    } else {

      g = g + geom_errorbar(aes(ymin = temp_min_adjust, ymax = temp_max_adjust), width = 0.2, color = "purple4") +
        geom_errorbar(data = dummy %>% filter(species %in% same_species), aes(ymin = temp_min_adjust, ymax = temp_max_adjust), width = 0.2, color = "tomato3")

      if (dodge_labels == F) {
        g = g + geom_label(aes(x = date, y = temp_max_adjust, label = species), fill = 'white', colour = "purple4") +
          geom_label(data = dummy %>% filter(species %in% same_species), aes(x = date, y = temp_max_adjust, label = species), fill = 'white', colour = "tomato3")
      } else {
        g = g + geom_label_repel(data = dummy %>% filter(!is.na(species), !species %in% same_species, name == 'SST'), aes(x = date, y = temp_max_adjust, label = species), fill = 'white', colour = "purple4", min.segment.length = 0, max.overlaps = 15, nudge_y = 1) +
          geom_label_repel(data = dummy %>% filter(species %in% same_species, name == 'SST'), aes(x = date, y = temp_max_adjust, label = species), fill = 'white', colour = "tomato3", min.segment.length = 0, max.overlaps = 15, nudge_y = 1)
      }
    }

  }

  print(g)
}

plot_species_biomass = function(sim_biomass, title_pick = 'Species biomass', min_date = NULL, max_date = NULL) {

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  ggplot(sim_biomass, aes(x = date, y = biomass, colour = species)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = as.numeric(as.Date('2010-12-01')), colour = 'grey51', linetype = 'dashed') +
    scale_x_date(breaks = seq(min(sim_biomass$date), max(sim_biomass$date), by = "5 years"), date_labels = "%Y") +
    scale_y_log10() +
    theme_classic() +
    labs(x = 'Year', y = 'Biomass [g]', title = title_pick, colour = 'Legend') +
    scale_colour_manual(values = species_palette(length(unique(sim_biomass$species))))

}

plot_biomass_variation = function(sim, dates_vec = dates_vec, min_date = '2005-01-01',
                                  title_pick = 'Range of biomass') {

  sim_biomass = getBiomass(sim) %>%
    as_tibble() %>%
    mutate(date = dates_vec, .before = 1) %>%
    filter(date >= min_date)

  biomass_variation = tibble(species = sim@params@species_params$species,
                             biomass_observed = sim@params@species_params$biomass_observed,
                             mean_biomass = colMeans(sim_biomass %>% select(where(is.numeric))),
                             sim_biomass %>%
                               summarise(across(everything(), min)) %>%
                               pivot_longer(-date) %>%
                               select(-c(date, name)) %>%
                               rename(lowest_biomass = value),
                             sim_biomass %>%
                               summarise(across(everything(), max)) %>%
                               pivot_longer(-date) %>%
                               select(-c(date, name)) %>%
                               rename(highest_biomass = value)) %>%
    mutate(mean_percent = mean_biomass*100/biomass_observed,
           min_percent = lowest_biomass*100/biomass_observed,
           max_percent = highest_biomass*100/biomass_observed,
           percent_range = max_percent - min_percent)

  g = ggplot(biomass_variation, aes(x = species, color = percent_range)) +
    geom_hline(yintercept = 100, colour = 'black', linetype = 'dashed') +
    geom_point(aes(y = mean_percent), size = 1.5) +
    geom_errorbar(aes(ymin = min_percent, ymax = max_percent), width = 0.5) +
    labs(x = 'Species', y = '% of observed biomass', title = title_pick,
         colour = 'Percent range') +
    theme_classic() +
    scale_color_gradient(low = 'navyblue', high = 'orangered2')

  print(g)
}

plot_compare_scenario = function(no_cc_sim, no_cc_fishing_sim, cc_sim, cc_fish_sim, dates_vec = dates_vec, min_date = '2056-01-01',
                                 col_pick = run_colours,
                                 title_pick = 'Difference in mean biomass') {

  # browser()

  no_cc_sim_biomass = getBiomass(no_cc_sim) %>%
    as_tibble() %>%
    mutate(date = dates_vec, .before = 1)

  no_cc_fishing_sim_biomass = getBiomass(no_cc_fishing_sim) %>%
    as_tibble() %>%
    mutate(date = dates_vec, .before = 1)

  cc_sim_biomass = getBiomass(cc_sim) %>%
    as_tibble() %>%
    mutate(date = dates_vec, .before = 1)

  cc_fish_sim_biomass = getBiomass(cc_fish_sim) %>%
    as_tibble() %>%
    mutate(date = dates_vec, .before = 1)

  no_cc_mean_biomass = no_cc_sim_biomass %>%
    filter(date >= min_date) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-date) %>%
    select(-date) %>%
    rename(species = name, no_cc_biomass = value)

  no_cc_fishing_mean_biomass = no_cc_fishing_sim_biomass %>%
    filter(date >= min_date) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-date) %>%
    select(-date) %>%
    rename(species = name, no_cc_fishing_biomass = value)

  cc_mean_biomass = cc_sim_biomass %>%
    filter(date >= min_date) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-date) %>%
    select(-date) %>%
    rename(species = name, cc_biomass = value)

  cc_fish_mean_biomass = cc_fish_sim_biomass %>%
    filter(date >= min_date) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-date) %>%
    select(-date) %>%
    rename(species = name, cc_fish_biomass = value)

  biomass_tibble = tibble(species = no_cc_sim@params@species_params$species,
                          biomass_observed = no_cc_sim@params@species_params$biomass_observed) %>%
    left_join(no_cc_mean_biomass, by = 'species') %>%
    left_join(no_cc_fishing_mean_biomass, by = 'species') %>%
    left_join(cc_mean_biomass, by = 'species') %>%
    left_join(cc_fish_mean_biomass, by = 'species') %>%
    mutate(no_cc_percent = no_cc_biomass*100/biomass_observed,
           no_cc_fishing_percent = no_cc_fishing_biomass*100/biomass_observed,
           cc_percent = cc_biomass*100/biomass_observed,
           cc_fish_percent = cc_fish_biomass*100/biomass_observed) %>%
    pivot_longer(c(no_cc_percent, no_cc_fishing_percent, cc_percent, cc_fish_percent)) %>%
    mutate(name = case_when(name == 'no_cc_percent' ~ 'Baseline',
                            name == 'no_cc_fishing_percent' ~ 'With fishing',
                            name == 'cc_percent' ~ 'With climate change',
                            T ~ 'With climate change and fishing'),
           name = factor(name, levels = c('Baseline', 'With fishing', 'With climate change', 'With climate change and fishing')))

  g = ggplot(biomass_tibble, aes(x = species, y = value, fill = name)) +
    geom_bar(stat='identity', position='dodge') +
    labs(x = 'Species', y = '% of observed biomass', title = title_pick, fill = 'Scenario') +
    scale_fill_manual(values = col_pick) +
    theme_classic() +
    theme(legend.position="bottom") # +
  # guides(fill = guide_legend(nrow = 1))

  print(g)
}


thermal_tolerance_correlation = function(data, min_date = NULL, max_date = NULL, col_pick = 'darkorchid4', title_pick = 'Comparing maximum thermal tolerance with biomass change') {

  if (is.null(min_date)) min_date = min(data$date)
  if (is.null(max_date)) max_date = max(data$date)

  dummy = data %>% filter(run == 'With climate change', date >= min_date, date <= max_date) %>%
    group_by(species) %>%
    mutate(mean_biomass = mean(biomass),
           mean_relative_biomass = mean(relative_biomass)) %>%
    select(species, mean_biomass, mean_relative_biomass, temp_max_adjust) %>%
    distinct()

  lin_model <- lm(mean_relative_biomass ~ temp_max_adjust, dummy,
                  weights = mean_biomass)

  g = ggplot(dummy, aes(x = temp_max_adjust, y = mean_relative_biomass)) +
    geom_point(aes(size = mean_biomass/1e6), colour = col_pick) +
    geom_hline(yintercept = 1, colour = 'grey20', linetype = 'dashed') +
    geom_smooth(method='lm', mapping = aes(weight = mean_biomass),
                colour = col_pick, fill = col_pick, alpha = 0.2) +
    geom_label_repel(aes(label = species), colour = col_pick, max.overlaps = 30) +
    theme_classic() +
    # scale_x_continuous(expand=c(0,0), limits=c(0,10)) +
    scale_y_continuous(limits = c(-0.5, 3.5)) +
    coord_cartesian(ylim=c(0, 3)) +
    labs(x = 'maximum thermal tolerance (²C)', y = 'Biomass relative from baseline',
         title = title_pick, size = 'Average biomass (t)') +
    annotate("text", x = min(dummy$temp_max_adjust), y = 3, hjust = 0,
              label = paste0("Gradient = ", round(coef(lin_model)[2], 2)))

  print(g)
}

# Updated plots

plot_biomass_community = function(sim_biomass, select_group = NULL, title_pick = 'Community biomass',
                                  min_date = NULL, max_date = NULL, proportion = F, smooth_lines = F,
                                  run_colours_pick = run_colours, min_y = NULL, max_y = NULL, plot_yield = F) {

  # browser()

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  if (proportion == T) sim_biomass = sim_biomass %>% filter(run != 'Baseline')

  # Plotting

  if (proportion == F) {
    if (is.null(select_group)) {
      if (plot_yield == F) { # plot total biomass
        g = ggplot(sim_biomass, aes(x = date, y = total_biomass, colour = run, fill = run, linetype = run)) +
          labs(y = 'Biomass [g]')
      } else { # plot total yield
        g = ggplot(sim_biomass, aes(x = date, y = total_yield, colour = run, fill = run, linetype = run)) +
          labs(y = 'Yield [g]')
      }
    } else {
      if (plot_yield == F) { # plot biomass by group
        g = ggplot(sim_biomass, aes(x = date, y = biomass, colour = run, fill = run, linetype = run)) +
          labs(y = 'Biomass [g]') +
          facet_wrap(as.formula(paste("~", select_group)), scales = 'free_y', ncol = 5)
      } else { # plot yield by group
        g = ggplot(sim_biomass, aes(x = date, y = yield, colour = run, fill = run, linetype = run)) +
          labs(y = 'Yield [g]') +
          facet_wrap(as.formula(paste("~", select_group)), scales = 'free_y')
      }
    }

  } else {
    if (is.null(select_group)) {
      g = ggplot(sim_biomass, aes(x = date, y = relative_total_biomass, colour = run, fill = run, linetype = run)) +
        geom_hline(yintercept = 1, colour = 'black', linetype = 'dotted') +
        labs(y = 'Biomass relative to baseline')
    } else {
      g = ggplot(sim_biomass, aes(x = date, y = relative_biomass, colour = run, fill = run, linetype = run)) +
        geom_hline(yintercept = 1, colour = 'black', linetype = 'dotted') +
        labs(y = 'Biomass relative to baseline') +
        facet_wrap(as.formula(paste("~", select_group)), scales = 'free_y', ncol = 5)
    }
  }

  if (smooth_lines == T) {
    g = g + geom_line(alpha = 0.5) + # make lines fainter
          geom_smooth() # add smoother
  } else {
    g = g + geom_line() # make lines more visible
  }

  g = g + labs(x = 'Year', title = title_pick, colour = 'Legend', fill = 'Legend', linetype = 'Legend') +
    geom_vline(xintercept = as.numeric(as.Date('2011-01-01')), colour = 'grey51', linetype = 'dashed') +
    scale_x_date(breaks = seq(min(sim_biomass$date), max(sim_biomass$date), by = "10 years"), date_labels = "%Y") +
    scale_y_continuous(limits = c(min_y, max_y)) +
    scale_color_manual(values = run_colours_pick) +
    scale_fill_manual(values = run_colours_pick) +
    scale_linetype_manual(values = 1:6) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.position = "bottom")

  print(g)

}

plot_biomass_range = function(sim_biomass, min_date = NULL, max_date = NULL,
                              title_pick = 'Range of biomass', proportion = F,
                              min_y = NULL, max_y = NULL, y_log = T) {

  # browser()

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  if ('weight_group' %in% names(sim_biomass)) {
    sim_biomass = sim_biomass %>%
      rename(species = weight_group)
    to_label_x = 'Weight group (g)'
  } else {
    to_label_x = 'Species'
  }

  if (proportion == F) {
    g = ggplot(sim_biomass, aes(x = species, y = biomass, colour = run, fill = run)) +
      labs(y = 'Biomass [g]')

  } else {
    g = ggplot(sim_biomass %>% filter(run != 'Baseline'),
               aes(x = species, y = relative_biomass, colour = run, fill = run)) +
      geom_hline(yintercept = 1, colour = 'black', linetype = 'dashed') +
      labs(y = 'Biomass relative to baseline')
  }

  g = g + geom_boxplot(alpha = 0.4) +
    labs(x = to_label_x, title = title_pick, colour = 'Legend', fill = 'Legend') +
    scale_y_log10(limits = c(min_y, max_y)) +
    scale_color_manual(values = run_colours) +
    scale_fill_manual(values = run_colours) +
    theme_classic() +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE),
           colour = guide_legend(nrow = 2, byrow = TRUE))

  if (y_log == F) g = g + scale_y_continuous(limits = c(min_y, max_y)) # use linear scale if needed

  print(g)

}

plot_biomass_range_species_winf = function(sim_biomass, min_date = NULL, max_date = NULL, title_pick = 'Range of biomass', proportion = F, w_inf_limits = c(2e3, 10e3)) {

  # browser()

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  # Grouping by weight

  temp = sim_biomass %>% mutate(cuts = cut(w_inf, c(0, w_inf_limits, Inf))) %>%
    group_by(date, run, cuts) %>%
    summarise(group_biomass = sum(biomass), .groups = 'keep') %>%
    group_by(date, cuts) %>%
    mutate(relative_group_biomass = group_biomass / first(group_biomass)) %>%
    ungroup()

  # Add for all weights

  to_add = sim_biomass %>% mutate(cuts = 'Total') %>%
    group_by(date, run, cuts) %>%
    summarise(group_biomass = sum(biomass), .groups = 'keep') %>%
    group_by(date, cuts) %>%
    mutate(relative_group_biomass = group_biomass / first(group_biomass)) %>%
    ungroup()

  # to_add = sim_biomass %>% select(date, run, total_biomass, relative_biomass) %>%
  #   distinct() %>%
  #   mutate(cuts = 'total') %>%
  #   rename(group_biomass = total_biomass, relative_group_biomass = relative_biomass) %>%
  #   relocate(names(temp))

  temp = temp %>% add_row(to_add) %>%
    mutate(cuts = factor(cuts, levels = c(levels(temp$cuts), 'Total'))) %>%
    arrange(date, run, cuts)

  if (proportion == F) {
    g = ggplot(temp, aes(x = cuts, y = group_biomass, colour = run, fill = run)) +
      labs(y = 'Biomass [g]')
  } else {
    g = ggplot(temp %>% filter(run != 'Baseline'), aes(x = cuts, y = relative_group_biomass, colour = run, fill = run)) +
      labs(y = 'Biomass relative to baseline') +
      geom_hline(yintercept = 1, colour = 'black', linetype = 'dashed')
  }

  g = g + geom_boxplot(alpha = 0.4) +
    labs(x = 'Asymptotic weight (g)', y = 'Biomass [g]', title = title_pick, colour = 'Legend', fill = 'Legend') +
    scale_y_log10() +
    scale_fill_manual(values = run_colours) +
    scale_colour_manual(values = run_colours) +
    theme_classic() +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE),
           colour = guide_legend(nrow = 2, byrow = TRUE))

  print(g)

}

plot_compare_same_species = function(sim_biomass, title_pick = 'Species biomass', species_pick = same_species, min_date = NULL, max_date = NULL, proportion = F, smooth_lines = F, run_colours_pick = run_colours, min_y = NULL, max_y = NULL) {

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  if (proportion == T) sim_biomass = sim_biomass %>% filter(run != 'Baseline')

  # Filter to species in both ecosystems

  sim_biomass = sim_biomass %>% filter(species %in% species_pick) %>%
    mutate(ecosystem = factor(ecosystem, levels = unique(sim_biomass$ecosystem))) # put TBGB first in legend

  # Plotting

  if (proportion == F) {
    g = ggplot(sim_biomass, aes(x = date, y = biomass, colour = run, linetype = ecosystem)) +
      labs(y = 'Biomass [g]')
  } else {
    g = ggplot(sim_biomass, aes(x = date, y = relative_biomass, colour = run, linetype = ecosystem)) +
      geom_hline(yintercept = 1, colour = 'black', linetype = 'dotted') +
      labs(y = 'Biomass relative to baseline')
  }

  if (smooth_lines == T) {
    g = g + geom_line(alpha = 0.5) + # make lines fainter
      geom_smooth() # add smoother
  } else {
    g = g + geom_line() # make lines more visible
  }

  g = g + facet_wrap(~species, scales = 'free_y') +
    labs(x = 'Year', title = title_pick, colour = 'Run', linetype = 'Ecosystem') +
    geom_vline(xintercept = as.numeric(as.Date('2011-01-01')), colour = 'grey51', linetype = 'dashed') +
    scale_x_date(breaks = seq(min(sim_biomass$date), max(sim_biomass$date), by = "10 years"), date_labels = "%Y") +
    scale_color_manual(values = run_colours_pick) +
    theme_classic() +
    theme(legend.position = "bottom")

  print(g)

}

plot_compare_same_species2 = function(sim_biomass, title_pick = 'Species biomass', species_pick = same_species, min_date = NULL, max_date = NULL, proportion = F, smooth_lines = F, run_colours_pick = run_colours, min_y = NULL, max_y = NULL) {

  if (!is.null(min_date)) sim_biomass = sim_biomass %>% filter(date >= min_date)
  if (!is.null(max_date)) sim_biomass = sim_biomass %>% filter(date <= max_date)

  if (proportion == T) sim_biomass = sim_biomass %>% filter(run != 'Baseline')

  # Filter to species in both ecosystems

  sim_biomass = sim_biomass %>% filter(species %in% species_pick) %>%
    mutate(ecosystem = factor(ecosystem, levels = unique(sim_biomass$ecosystem))) # put TBGB first in legend

  # Plotting - do in loop

  p = list()

  for (j in 1:length(same_species)) {

    dummy = sim_biomass %>% filter(species == same_species[j])

    if (proportion == F) {
      g = ggplot(dummy, aes(x = date, y = biomass, colour = run)) +
        labs(y = 'Biomass [g]')

    } else {
      g = ggplot(dummy, aes(x = date, y = relative_biomass, colour = run)) +
        geom_hline(yintercept = 1, colour = 'black', linetype = 'dotted') +
        labs(y = 'Biomass relative to baseline')
    }

    if (smooth_lines == T) {
      g = g + geom_line(alpha = 0.5) + # make lines fainter
        geom_smooth() # add smoother
    } else {
      g = g + geom_line() # make lines more visible
    }

    g = g + facet_wrap(~species + ecosystem, scales = 'free_x', nrow = 1) +
      labs(x = 'Year', colour = 'Run', linetype = 'Ecosystem') +
      geom_vline(xintercept = as.numeric(as.Date('2011-01-01')), colour = 'grey51', linetype = 'dashed') +
      scale_x_date(breaks = seq(min(sim_biomass$date), max(sim_biomass$date), by = "10 years"), date_labels = "%Y") +
      scale_color_manual(values = run_colours_pick) +
      theme_classic() +
      theme(legend.position = "bottom") +
      annotate("segment", x=min(dummy$date), xend=max(dummy$date), y=-Inf, yend=-Inf) +
      annotate("segment", x=min(dummy$date), xend=min(dummy$date), y=-Inf, yend=Inf)

    if (j == 1) g = g + labs(title = title_pick)

    p[[j]] = g # assign into list
  }

  print(p)

}

plot_fishing_effort_heat_map = function(effort_data, title_pick = 'Fishing effort') {

  g = ggplot(effort_data %>% filter(date <= '2010-12-01'), aes(date, species)) +
    geom_tile(aes(fill = log10(effort))) +
    labs(x = 'Year', y = 'Species', title = title_pick, fill = 'log10(Effort)') +
    scale_x_date(breaks = seq(min(effort_data$date), max(effort_data$date), by = "10 years"), date_labels = "%Y") +
    scale_fill_distiller(palette = 'YlOrRd', direction = 1, na.value = 'grey67') +
    theme_classic() +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 10))

  print(g)

}
