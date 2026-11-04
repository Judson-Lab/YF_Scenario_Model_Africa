
# SOURCED FROM YF_MODEL_SOUTH_SUDAN OR YF_MODEL_GHANA, NOT DIRECTLY

#Viable Runs

# Run viable runs from the baseline scenario for other scenarios

# Now we have the N viable runs from the lhs of 10,000.  Run each set of viable
# params through the different scenarios.

# Non-heatmap version scenarios
vacc_delay_scenarios <- c(20, 60, 120) # Days until vaccination starts
vacc_coverage_scenarios <- c(0.5, 0.6, 0.8)

# Options needed to run baseline scenarios
# vacc_delay_scenarios <- c(22)
# vacc_delay_scenarios <- c(61)
# vacc_coverage_scenarios <- c(0.7)
# vacc_coverage_scenarios <- c(0.5)

 #Create a grid of all combinations
scenario_grid <- expand.grid(delay = vacc_delay_scenarios, 
                             coverage = vacc_coverage_scenarios)

scenario.count = nrow(scenario_grid)
# LHS to examine during posterior analysis
posterior.params <- c("rho", "Iv0", "v_mult")


# We want to pull from the sampled weights, rather than all the returned values.
param.values <- params_10k

good.param.count <- nrow(param.values)

# For the cumulative plots, create a cumulative data set

dat_cum <- dat %>%
  arrange(time) %>%
  mutate(Cobs = cumsum(onset)) 

# We store the trajectory/nll information in `out`
run.results = list()

for (idx in 1:scenario.count) {
  
  # For each scenario, I have to run the set of viable parameters (param.values)
  
  iter.str = sprintf("Running scenario %g/%g (%g runs)\n", idx, scenario.count, good.param.count)
  
  cat(iter.str)
  
  vals <- scenario_grid[idx,]
  params_copy <- model_pars
  
  # Grid params
  params_copy$vacc_start_day <- vals$delay
  params_copy$vacc_coverage <- vals$coverage
  
  # param.values can replace lhs_params directly
  
  # Scenario updated, run model:
  res_list <- run.parallel.params(params_copy, param.values, ncores)
  
  names(res_list$traj) <- paste0("params_", seq_along(res_list$res))
  
  run.results[[idx]] <- res_list
}

# Plots

out = list()
traj.plots = list()
cum.plots = list()
death.plots = list()
traj.bands = list()

coverage.cols = c(rep("blue",3),rep("green",3),rep("orange",3))

for (idx in 1:scenario.count) {
  
  traj <- run.results[[idx]]$traj
  res <- run.results[[idx]]$res
  
  title = sprintf("Delay %g, Coverage %.0f%%", 
                  scenario_grid[idx,]$delay, scenario_grid[idx,]$coverage * 100)
  
  ## build a long data.frame of trajectories
  traj_df <- do.call(rbind, lapply(seq_along(traj), function(i) {
    tr <- traj[[i]]
    data.frame(
      traj_id  = i,
      time     = tr$time,
      Hsymp    = tr$Hsymp,
      obs_onset = tr$Hsymp * res$rho[[i]],
      TotDeath = tr$TotDeath,
      stringsAsFactors = FALSE
    )
  }))
  
  band <- traj_df %>%
    group_by(time) %>%
    summarise(
      mean = mean(obs_onset, na.rm = TRUE),
      lo   = quantile(obs_onset, 0.025, na.rm = TRUE),
      hi   = quantile(obs_onset, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  delay_labels    <- c(" 20d", " 60d", "120d")
  coverage_labels <- c("50%", "60%", "80%")
  
  row <- (idx - 1) %/% 3 + 1
  col <- (idx - 1) %% 3 + 1
  
  traj.bands[[idx]] <- band %>%
    mutate(
      delay = factor(delay_labels[col], levels = delay_labels),
      coverage = factor(coverage_labels[row], levels = coverage_labels)
    )
  
  traj_cum <- traj_df %>%
    arrange(traj_id, time) %>%
    group_by(traj_id) %>%
    mutate(Cobs_onset = cumsum(obs_onset)) %>%
    ungroup()
  
  band_cum <- traj_cum %>%
    group_by(time) %>%
    summarise(
      mean = mean(Cobs_onset, na.rm = TRUE),
      lo   = quantile(Cobs_onset, 0.025, na.rm = TRUE),
      hi   = quantile(Cobs_onset, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  cum.plots[[idx]] <- ggplot() +
    geom_ribbon(data = band_cum, aes(x = time, ymin = lo, ymax = hi), fill = coverage.cols[idx], alpha = 0.25) +
    geom_line(data = band_cum, aes(x = time, y = mean), color = coverage.cols[idx], linewidth = 0.8) +
    labs(x = "day", y = "cumulative cases", title = title) +
    theme_minimal()
  
  # Total Deaths
  death_totals <- traj_df %>%
    group_by(traj_id) %>%
    summarise(total_deaths = max(TotDeath, na.rm = TRUE), .groups = "drop")

  death.plots[[idx]] <- ggplot(death_totals, aes(x = total_deaths)) +
    geom_histogram(bins = 30) +
    labs(x = "total deaths (end of simulation)", y = "count", title = title) +
    theme_minimal()

  out[[idx]] <- list(traj_kept=traj, res_kept=res)
}

# Create the trajectories plot; just the data right now, combine to make the stacked 3 plot
# It needs data from all scenarios,, so must be created afterwards
traj_data_frame <- bind_rows(traj.bands) %>%
  filter(time <= 200)

row_plots <- lapply(levels(traj_data_frame$coverage), function(cov_lab) {

  df_cov <- filter(traj_data_frame, coverage == cov_lab)

  ymax_row <- max(df_cov$hi, na.rm = TRUE)

  ggplot(df_cov, aes(x = time, y = mean, color = delay)) +

    ## CI ribbon only for 120-day delay
    geom_ribbon(
      data = filter(df_cov, delay == "120d"),
      aes(x = time, ymin = lo, ymax = hi, fill = "ci120", group = 1),
      inherit.aes = FALSE,
      # fill = "grey70",
      alpha = 0.30,
      color = NA
    ) +
    
    scale_fill_manual(
      values = c("ci120" = "grey70"),
      breaks = "ci120",
      labels = "120-day 95% Interval",
      name = NULL
    ) +

    geom_line(linewidth = 1.6) +

    coord_cartesian(
      xlim = c(0, 200),
      ylim = c(0, ymax_row + 2)
    ) +

    labs(
      title = paste("Vaccination coverage:", cov_lab),
      x = NULL,
      y = "daily cases",
      color = "Vaccination delay"
    ) +

    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right"
    )
})

row_plots[[3]] <- row_plots[[3]] + labs(x = "day")

## Extract the legend from the first plot
legend_grob <- cowplot::get_legend(
  row_plots[[1]] +
    theme(legend.position = "right")
)

## Remove the embedded legends from all three plots
row_plots_no_legend <- lapply(
  row_plots,
  function(p) p + theme(legend.position = "none")
)

plot.col.count = 3 # Constant column width (11x11 HEATMAP)
plot.row.count = ceiling(scenario.count / plot.col.count)

# Arrange the plots into one image
## Arrange plots in the first column and reserve a legend column
updated.trajectories <- grid.arrange(
  grobs = list(
    row_plots_no_legend[[1]], legend_grob,
    row_plots_no_legend[[2]], nullGrob(),
    row_plots_no_legend[[3]], nullGrob()
  ),
  layout_matrix = matrix(
    1:6,
    nrow = 3,
    byrow = TRUE
  ),
  widths = c(1, 0.22),
  heights = c(1, 1, 1)
)

all.cumulative.plots <- grid.arrange(grobs = cum.plots, ncol = plot.col.count, nrow=plot.row.count)
all.death.plots <- grid.arrange(grobs = death.plots, ncol = plot.col.count, nrow=plot.row.count)

ggsave(filename=sprintf("Plots/%s_trajectories.png", country),
       plot = updated.trajectories,
       width=plot.col.count * 4,
       height=plot.row.count * 4,
       dpi=300)

ggsave(filename=sprintf("Plots/%s_cumulative_trajectories.png", country),
       plot = all.cumulative.plots,
       width=plot.col.count * 4,
       height=plot.row.count * 4,
       dpi=300)

ggsave(filename=sprintf("Plots/%s_deaths.png", country),
       plot = all.death.plots,
       width=plot.col.count * 4,
       height=plot.row.count * 4,
       dpi=300)