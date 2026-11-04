
# SOURCED FROM YF_MODEL_SOUTH_SUDAN OR YF_MODEL_GHANA, NOT DIRECTLY

## build a long data.frame of trajectories
traj_df <- do.call(rbind, lapply(seq_along(traj_kept), function(i) {
  tr <- traj_kept[[i]]
  data.frame(
    traj_id  = i,
    time     = tr$time,
    Hsymp    = tr$Hsymp,
    obs_onset = tr$Hsymp * res_kept$rho[i],
    TotDeath = tr$TotDeath,
    nll      = tr$nll[1],   # assumes nll is per-trajectory
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

display.prefix = if (country == "sudan") { "South " } else { "" }

baseline.traj.plot <- ggplot() +
  geom_point(data = dat, aes(x = time, y = onset), shape = 16) +
  geom_ribbon(data = band, aes(x = time, ymin = lo, ymax = hi), alpha = 0.25) +
  geom_line(data = band, aes(x = time, y = mean), linewidth = 0.8) +
  coord_cartesian(ylim = c(0, 8)) +
  theme(plot.title = element_text(size = 10)) +
  labs(x = "day", y = "cases", title = paste0(display.prefix,tools::toTitleCase(country),": ",title)) +
  theme_minimal()


# Cumulative plot data
traj_cum <- traj_df %>%
  arrange(traj_id, time) %>%
  group_by(traj_id) %>%
  mutate(
    Csymp = cumsum(Hsymp),
    Cobs_onset = cumsum(obs_onset)
  ) %>%
  ungroup()

band_cum <- traj_cum %>%
  group_by(time) %>%
  summarise(
    mean = mean(Cobs_onset, na.rm = TRUE),
    lo   = quantile(Cobs_onset, 0.025, na.rm = TRUE),
    hi   = quantile(Cobs_onset, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

baseline.cum.plot <- ggplot() +
  geom_ribbon(
    data = band_cum,
    aes(x = time, ymin = lo, ymax = hi, fill = "95% interval"),
    alpha = 0.25
  ) +
  geom_line(
    data = band_cum,
    aes(x = time, y = mean, colour = "Model mean"),
    linewidth = 0.8
  ) +
  geom_step(
    data = dat_cum,
    aes(x = time, y = Cobs, colour = "Observed cumulative"),
    linewidth = 0.8 
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Observed cumulative" = "black",
               "Model mean" = "blue")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("95% interval" = "grey60")
  ) +
  labs(x = "day", y = "cumulative cases", title = paste0(display.prefix,tools::toTitleCase(country),": ", title)) +
  theme_minimal() + 
  theme(
    legend.position = if (country == "ghana") { "inside" } else { "none" },
    legend.position.inside = c(0.98, 0.02),      # bottom-right
    legend.justification.inside = c(1, 0),       # anchor legend's bottom-right corner
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.size = grid::unit(0.4, "lines"),
    legend.text = element_text(size = 11),
    legend.spacing.y = grid::unit(0.1, "lines"),
    legend.margin = margin(2, 2, 2, 2)
  )

# Total Deaths
death_totals <- traj_df %>%
  group_by(traj_id) %>%
  summarise(total_deaths = max(TotDeath, na.rm = TRUE), .groups = "drop")

baseline.death.plots <- ggplot(death_totals, aes(x = total_deaths)) +
  geom_histogram(bins = 30) +
  labs(x = "total deaths (end of simulation)", y = "count", title = paste0(tools::toTitleCase(country),": ",title)) +
  theme_minimal()

baseline.plot <- plot.posteriors(res_kept, posterior.params, bin.count = 50, paste0(tools::toTitleCase(country),": ", title))
  
baseline.out <- traj_results

# output plots
ggsave(filename=sprintf("Plots/%s_baseline_posteriors.png", country), 
       plot=baseline.plot, width=4,height=4,dpi=300, bg="white")
ggsave(filename=sprintf("Plots/%s_baseline_trajectories.png", country),
       plot=baseline.traj.plot, width=4,height=4,dpi=300, bg="white")
ggsave(filename=sprintf("Plots/%s_baseline_cumulative_trajectories.png", country),
       plot=baseline.cum.plot, width=4,height=4,dpi=300, bg="white")
ggsave(filename=sprintf("Plots/%s_baseline_deaths.png", country), 
       plot=baseline.death.plots, width=4,height=4,dpi=300, bg="white")
