
# SOURCED FROM YF_MODEL_SOUTH_SUDAN OR YF_MODEL_GHANA, NOT DIRECTLY


# Heatmap Only

  # We want the heatmaps to have the same results as the data tables, but that will
  # only happen if we have the same number of iterations; making the heatmap with all
  # the extraneous data is taxing for RAM, write a secondary scenario executor that
  # runs the heatmap code.  Bypasses 'run_results' and goes right to 'out'
  
vacc_delay_scenarios <- c(20, 40, 60, 80, 100, 120)
vacc_coverage_scenarios <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

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

out <- vector("list", scenario.count)

for (idx in seq_len(scenario.count)) {

  cat(sprintf(
    "Running scenario %g/%g (%g runs)\n",
    idx,
    scenario.count,
    good.param.count
  ))

  vals <- scenario_grid[idx, ]
  params_copy <- model_pars

  params_copy$vacc_start_day <- vals$delay
  params_copy$vacc_coverage  <- vals$coverage

  res_list <- run.parallel.params(
    params_copy,
    param.values,
    ncores
  )

  out[[idx]] <- list(
    res_kept = res_list$res[, c(
      "rho",
      "total_Hsymp"
    ), drop = FALSE]
  )

  rm(res_list)
  gc(FALSE)
}
#Heatmap

# We want to produce a single figure that can represent multiple data points across the simulations.

scenario_grid$cluster_id <- seq_len(nrow(scenario_grid))

# How do we interpret the values from the runs? Accept the values, the nll values
# and the method: 
#    'mean' - mean of all the values
reduce_vals <- function(vals) {
  # Match.arg will return the user specified value or the first value in the method vector above
  # This would be the "best" value, ie the one with the lowest nll
  mean(vals)
}


# Extract a numeric vector from res_kept column (one value per accepted run)
extract_from_rescol <- function(res_col) {
  force(res_col)

  function(cluster) {
    rk <- cluster$res_kept

    list(
      vals = as.numeric(rk[[res_col]] * rk$rho)
    )
  }
}

# Extract a numeric vector from each trajectory (one value per accepted run)
extract_from_traj <- function(traj_fun) {
  # traj_fun is a function that operates over the entire column, reducing it to one value.
  # So "sum", "max", "min", etc.
  force(traj_fun)
  function(cluster) {
    rk <- cluster$res_kept
    tk <- cluster$traj_kept
    
    vals <- vapply(tk, function(tr) as.numeric(traj_fun(tr)), numeric(1))
    list(vals = vals)
  }
}

make_heatmap <- function(out,
                         # value.clamp = 1300,
                         extractor,
                         delays = vacc_delay_scenarios,
                         coverages = vacc_coverage_scenarios,
                         fill_label = "value",
                         title = NULL,
                         subtitle = NULL,
                         na_fill = "grey90") {

  grid <- scenario_grid

  # Use the extractor and the reducer to compute one scalar per cluster
  grid$value <- vapply(out, function(cluster) {
    ex <- extractor(cluster)         
    reduce_vals(ex$vals)
    # min (reduce_vals(ex$vals), value.clamp)
  }, numeric(1))
  
  ok <- is.finite(grid$value) & grid$value > 0
  
  grid$text_col <- "#111111"
  
  scaled_value <- rep(NA_real_, nrow(grid))
  
  scaled_value[ok] <- scales::rescale(
    log10(grid$value[ok]),
    to = c(0, 1)
  )

  grid$text_col[ok & scaled_value < 0.55] <- "#F2F2F2"
  
  # We also want to display the actual value in the cell, not the clamped value
  grid$orig <- vapply(out, function(cluster) {
    ex <- extractor(cluster)
    floor(reduce_vals(ex$vals))
  }, numeric(1))
  
  tile_width  <- min(diff(sort(unique(grid$delay))))
  tile_height <- min(diff(sort(unique(grid$coverage))))
  
  highlight_cell <- grid %>%
    dplyr::filter(
      delay == if (country == "ghana") { 20 } else { 60 },
      coverage == if (country == "ghana") { 0.70 } else { 0.50 }
    )
  
  p <- ggplot(grid, aes(x = delay, y = coverage, fill = value)) +
    geom_tile(color = NA) +
    
    geom_tile(
      data = highlight_cell,
      aes(x = delay, y = coverage),
      width = tile_width,
      height = tile_height,
      fill = NA,
      color = "#E6E6E6",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    
    geom_text(
      aes(label = round(orig, 1), color=text_col),
      # color = "black",
      # color = text_col,
      fontface = "bold",
      size = 2
    ) + scale_color_identity() +
    viridis::scale_fill_viridis(
      option = "C",
      na.value = na_fill,
      trans="log10",
      labels = function(x) format(x, scientific = FALSE, trim = TRUE)
    ) +
    scale_x_continuous(breaks = delays) +
    scale_y_reverse(breaks = coverages) + 
    labs(x = "Vaccination start delay (days)",
         y = "Vaccination coverage",
         fill = fill_label,
         title = title,
         subtitle = subtitle) +
    theme_minimal()

  list(grid = grid, plot = p)
}

# We're scaling all returned values (in this case total_Hsymp) by rho
hsymp_m <- make_heatmap(
  out,
  extractor = extract_from_rescol("total_Hsymp"),
  fill_label = "Cumulative Cases",
)

ggsave(filename=sprintf("Plots/%s_heatmap.png", country),
       plot = hsymp_m$plot,
       width=4,
       height=4,
       dpi=300,
       bg='white')