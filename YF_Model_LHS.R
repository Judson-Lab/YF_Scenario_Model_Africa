
# SOURCED FROM YF_MODEL_SOUTH_SUDAN OR YF_MODEL_GHANA, NOT DIRECTLY


## State definition

# The skeleton model states.  We now add in an accumulator for new toxic cases
# between observations

## State names (match your model)
state_names <- c(
  "Sh","Eh","Ih","Ah","Th","Rh",
  "TotInf","TotDeath",
  "Sv","Ev","Iv",
  "Hsymp"
)

## Model Parameters

# These are the pomp model parameters from Seth's SEIR model


model_pars <- list(
  Nh = region.pop, 
  
  # LHS PARAMETERS; not used in the parameters directly
  v_mult = NA_real_, # Multiplied by the human population to get the vector population
  Iv0 = NA_real_, # Initial number of infected mosquitoes
  # rho = 0.75,   # Case reporting probability
  rho = NA_real_,   # Case reporting probability
  
  # Varying as the grid before the LHS
  vacc_coverage = NA_real_, # What percentage of the population is previously covered
                       # by the vaccine.  In the SEIR: c(0.5, 0.6, 0.8)
  vacc_start_day = NA_real_,  # The day at which the vaccination campaign starts
                        # In the SEIR model, this varies: c(20, 60, 120)
  Ih0 = 1, # Initial number of infected humans
  
  # Model Fixed Parameters
  sigma_h=1/4,  # Human incubation rate (1/sigma_h = mean latent duration)
  gamma_h=1/4,  # Recovery rate from I and A
  eta=1/2,      # Progression rate from I -> T (Toxic)
  k=1/8,        # Toxic exit rate; how long do patients stay toxic? 1/8 =
                # mean toxic duration
  delta=0.12, #Ratio of symptomatic to asymptomatic individuals
  
  ########### THETA (Case fatality rate) ##############
  # General
  theta=0.095,  # Case fatality rate of those leaving T (Carlos)
  ##################### THETA #########################
  
  b=0.4,       # Probability that a bite from an infected mosquito infects
               # a human
  c=0.5,       # Probability that biting an infected human infects the mosquito
  xi=0.5,   # Relative infectiousness of asymptomatic vs symptomatic
  
  vacc_rate_per_person = 0.079 # Daily vaccination rate per person (Carlos)
)

## The deterministic skeleton

yf_skeleton_C <- vectorfield(
  Csnippet("
    double Nh_current;
    double Nv_current;
    double lambda_h;
    double lambda_v;
    double current_vacc_flow;
    double new_infections;
    double new_deaths;
    double new_sympto;

    /* Current population totals */
    Nh_current = Sh + Eh + Ih + Ah + Th + Rh;
    Nv_current = Sv + Ev + Iv;

    /* Forces of infection */
    lambda_h = (a_t * b * Iv) / Nh_current;
    lambda_v = (a_t * c * (Ih + xi * Ah)) / Nh_current;

    /* Vaccination flow */
    if (t >= vacc_start_day + 10) {
      current_vacc_flow = vacc_rate_per_person * Sh;
    } else {
      current_vacc_flow = 0.0;
    }

    /* Incidence and death flows */
    new_infections = sigma_h * Eh;
    new_deaths     = theta * k * Th;
    new_sympto     = sigma_h * delta * Eh;

    /* Human derivatives */
    DSh = -lambda_h * Sh
          - current_vacc_flow;

    DEh = lambda_h * Sh
          - sigma_h * Eh;

    DIh = delta * sigma_h * Eh
          - eta * Ih;

    DAh = (1.0 - delta) * sigma_h * Eh
          - gamma_h * Ah;

    DTh = eta * Ih
          - k * Th;

    DRh = gamma_h * Ah
          + (1.0 - theta) * k * Th
          + current_vacc_flow;

    /* Mosquito derivatives */
    DSv = mu_v_t * Nv_current
          - lambda_v * Sv
          - mu_v_t * Sv;

    DEv = lambda_v * Sv
          - (sigma_v_t + mu_v_t) * Ev;

    DIv = sigma_v_t * Ev
          - mu_v_t * Iv;

    /* Cumulative-counter derivatives */
    DTotInf   = new_infections;
    DTotDeath = new_deaths;
    DHsymp    = new_sympto;
  ")
)

## The initial state of the skeleton

yf_rinit <- function(t0, ...) {
  p <- list(...)

  with(p, {
    Rh0 <- Nh * vacc_coverage
    Sh0 <- Nh - Rh0 - Ih0

    c(
      Sh = Sh0, Eh = 0, Ih = Ih0, Ah = 0, Th = 0, Rh = Rh0,
      TotInf = 0, TotDeath = 0,
      Sv = Nh * v_mult - Iv0, Ev = 0, Iv = Iv0, Hsymp = 0
    )
  })
}

## Pomp functions

# Create the covariates here; mu, a and sigma, that depend on time

# Important here; both the ghana and the south sudan data have the same names
# and the same format; make sure we've properly selected the correct Data above
# before running

times <- 1:365

# Create the covariates
model_t0 <- 0

covariate_times <- seq(
  from = model_t0,
  to   = max(times) + 1,
  by   = 1
)

yf_covariates <- data.frame(
  time = covariate_times,

  mu_v_t = vapply(
    covariate_times,
    mu_step_fun,
    numeric(1)
  ),

  a_t = vapply(
    covariate_times,
    a_step_fun,
    numeric(1)
  ),

  sigma_v_t = vapply(
    covariate_times,
    sigma_step_fun,
    numeric(1)
  )
)

yf_covar_table <- covariate_table(
  yf_covariates,
  times = "time",
  order = "constant"
)

dat <- data.frame(time=times, onset=onset_data)

## Measurement Model

# Now that we are introducing the case data, we need a measurement model
# Previously I was clamping the value of mu, but this was causing nll problems
# further down the line.

dmeasure <- function(onset, Hsymp, rho, ..., log) {
  if (is.na(onset)) return(if (log) 0 else 1)

  mu <- rho * Hsymp

  # exact handling when mu is 0 or negative (shouldn't be negative, but be safe)
  if (mu <= 0) {
    if (log) return(if (onset == 0) 0 else -Inf)
    else     return(if (onset == 0) 1 else 0)
  }

  dpois(onset, lambda = mu, log = log)
}

rmeasure <- function(Hsymp, rho, ...) {
  mu <- pmax(rho * Hsymp, 0)
  c(onset = rpois(1, lambda = mu))
}

## Define the pomp object

po <- pomp(
  data = dat,
  times = "time",
  t0 = model_t0,
  skeleton = yf_skeleton_C,
  rinit = yf_rinit,
  rmeasure = rmeasure,
  dmeasure = dmeasure,
  covar = yf_covar_table,
  statenames = state_names,
  paramnames = names(unlist(model_pars)),
  accumvars = "Hsymp",
  compile = TRUE
)

## LHS

# We are looking to find sets of parameters that create outbreaks of YF that match
# the epi curves of established outbreaks.  We will vary the parameters using
# a latin hypercube sampling method.  Those variables we are looking at are defined
# here

lhs.run.count <- 10000
random.seed <- 6232932 #abstract 10k seed

lhs.scenario.sample.count = 4400
lhs.scenario.sample.count.heatmap = 4400

lhs.spec <- list(
  rho     = list(type="unif", min=0.2,  max=1.0), #Reporting Rate
  Iv0     = list(type="int",  min=20,   max=500), #Initial Infected Vectors
  v_mult  = list(type="unif",  min=1,    max=8)   # Initial size of the vector population (times Human Pop)
)

## LHS Parameter generation

# Now that we have identified all the parameters we would like to modify using 
# the lhs, we have to create the new parameter sets

lhs_sample <- function(n, spec, seed=1) {
  set.seed(seed)
  pnames <- names(spec)
  U <- randomLHS(n, length(pnames))
  colnames(U) <- pnames

  out <- as.data.frame(U)
  for (nm in pnames) {
    s <- spec[[nm]]
    u <- out[[nm]]
    if (s$type == "unif") {
      out[[nm]] <- s$min + u*(s$max - s$min)
    } else if (s$type == "log") {
      out[[nm]] <- exp(log(s$min) + u*(log(s$max) - log(s$min)))
    } else if (s$type == "int") {
      out[[nm]] <- as.integer(round(s$min + u*(s$max - s$min)))
    } else {
      stop("Unknown type: ", s$type)
    }
  }
  out
}

# Now we create a function that runs the model once:
# For the first stage LHS
run.single <- function(base_params, varied_row) {
  params <- base_params
  params[names(varied_row)] <- as.list(varied_row)

  tr <- trajectory(po, params=params, format="data.frame")

  total <- sum(tr$Hsymp, na.rm=TRUE)
  peak  <- max(tr$Hsymp, na.rm=TRUE)
  tpeak <- tr$time[which.max(tr$Hsymp)]

  c(total_Hsymp=total, peak_Hsymp=peak, t_peak=tpeak)
}

# This is the second stage scenario runner
run.single.params <- function(base_params, varied_row) {
  params <- base_params
  params[names(varied_row)] <- as.list(varied_row)

  tr <- trajectory(po, params=params, format="data.frame")

  total <- sum(tr$Hsymp, na.rm=TRUE)
  peak  <- max(tr$Hsymp, na.rm=TRUE)
  tpeak <- tr$time[which.max(tr$Hsymp)]
  total_death <- sum(tr$TotDeath, na.rm=TRUE)

  list("results" = c(total_Hsymp=total, peak_Hsymp=peak, t_peak=tpeak, total_death=total_death, rho=params$rho), "tr"=tr)
}

# Running the LHS model parameter exploration

# Using r parallel for parallelization

# Now, we have the ability to generate parameters and run one simulation, put
# them together and run multiple simulations

# how many cores do we went to use?
ncores <- max(1L, parallel::detectCores() - 1L)

# This is run stage two LHS
run.parallel.params <- function(run_params, lhs_params, core.count) {
  stopifnot(is.data.frame(lhs_params) || is.matrix(lhs_params))
  idx <- seq_len(nrow(lhs_params))

  res_list <- parallel::mclapply(
    idx,
    function(i) run.single.params(run_params, lhs_params[i, ]),
    mc.cores = core.count
  )

  # Extract the per-run named vectors
  results_mat <- do.call(rbind, lapply(res_list, `[[`, "results"))
  res <- as.data.frame(results_mat, stringsAsFactors = FALSE)

  # Keep LHS parameters alongside results (assumes lhs_params is data.frame-like)

  # Extract trajectories (list of data.frames), preserving order
  traj_list <- lapply(res_list, `[[`, "tr")

  list(res = res, traj = traj_list)
}


# This is the stage 1 lhs models executor
run.lhs.parallel <- function(run_params, lhs_params, core.count) {

  # create the sequential index for the loop
  idx <- seq_len(nrow(lhs_params))
  
  res_list <- mclapply(
    idx,
    function(i) {
      run.single(run_params, lhs_params[i, ])
    },
    mc.cores = core.count
  )
  
  res <- as.data.frame(do.call(rbind, res_list))
  res <- cbind(res, lhs_params)
  
  res
}

# Similarity to baseline

# Now, we are looking for results within a certain threshold. Use the NLL 
# (negative log likelihood) function as a vehicle to rank the results

calc.nll.parallel <- function(run_params, lhs_params, core.count) {

  est <- colnames(lhs_params)

  # Complete parameter vector needed to construct traj_objfun
  init_params <- unlist(run_params)

  # Supply valid initial values for all LHS-estimated parameters.
  # These are subsequently replaced on every call to ofun().
  init_params[est] <- as.numeric(lhs_params[1, est])

  ofun <- traj_objfun(
    po,
    est = est,
    params = init_params
  )

  idx <- seq_len(nrow(lhs_params))

  nll <- unlist(mclapply(
    idx,
    function(i) {

      par_i <- setNames(
        as.numeric(lhs_params[i, est]),
        est
      )

      ofun(par = par_i)
    },
    mc.cores = core.count
  ))

  nll

}

low.runs.parallel <- function(results, run_params, lhs_params, core.count, nll.prox = 10) {
  est <- colnames(lhs_params)
  # Keep those within nll.prox NLL units of best
  keep <- results$nll <= (min(results$nll, na.rm=TRUE) + nll.prox)
  res_kept <- results[keep, ]
  table(keep)
  
  idx_kept <- which(keep)   # or keep, depending what you’re saving
  
  # Re-calculate the trajectory of the runs that were successful
  traj_kept <- mclapply(
    idx_kept,
    function(i) {
      params <- unlist(run_params)
      params[est] <- as.numeric(lhs_params[i, est])
  
      tr <- trajectory(po, params = params, format = "data.frame")
      tr$run_id <- i
      tr$nll    <- results$nll[i]
      tr
    },
    mc.cores = ncores
  )
  
  list(idx_kept = idx_kept, traj_kept = traj_kept, res_kept = res_kept)
}

# Posterior examiniation

# We look at res_kept for the best runs, and examine their parameter values
# to check whether we can narrow our hypercube parameter ranges

# Calculate bins for each parameter
bin_one_param <- function(df, bins) {
  v <- df$value
  v_max <- max(v)
  v_min <- 0

  brks <- seq(v_min, v_max, length.out = bins + 1)
  if (v_max == 0) brks <- seq(0, 1, length.out = bins + 1)  # degenerate guard

  b <- findInterval(v, brks, rightmost.closed = TRUE)
  b <- pmin(pmax(b, 1), length(brks) - 1)

  df %>%
    mutate(bin = b) %>%
    group_by(bin) %>%
    summarise(
      xmin = brks[first(bin)],
      xmax = brks[first(bin) + 1],
      count = n(),
      nll_med = median(nll, na.rm = TRUE),
      .groups = "drop"
    )
}

plot.posteriors <- function(results, params, bin.count = 50, title.text) {

  # Clean up and pivot the data
  plot_df <- results %>%
    select(any_of(c(params, "nll"))) %>%
    pivot_longer(cols = any_of(params), names_to = "param", values_to = "value") %>%
    filter(!is.na(value), !is.na(nll))
  
  # Apply the bin calculating function above to each parameter
  binned <- plot_df %>%
    group_by(param) %>%
    group_modify(~ bin_one_param(.x, bins = bin.count)) %>%
    ungroup()
  
  # Plot the histograms
  plot <- ggplot(binned) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = nll_med)) +
    facet_wrap(~ param, scales = "free_x", ncol = 1) +
    # Start at zero
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
    # Legend
    scale_fill_gradient(low = "grey10", high = "grey90", name = "median NLL") +
    theme(plot.title = element_text(size = 10)) +
    labs(x = NULL, y = "count", title = title.text)
  # print(plot)
  invisible(plot)
}

# Baseline-LHS

# We now want to produce a single, large LHS for the scenario we are mapping

# Single version params
vacc_delay_scenario <- vacc_delay
vacc_coverage_scenario <- vacc_cov

 #Create a grid of all combinations
scenario_grid <- expand.grid(#rr = relative_infectivity,
                             delay = vacc_delay_scenario,
                             coverage = vacc_coverage_scenario)

# scenario.count = nrow(scenario_grid)
# LHS to examine during posterior analysis
posterior.params <- c("rho", "Iv0", "v_mult")

# For the cumulative plots, create a cumulative data set

dat_cum <- dat %>%
  arrange(time) %>%
  mutate(Cobs = cumsum(onset)) 

iter.str = sprintf("%s Baseline LHS, Size %g", tools::toTitleCase(country), lhs.run.count)
 
cat(iter.str)
  
vals <- list(delay=vacc_delay_scenario, coverage=vacc_coverage_scenario)
params_copy <- model_pars
  
# Grid params
params_copy$vacc_start_day <- vals$delay
params_copy$vacc_coverage <- vals$coverage
  
# Sample the LHS
lhs_params <- lhs_sample(n=lhs.run.count, spec=lhs.spec, seed=random.seed)
  
# Scenario updated, run model:
res <- run.lhs.parallel(params_copy, lhs_params, ncores)
res$nll <- calc.nll.parallel(params_copy, lhs_params, ncores)
  
# Calculate within nll.prox points, and retrieve those trajectories
traj_results = low.runs.parallel(res, params_copy, lhs_params, ncores, nll.prox = 10)

traj_kept <- traj_results$traj_kept
idx_kept <- traj_results$idx_kept
res_kept <- traj_results$res_kept

names(traj_kept) <- paste0("lhs_", idx_kept)

title = sprintf("Delay %g, Coverage %.0f%%",
                vacc_delay, vacc_cov * 100)
  