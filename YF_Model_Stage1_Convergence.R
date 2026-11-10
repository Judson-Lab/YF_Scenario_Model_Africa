# Should not be run directly; run through YF_Model_Ghana or YF_Model_South_Sudan

# Convergence Analysis

# Stage 2 sample-size convergence analysis

set.seed(6232932)

# First-stage likelihood weights
nll <- res_kept$nll

w <- exp(-(nll - min(nll, na.rm = TRUE)))
w <- w / sum(w)

# Effective sample size of retained first-stage distribution
ess <- 1 / sum(w^2)

# Stage-2 scenario outcome
severe_cases <- res_kept$rho * res_kept$total_Hsymp

# Remove invalid values if necessary
ok <- is.finite(severe_cases) &
      is.finite(w) &
      w > 0

x <- severe_cases[ok]
w_use <- w[ok]

w_use <- w_use / sum(w_use)

# ------------------------------------------------------------
# Full weighted result
#
# This is the reference value toward which weighted
# resampling should converge.
# ------------------------------------------------------------

target <- weighted.mean(
  x,
  w = w_use
)


# Candidate numbers of parameter sets propagated to stage 2
sample_sizes <- c(
  250,
  500,
  1000,
  1500,
  2000,
  2500,
  3000,
  3500,
  4400
)

n_rep <- 50000


# ------------------------------------------------------------
# Weighted resampling
# ------------------------------------------------------------

resample_mean <- function(M, x, w) {

  idx <- sample(
    seq_along(x),
    size = M,
    replace = TRUE,
    prob = w
  )

  mean(x[idx])
}


estimates <- lapply(
  sample_sizes,
  function(M) {
    replicate(
      n_rep,
      resample_mean(
        M = M,
        x = x,
        w = w_use
      )
    )
  }
)


# ------------------------------------------------------------
# Summarize convergence
# ------------------------------------------------------------

convergence <- do.call(
  rbind,
  lapply(seq_along(sample_sizes), function(i) {

    vals <- estimates[[i]]

    mcse <- sd(vals)

    data.frame(
      sample_size = sample_sizes[i],

      mean_estimate = mean(vals),

      target = target,

      mcse = mcse,

      relative_mcse_pct =
        100 * mcse / target,

      bias =
        mean(vals) - target,

      relative_bias_pct =
        100 * (mean(vals) - target) / target,

      lower_95 =
        unname(quantile(vals, 0.025)),

      upper_95 =
        unname(quantile(vals, 0.975))
    )
  })
)


cat(
  "Number of retained stage-1 parameter sets:",
  nrow(res_kept),
  "\n"
)

cat(
  "Stage-1 effective sample size:",
  round(ess, 1),
  "\n"
)

cat(
  "Full weighted severe-case estimate:",
  round(target, 2),
  "\n\n"
)

# ------------------------------------------------------------
# Convergence plot
# ------------------------------------------------------------

library(ggplot2)

ggplot(
  convergence,
  aes(
    x = sample_size,
    y = relative_mcse_pct
  )
) +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = 4400,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  labs(
    x = "Stage 2 sample size",
    y = "Monte Carlo error (%)",
    title = "Stage 2 sample-size convergence",
    subtitle = paste0(
      "Severe-case outcome; stage 1 ESS = ",
      round(ess)
    )
  ) +
  theme_minimal()