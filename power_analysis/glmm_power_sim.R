# Simulation-based power analysis for the planned GLMM (Q3), per
# Dokumentation/POC Workflow - HGT Machbarkeitsstudie M. oryzae.md, Section 4.
#
# Data-independent: runnable now, using placeholder effect sizes (see the
# scenario table below). Update baseline_rate/effect_rr with observed
# values once WP1.3 delivers real HGT candidate rates, and re-run.
#
# Requires: R with glmmTMB, simr (conda env: envs/r.yaml).

library(glmmTMB)
library(simr)

set.seed(1)

simulate_scenario <- function(n_fields, n_years, n_iso_per_fy,
                               baseline_rate, effect_rr) {
  df <- expand.grid(field = factor(1:n_fields), year = factor(1:n_years))
  df$log_n_isolates <- log(n_iso_per_fy)
  df$donor     <- rbinom(nrow(df), 1, 0.5)
  df$recipient <- rbinom(nrow(df), 1, 0.5)
  df$vector    <- factor(sample(c("mChr", "Starship", "both"),
                                 nrow(df), replace = TRUE))

  beta <- c(intercept = log(baseline_rate),
            donor = 0.3, recipient = 0.2,
            vectorStarship = 0.2, vectorboth = 0.4,
            donor_vectorStarship = log(effect_rr))

  # Linear predictor incl. random intercepts (field, year) + Poisson draw.
  field_re <- rnorm(n_fields, 0, 0.4)[df$field]
  year_re  <- rnorm(n_years, 0, 0.3)[df$year]
  lp <- beta["intercept"] + beta["donor"] * df$donor +
        beta["recipient"] * df$recipient +
        ifelse(df$vector == "Starship", beta["vectorStarship"],
               ifelse(df$vector == "both", beta["vectorboth"], 0)) +
        ifelse(df$donor == 1 & df$vector == "Starship",
               beta["donor_vectorStarship"], 0) +
        field_re + year_re + df$log_n_isolates

  df$hgt_count <- rpois(nrow(df), exp(lp))
  df
}

power_for_scenario <- function(n_fields, n_years, n_iso_per_fy,
                                baseline_rate, effect_rr, nsim = 200) {
  df <- simulate_scenario(n_fields, n_years, n_iso_per_fy, baseline_rate, effect_rr)
  fit <- glmmTMB(hgt_count ~ donor * vector + recipient * vector +
                   offset(log_n_isolates) + (1 | field) + (1 | year),
                 family = poisson, data = df)
  powerSim(fit, test = fixed("donor:vectorStarship"), nsim = nsim)
}

# Scenarios from the guideline's Section 4 table. These are explicit
# placeholders, not published effect sizes - HGT event rates are unknown
# per the exposé's own risk assessment (Q3).
scenarios <- list(
  konservativ  = list(n_fields = 5,  n_years = 2, n_iso_per_fy = 10, baseline_rate = 0.05, effect_rr = 1.5),
  moderat      = list(n_fields = 10, n_years = 2, n_iso_per_fy = 20, baseline_rate = 0.08, effect_rr = 1.8),
  optimistisch = list(n_fields = 15, n_years = 3, n_iso_per_fy = 30, baseline_rate = 0.12, effect_rr = 2.0)
)

results <- lapply(scenarios, function(s) do.call(power_for_scenario, s))

cat("\nPower for the Donor x Starship interaction effect, by scenario:\n")
for (name in names(results)) {
  cat(sprintf("  %-13s %s\n", name, summary(results[[name]])$Power))
}
cat("\nFeasibility checklist (Section 4 of the guideline):\n")
cat("  - Power >= 80% for main effects in the moderate scenario?\n")
cat("  - Power >= 80% for at least one interaction effect in a realistic scenario?\n")
cat("  - If the optimistic scenario is still < 80% power for the interaction:\n")
cat("    simplify the model (interaction as exploratory only) or resize the full study.\n")
