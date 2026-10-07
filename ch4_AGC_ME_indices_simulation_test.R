############################################################
# CH4 manuscript simulation-based reproducibility script
# Purpose:
#   1) Generate simulated data with the same column structure as real evaluation data
#   2) Demonstrate genetic-trend adjustment
#   3) Calculate AGC with bootstrap SE and 95% CI
#   4) Calculate MEP/MEF/MEG and RMEP/RMEF/RMEG
#   5) Demonstrate a simulated daughter-performance workflow using latent methane merit
#   6) Calculate LR validation statistics
#
# Notes:
#   - This is a simulation/testing workflow, not a reproduction of confidential empirical results.
#   - Use B_AGC=1000 for the full simulation run; smaller values are for quick code checks only.
#   - RME reliability values below are explicitly labelled simulation proxies. Real analyses
#     must use the official index reliability procedure/results.
############################################################

rm(list = ls())

required_packages <- c("data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package: ", pkg)
  }
}
library(data.table)

set.seed(20260606)

# ==========================================================
# Module 0. Basic settings
# ==========================================================
# This module defines shared constants, output locations, and trait groups.
# The bootstrap count can be controlled without editing the script:
#   Sys.setenv(B_AGC = 50)       # quick test
#   Sys.setenv(B_AGC = 1000)     # full simulation run

get_bootstrap_B <- function(default = 1000L) {
  env_B <- Sys.getenv("B_AGC", unset = NA_character_)
  if (is.na(env_B) || !nzchar(env_B)) return(default)

  B <- suppressWarnings(as.integer(env_B))
  if (is.na(B) || B < 1L) {
    stop("B_AGC must be a positive integer. Current value: ", env_B)
  }
  B
}

outdir <- "CH4_simulation_test_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# For the full simulation run, use the default B_agc = 1000.
# For quick local testing, set environment variable B_AGC to 50, 200, or 300.
B_agc <- get_bootstrap_B(default = 1000L)

traits37 <- c(
  "MY", "FY", "PY", "FP", "PP",
  "UH", "LONG", "FF", "BCS", "DCE", "MCE",
  "STA", "CWI", "BDE", "RAN", "RWI", "FAN", "RLS", "RLR",
  "UDE", "USU", "FUD", "FTP", "TLE", "RUH", "RTP", "ANG",
  "OFL", "OUS", "OCS",
  "VEL", "VEM", "VEC", "VEP", "VET", "VEF", "VEG"
)
all_traits <- c("PCH4", traits37)
prod_traits <- c("MY", "FY", "PY")
functional_index <- "VEF"
global_index <- "VEG"

# ==========================================================
# Module 1. Generate simulated GEBV + GREL data
# ==========================================================
# Purpose:
#   - Create a synthetic animal-level dataset with the same column structure
#     expected from real GEBV/GREL input files.
#   - Include genetic trends by birth year so the adjustment step can be tested.
#   - Force the first n_selected_like bulls to meet reliability-like criteria so
#     the downstream selected-bull workflow has a stable test sample.

simulate_gebv_data <- function(n_bulls = 1200, n_cows = 5000, n_selected_like = 1020) {
  n <- n_bulls + n_cows
  dt <- data.table(
    animal_id = sprintf("A%06d", seq_len(n)),
    animal_type = c(rep("bull", n_bulls), rep("cow", n_cows))
  )

  dt[animal_type == "bull", birth_year := sample(1995:2017, .N, replace = TRUE)]
  dt[animal_type == "cow",  birth_year := sample(2018:2022, .N, replace = TRUE,
                                                 prob = c(0.15, 0.20, 0.35, 0.20, 0.10))]

  dt[, n_pch4_daughters := 0L]
  dt[animal_type == "bull", n_pch4_daughters := sample(0:70, .N, replace = TRUE)]
  dt[1:n_selected_like, n_pch4_daughters := sample(30:300, n_selected_like, replace = TRUE)]

  dt[, has_pch4_record := as.integer(animal_type == "cow" & runif(.N) < 0.75)]

  # Latent genetic factors
  prod <- rnorm(n)
  fat  <- rnorm(n)
  func <- rnorm(n)
  conf <- rnorm(n)
  misc <- rnorm(n)

  # Adjusted genetic values before adding long-term trend
  g_adj <- list()
  g_adj$PCH4 <- 7.0 * prod + 6.0 * fat + 2.5 * func + rnorm(n, 0, 5.5)
  # Keep the latent trend-free methane genetic merit for the simulated daughter workflow.
  # It is not used to construct RME indices and is available only because this is a simulation.
  dt[, PCH4_TRUE_BV_adj := g_adj$PCH4]

  g_adj$MY <-  8.0 * prod - 1.5 * fat + rnorm(n, 0, 5.0)
  g_adj$FY <-  5.0 * prod + 4.5 * fat + rnorm(n, 0, 4.5)
  g_adj$PY <-  6.0 * prod + 1.5 * fat + rnorm(n, 0, 4.8)
  g_adj$FP <-  7.5 * fat - 2.0 * prod + rnorm(n, 0, 4.0)
  g_adj$PP <-  4.0 * fat - 1.0 * prod + rnorm(n, 0, 4.5)

  func_traits <- c("UH", "LONG", "FF", "BCS", "DCE", "MCE")
  for (tr in func_traits) {
    g_adj[[tr]] <- 5.5 * func + rnorm(n, 0, 6.0)
  }

  conf_traits <- c("STA", "CWI", "BDE", "RAN", "RWI", "FAN", "RLS", "RLR",
                   "UDE", "USU", "FUD", "FTP", "TLE", "RUH", "RTP", "ANG",
                   "OFL", "OUS", "OCS")
  for (tr in conf_traits) {
    g_adj[[tr]] <- 4.5 * conf + 1.5 * misc + rnorm(n, 0, 6.0)
  }

  # Simulated economic indices
  g_adj$VEL <- 0.45 * g_adj$MY + 0.35 * g_adj$FY + 0.30 * g_adj$PY + rnorm(n, 0, 4.0)
  g_adj$VEM <- 0.50 * g_adj$STA + 0.30 * g_adj$CWI + rnorm(n, 0, 5.0)
  g_adj$VEC <- 0.50 * g_adj$BDE + 0.30 * g_adj$RWI + rnorm(n, 0, 5.0)
  g_adj$VEP <- 0.35 * g_adj$UDE + 0.35 * g_adj$USU + 0.30 * g_adj$OUS + rnorm(n, 0, 5.0)
  g_adj$VET <- 0.35 * g_adj$OFL + 0.35 * g_adj$OUS + 0.30 * g_adj$OCS + rnorm(n, 0, 5.0)
  g_adj$VEF <- 0.25 * g_adj$UH + 0.30 * g_adj$LONG + 0.25 * g_adj$FF + 0.20 * g_adj$BCS + rnorm(n, 0, 4.0)
  g_adj$VEG <- 0.40 * g_adj$VEL + 0.35 * g_adj$VEF + 0.15 * g_adj$VEP + 0.10 * g_adj$VET + rnorm(n, 0, 4.0)

  # Add trait-specific genetic trends by birth year
  center_year <- dt$birth_year - 2005
  for (tr in all_traits) {
    # Stronger positive trend for PCH4 and production traits, weaker/random for others
    slope <- if (tr == "PCH4") 0.45 else if (tr %in% c("MY", "FY", "PY", "VEL", "VEG")) 0.35 else runif(1, -0.10, 0.18)
    quad  <- if (tr %in% c("PCH4", "MY", "FY", "PY", "VEL", "VEG")) 0.012 else runif(1, -0.004, 0.006)
    trend <- slope * center_year + quad * center_year^2
    dt[, paste0(tr, "_GEBV") := as.numeric(g_adj[[tr]] + trend)]
  }

  # Simulate GREL. First 1,020 bulls are designed to pass selection criteria.
  for (tr in all_traits) {
    rel <- numeric(n)
    rel[seq_len(n_selected_like)] <- runif(n_selected_like, 0.55, 0.96)
    rel[(n_selected_like + 1):n_bulls] <- runif(n_bulls - n_selected_like, 0.25, 0.90)
    rel[(n_bulls + 1):n] <- runif(n_cows, 0.20, 0.80)
    dt[, paste0(tr, "_GREL") := rel]
  }

  dt[]
}

gebv <- simulate_gebv_data()
fwrite(gebv, file.path(outdir, "sim_gebv_all_animals.csv"))

# ==========================================================
# Module 2. Select highly reliable bulls
# ==========================================================
# Purpose:
#   - Mimic the manuscript selection rule for bulls with enough PCH4 daughters
#     and acceptable GREL for every trait.
#   - Save the selected reference set used for trend fitting and AGC estimates.

gebv_cols <- paste0(all_traits, "_GEBV")
grel_cols <- paste0(all_traits, "_GREL")

selected_bulls <- copy(gebv[animal_type == "bull" & n_pch4_daughters >= 30])
for (gc in grel_cols) {
  selected_bulls <- selected_bulls[get(gc) >= 0.50]
}
cat("No. of selected bulls =", nrow(selected_bulls), "\n")
fwrite(selected_bulls, file.path(outdir, "selected_bulls.csv"))

# ==========================================================
# Module 3. Genetic trend adjustment
# ==========================================================
# Purpose:
#   - Estimate birth-year trends from selected bulls.
#   - Subtract predicted trend components from every animal's trait GEBV.
#   - Save trend tables so the adjustment can be inspected outside R.

fit_trend_one <- function(ref_data, trait) {
  gebv_col <- paste0(trait, "_GEBV")
  trend_by_year <- ref_data[!is.na(get(gebv_col)) & !is.na(birth_year), .(
    n_bulls = .N,
    mean_gebv = mean(get(gebv_col), na.rm = TRUE)
  ), by = birth_year]
  fit <- lm(mean_gebv ~ birth_year + I(birth_year^2), data = trend_by_year, weights = n_bulls)
  list(trait = trait, data = trend_by_year, fit = fit)
}

apply_trend_adjustment <- function(data, trend_fits, traits) {
  out <- copy(data)
  for (tr in traits) {
    gebv_col <- paste0(tr, "_GEBV")
    adj_col <- paste0(tr, "_GEBV_adj")
    pred <- predict(trend_fits[[tr]]$fit, newdata = data.frame(birth_year = out$birth_year))
    out[, (adj_col) := get(gebv_col) - pred]
  }
  out[]
}

trend_fits <- setNames(lapply(all_traits, function(tr) fit_trend_one(selected_bulls, tr)), all_traits)
gebv_adj <- apply_trend_adjustment(gebv, trend_fits, all_traits)
selected_bulls_adj <- gebv_adj[animal_id %in% selected_bulls$animal_id]

trend_table <- rbindlist(lapply(all_traits, function(tr) {
  tmp <- copy(trend_fits[[tr]]$data)
  tmp[, trait := tr]
  tmp[, expected_trend := predict(trend_fits[[tr]]$fit, newdata = tmp)]
  tmp
}), fill = TRUE)
fwrite(trend_table, file.path(outdir, "genetic_trend_by_birth_year.csv"))

adj_trend_check <- rbindlist(lapply(all_traits, function(tr) {
  selected_bulls_adj[, .(
    n = .N,
    mean_adj = mean(get(paste0(tr, "_GEBV_adj")), na.rm = TRUE)
  ), by = birth_year][, trait := tr]
}), fill = TRUE)
fwrite(adj_trend_check, file.path(outdir, "adjusted_trend_check.csv"))

# ==========================================================
# Module 4. AGC with bootstrap SE and 95% confidence interval
# ==========================================================
# Purpose:
#   - Estimate approximate genetic correlations from GEBV correlations and
#     individual-animal reliability weights.
#   - Use non-parametric bootstrap resampling to obtain SE and percentile CI.
#   - Do NOT truncate estimates to [-1, 1]. The attenuation-correction
#     approximation can occasionally exceed the theoretical correlation bounds,
#     especially when reliabilities are simulated or approximated. Such cases are
#     flagged as diagnostics rather than silently clipped, because clipping biases
#     the bootstrap SE and confidence interval.

agc_one <- function(data, x_col, y_col, rel_x_col, rel_y_col) {
  dt <- data[!is.na(get(x_col)) & !is.na(get(y_col)) &
               !is.na(get(rel_x_col)) & !is.na(get(rel_y_col))]
  if (nrow(dt) < 3) {
    return(data.table(
      n = nrow(dt), r_gebv = NA_real_, sum_rel_x = NA_real_,
      sum_rel_y = NA_real_, sum_rel_xy = NA_real_,
      agc_factor = NA_real_, AGC = NA_real_, outside_unit_interval = NA
    ))
  }

  r_gebv <- cor(dt[[x_col]], dt[[y_col]], use = "complete.obs")
  sum_rel_x <- sum(dt[[rel_x_col]], na.rm = TRUE)
  sum_rel_y <- sum(dt[[rel_y_col]], na.rm = TRUE)
  sum_rel_xy <- sum(dt[[rel_x_col]] * dt[[rel_y_col]], na.rm = TRUE)

  if (sum_rel_x <= 0 || sum_rel_y <= 0 || sum_rel_xy <= 0) {
    return(data.table(
      n = nrow(dt), r_gebv = r_gebv, sum_rel_x = sum_rel_x,
      sum_rel_y = sum_rel_y, sum_rel_xy = sum_rel_xy,
      agc_factor = NA_real_, AGC = NA_real_, outside_unit_interval = NA
    ))
  }

  agc_factor <- sqrt(sum_rel_x * sum_rel_y) / sum_rel_xy
  agc <- agc_factor * r_gebv

  data.table(
    n = nrow(dt), r_gebv = r_gebv,
    sum_rel_x = sum_rel_x, sum_rel_y = sum_rel_y,
    sum_rel_xy = sum_rel_xy, agc_factor = agc_factor,
    AGC = agc,
    outside_unit_interval = abs(agc) > 1
  )
}

agc_bootstrap_ci <- function(data, x_col, y_col, rel_x_col, rel_y_col,
                             B = 1000, conf_level = 0.95, seed = 20260606) {
  set.seed(seed)
  dt <- data[!is.na(get(x_col)) & !is.na(get(y_col)) &
               !is.na(get(rel_x_col)) & !is.na(get(rel_y_col))]
  n <- nrow(dt)
  if (n < 3) stop("Too few complete records for bootstrap AGC: ", x_col, " vs ", y_col)

  est <- agc_one(dt, x_col, y_col, rel_x_col, rel_y_col)

  boot_agc <- replicate(B, {
    idx <- sample.int(n, size = n, replace = TRUE)
    agc_one(dt[idx], x_col, y_col, rel_x_col, rel_y_col)$AGC
  })
  boot_agc <- as.numeric(boot_agc)
  valid_boot <- is.finite(boot_agc)
  if (sum(valid_boot) < 2L) {
    stop("Too few finite bootstrap estimates for AGC: ", x_col, " vs ", y_col)
  }

  alpha <- 1 - conf_level
  ci <- as.numeric(quantile(
    boot_agc[valid_boot],
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  ))

  est[, `:=`(
    SE = sd(boot_agc[valid_boot], na.rm = TRUE),
    CI_low = ci[1],
    CI_high = ci[2],
    conf_level = conf_level,
    CI_method = "percentile_bootstrap",
    B = B,
    n_boot_finite = sum(valid_boot),
    bootstrap_fraction_outside_unit_interval = mean(abs(boot_agc[valid_boot]) > 1)
  )]
  est[]
}

agc_pch4_adj <- rbindlist(lapply(traits37, function(tr) {
  res <- agc_bootstrap_ci(
    data = selected_bulls_adj,
    x_col = "PCH4_GEBV_adj",
    y_col = paste0(tr, "_GEBV_adj"),
    rel_x_col = "PCH4_GREL",
    rel_y_col = paste0(tr, "_GREL"),
    B = B_agc
  )
  res[, `:=`(target = "PCH4", trait = tr, scale = "trend_adjusted")]
  res
}), fill = TRUE)
fwrite(agc_pch4_adj, file.path(outdir, "AGC_PCH4_vs_37traits_adjusted_with_95CI.csv"))

# ==========================================================
# Module 5. Methane efficiency indices: MEP/MEF/MEG and RMEP/RMEF/RMEG
# ==========================================================
# Purpose:
#   - Fit expected PCH4 from production, functional, and global index predictors.
#   - Define ME as expected PCH4 minus observed adjusted PCH4, so larger ME
#     indicates lower-than-expected methane emissions.
#   - Standardize ME to RME with reference cows born in 2020 that have PCH4
#     records; by construction, RME has mean 100 and SD 10 in that reference set.

fit_expected_pch4 <- function(data, predictors) {
  form <- as.formula(paste("PCH4_GEBV_adj ~", paste(paste0(predictors, "_GEBV_adj"), collapse = " + ")))
  lm(form, data = data)
}

fit_MEP <- fit_expected_pch4(selected_bulls_adj, prod_traits)
fit_MEF <- fit_expected_pch4(selected_bulls_adj, functional_index)
fit_MEG <- fit_expected_pch4(selected_bulls_adj, global_index)

coef_table <- rbindlist(list(
  data.table(index = "MEP", term = names(coef(fit_MEP)), coefficient = as.numeric(coef(fit_MEP))),
  data.table(index = "MEF", term = names(coef(fit_MEF)), coefficient = as.numeric(coef(fit_MEF))),
  data.table(index = "MEG", term = names(coef(fit_MEG)), coefficient = as.numeric(coef(fit_MEG)))
), fill = TRUE)
fwrite(coef_table, file.path(outdir, "ME_expected_PCH4_regression_coefficients.csv"))

add_ME <- function(data, fit_MEP, fit_MEF, fit_MEG) {
  out <- copy(data)
  out[, expPCH4_MEP := predict(fit_MEP, newdata = out)]
  out[, expPCH4_MEF := predict(fit_MEF, newdata = out)]
  out[, expPCH4_MEG := predict(fit_MEG, newdata = out)]

  # ME = expected PCH4 - observed PCH4.
  # Larger ME means lower-than-expected methane emission and greater methane efficiency.
  out[, MEP := expPCH4_MEP - PCH4_GEBV_adj]
  out[, MEF := expPCH4_MEF - PCH4_GEBV_adj]
  out[, MEG := expPCH4_MEG - PCH4_GEBV_adj]
  out[]
}

gebv_me <- add_ME(gebv_adj, fit_MEP, fit_MEF, fit_MEG)

ref2020 <- gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]
standardize_ME <- function(x, ref_x) {
  if (length(ref_x) < 2L || is.na(sd(ref_x, na.rm = TRUE)) || sd(ref_x, na.rm = TRUE) == 0) {
    stop("Reference group is too small or has zero variance; cannot standardize ME.")
  }
  (x - mean(ref_x, na.rm = TRUE)) / sd(ref_x, na.rm = TRUE) * 10 + 100
}

gebv_me[, RMEP := standardize_ME(MEP, ref2020$MEP)]
gebv_me[, RMEF := standardize_ME(MEF, ref2020$MEF)]
gebv_me[, RMEG := standardize_ME(MEG, ref2020$MEG)]

ref_stats <- data.table(
  index = c("MEP", "MEF", "MEG", "RMEP", "RMEF", "RMEG"),
  n_ref = nrow(ref2020),
  mean_ref = c(mean(ref2020$MEP), mean(ref2020$MEF), mean(ref2020$MEG),
               mean(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEP),
               mean(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEF),
               mean(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEG)),
  sd_ref = c(sd(ref2020$MEP), sd(ref2020$MEF), sd(ref2020$MEG),
             sd(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEP),
             sd(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEF),
             sd(gebv_me[animal_type == "cow" & birth_year == 2020 & has_pch4_record == 1]$RMEG))
)
fwrite(ref_stats, file.path(outdir, "ME_RME_standardization_check.csv"))
fwrite(gebv_me, file.path(outdir, "all_animals_with_ME_RME_indices.csv"))

selected_bulls_me <- gebv_me[animal_id %in% selected_bulls$animal_id]

# Simulation-only reliability proxies for RME indices.
# IMPORTANT: these are NOT the official index reliabilities used in the empirical study.
# They are retained only so the AGC workflow can be exercised end-to-end without
# confidential/official reliability inputs. The *_GREL_PROXY names are deliberate.
selected_bulls_me[, RMEP_GREL_PROXY := pmin(pmax(rowMeans(.SD, na.rm = TRUE), 0.001), 0.999),
                  .SDcols = c("PCH4_GREL", "MY_GREL", "FY_GREL", "PY_GREL")]
selected_bulls_me[, RMEF_GREL_PROXY := pmin(pmax(rowMeans(.SD, na.rm = TRUE), 0.001), 0.999),
                  .SDcols = c("PCH4_GREL", "VEF_GREL")]
selected_bulls_me[, RMEG_GREL_PROXY := pmin(pmax(rowMeans(.SD, na.rm = TRUE), 0.001), 0.999),
                  .SDcols = c("PCH4_GREL", "VEG_GREL")]
selected_bulls_me[, RME_reliability_source := "simulation_proxy_mean_component_GREL"]

fwrite(selected_bulls_me, file.path(outdir, "selected_bulls_with_ME_RME_indices.csv"))

# AGC between RME indices and 37 traits, including 95% CI.
agc_RME_vs_37 <- rbindlist(lapply(c("RMEP", "RMEF", "RMEG"), function(idx) {
  rbindlist(lapply(traits37, function(tr) {
    res <- agc_bootstrap_ci(
      data = selected_bulls_me,
      x_col = idx,
      y_col = paste0(tr, "_GEBV_adj"),
      rel_x_col = paste0(idx, "_GREL_PROXY"),
      rel_y_col = paste0(tr, "_GREL"),
      B = B_agc
    )
    res[, `:=`(target = idx, trait = tr, scale = "trend_adjusted")]
    res
  }), fill = TRUE)
}), fill = TRUE)
fwrite(agc_RME_vs_37, file.path(outdir, "AGC_RME_indices_vs_37traits_with_95CI.csv"))

# AGC among RME indices, including 95% CI.
agc_RME_among <- rbindlist(list(
  cbind(pair = "RMEP_RMEF", agc_bootstrap_ci(selected_bulls_me, "RMEP", "RMEF", "RMEP_GREL_PROXY", "RMEF_GREL_PROXY", B = B_agc)),
  cbind(pair = "RMEP_RMEG", agc_bootstrap_ci(selected_bulls_me, "RMEP", "RMEG", "RMEP_GREL_PROXY", "RMEG_GREL_PROXY", B = B_agc)),
  cbind(pair = "RMEF_RMEG", agc_bootstrap_ci(selected_bulls_me, "RMEF", "RMEG", "RMEF_GREL_PROXY", "RMEG_GREL_PROXY", B = B_agc))
), fill = TRUE)
fwrite(agc_RME_among, file.path(outdir, "AGC_among_RME_indices_with_95CI.csv"))

# ==========================================================
# Module 6. Simulated daughter validation
# ==========================================================
# Purpose:
#   - Create a sire-level daughter mean PCH4 validation dataset.
#   - Test whether higher RME predicts lower daughter PCH4, as intended.
#   - Summarize both continuous regressions and RME class means.

# The simulated daughter phenotype is generated from the sire's latent methane
# genetic merit, NOT from RMEP/RMEF/RMEG. This avoids building the validation
# outcome directly from the predictor being tested. It remains a simulation-only
# workflow check and must not be interpreted as independent biological validation.
val_dau <- selected_bulls_me[, .(
  sire_id = animal_id,
  n_daughters = n_pch4_daughters,
  PCH4_TRUE_BV_adj,
  RMEP, RMEF, RMEG
)]
val_dau[, daughter_mean_PCH4 :=
          345 + 0.60 * PCH4_TRUE_BV_adj +
          rnorm(.N, mean = 0, sd = 45 / sqrt(pmax(n_daughters, 1L)))]
fwrite(val_dau, file.path(outdir, "sim_daughter_pch4_by_sire.csv"))

reg_RME <- rbindlist(lapply(c("RMEP", "RMEF", "RMEG"), function(idx) {
  fit <- lm(as.formula(paste("daughter_mean_PCH4 ~", idx)), data = val_dau)
  sm <- summary(fit)
  data.table(
    index = idx,
    n_bulls = nobs(fit),
    intercept = coef(fit)[1],
    slope = coef(fit)[2],
    R2 = sm$r.squared,
    p_value = coef(sm)[2, 4]
  )
}))
fwrite(reg_RME, file.path(outdir, "daughter_PCH4_regression_on_RME.csv"))

class_breaks <- c(-Inf, 85, 95, 105, 115, Inf)
class_labels <- c("<85", "85-95", "95-105", "105-115", ">115")

class_summary <- rbindlist(lapply(c("RMEP", "RMEF", "RMEG"), function(idx) {
  tmp <- copy(val_dau)
  tmp[, RME_class := cut(get(idx), breaks = class_breaks, labels = class_labels, right = FALSE)]
  tmp[, .(
    n_bulls = .N,
    mean_daughter_PCH4 = mean(daughter_mean_PCH4),
    sd_daughter_PCH4 = sd(daughter_mean_PCH4),
    se_daughter_PCH4 = sd(daughter_mean_PCH4) / sqrt(.N)
  ), by = RME_class][, index := idx][]
}), fill = TRUE)
fwrite(class_summary, file.path(outdir, "daughter_PCH4_by_RME_class_summary.csv"))

# ==========================================================
# Module 7. Simulated LR validation statistics
# ==========================================================
# Purpose:
#   - Create partial and whole GEBV values for validation animals.
#   - Calculate LR-style accuracy, dispersion, and bias diagnostics.

lr_validation <- function(data, gebv_partial = "GEBV_partial", gebv_whole = "GEBV_whole",
                          sigma_a2, f_col = NULL) {
  dt <- data[!is.na(get(gebv_partial)) & !is.na(get(gebv_whole))]
  if (nrow(dt) < 3L) stop("Too few complete validation animals for LR statistics.")

  cov_pw <- cov(dt[[gebv_partial]], dt[[gebv_whole]], use = "complete.obs")
  var_p <- var(dt[[gebv_partial]], na.rm = TRUE)
  cor_pw <- cor(dt[[gebv_partial]], dt[[gebv_whole]], use = "complete.obs")
  f_bar <- if (!is.null(f_col) && f_col %in% names(dt)) {
    mean(dt[[f_col]], na.rm = TRUE)
  } else {
    0
  }

  denom <- (1 - f_bar) * sigma_a2
  accuracy <- if (is.finite(cov_pw) && cov_pw >= 0 && is.finite(denom) && denom > 0) {
    sqrt(cov_pw / denom)
  } else {
    NA_real_
  }
  dispersion <- if (is.finite(var_p) && var_p > 0) cov_pw / var_p else NA_real_
  bias <- mean(dt[[gebv_partial]] - dt[[gebv_whole]], na.rm = TRUE)

  # In the LR method, cor(partial, whole) estimates acc_partial / acc_whole.
  # Do not divide this correlation by mean whole-evaluation reliability.
  ratio_accuracy <- cor_pw

  data.table(
    n = nrow(dt),
    mean_F = f_bar,
    cov_partial_whole = cov_pw,
    cor_partial_whole = cor_pw,
    accuracy_LR = accuracy,
    dispersion = dispersion,
    bias_partial_minus_whole = bias,
    ratio_accuracy_partial_to_whole = ratio_accuracy
  )
}

n_val <- 2038
sigma_a2 <- 500
true_bv <- rnorm(n_val, 0, sqrt(sigma_a2))
lr_dt <- data.table(
  animal_id = sprintf("VAL%05d", seq_len(n_val)),
  birth_year = sample(2018:2021, n_val, replace = TRUE),
  GEBV_whole = true_bv + rnorm(n_val, 0, 8),
  GEBV_partial = 0.90 * true_bv + rnorm(n_val, 0, 10) + 2,
  F_inbreeding = runif(n_val, 0.00, 0.06)
)
fwrite(lr_dt, file.path(outdir, "sim_lr_validation_animals.csv"))
lr_stats <- lr_validation(lr_dt, sigma_a2 = sigma_a2, f_col = "F_inbreeding")
fwrite(lr_stats, file.path(outdir, "LR_validation_statistics.csv"))

# ==========================================================
# Module 8. Save short test summary
# ==========================================================
# Purpose:
#   - Write a compact, R-generated summary of the current run.
#   - The summary records the actual B_agc used so quick tests and manuscript
#     runs cannot be confused.

summary_lines <- c(
  paste0("Selected bulls: ", nrow(selected_bulls)),
  paste0("Reference cows born in 2020 with PCH4 record: ", nrow(ref2020)),
  paste0("Bootstrap replicates used in this R run (B_agc): ", B_agc),
  "",
  "RME standardization check; RMEP/RMEF/RMEG should be close to mean 100 and SD 10 in reference cows:",
  capture.output(print(ref_stats)),
  "",
  "Example AGC with 95% CI for PCH4 vs selected traits:",
  capture.output(print(agc_pch4_adj[trait %in% c("MY", "FY", "FP", "VEF", "VEG"),
                                    .(target, trait, AGC, SE, CI_low, CI_high, B)])),
  "",
  "AGC estimates are not clipped to [-1, 1]; outside-range values are flagged as diagnostics.",
  "RME AGC calculations use simulation-only reliability proxies, not official empirical index reliabilities.",
  "",
  "Example AGC with 95% CI for RME indices vs selected indices:",
  capture.output(print(agc_RME_vs_37[target %in% c("RMEP", "RMEF", "RMEG") & trait %in% c("VEL", "VEF", "VEG"),
                                     .(target, trait, AGC, SE, CI_low, CI_high, B)])),
  "",
  "Simulated daughter-workflow regression (not empirical validation):",
  capture.output(print(reg_RME)),
  "",
  "LR validation statistics (ratio of accuracies = cor(partial, whole)):",
  capture.output(print(lr_stats))
)
writeLines(summary_lines, file.path(outdir, "test_summary.txt"))

cat("Simulation test finished. Results saved to:", outdir, "\n")
