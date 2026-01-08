#!/usr/bin/env Rscript

# -----------------------------
# Parse arguments
# -----------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(pammtools)
  library(mgcv)
  library(survival)
})

parse_args <- function(args) {
  arg_list <- list()
  for (arg in args) {
    s <- strsplit(arg, "=")[[1]]
    if (length(s) == 2) arg_list[[s[1]]] <- s[2]
  }
  arg_list
}

args <- commandArgs(trailingOnly = TRUE)
arg_list <- parse_args(args)

data_path    <- arg_list[["--data_path"]]
duration_col <- arg_list[["--duration_col"]]
event_col    <- arg_list[["--event_col"]]
feature      <- arg_list[["--feature"]]
model_dir    <- arg_list[["--model_dir"]]

ped_by <- if (!is.null(arg_list[["--ped_by"]])) {
  as.numeric(arg_list[["--ped_by"]])
} else {
  1
}

n_grid <- if (!is.null(arg_list[["--n_grid"]])) {
  as.numeric(arg_list[["--n_grid"]])
} else {
  100
}

# --------------------------------
# Reformat data from wide to long.
# --------------------------------
data_wide <- read.csv(data_path, stringsAsFactors = FALSE)
data_wide$id <- seq_len(nrow(data_wide))

surv_formula <- as.formula(
  paste0("Surv(", duration_col, ", ", event_col, ") ~ ", feature, " + age")
)

data_long <- as_ped(
  data_wide,
  formula = surv_formula,
  cut = seq(0, max(data_wide[[duration_col]], na.rm = TRUE), by = ped_by)
)

# Save wide-form and long-form data.
dir.create(file.path(model_dir, "data"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  data_long,
  file.path(model_dir, "data", "data_long.csv"),
  row.names = FALSE
)
write.csv(
  data_wide,
  file.path(model_dir, "data", "data_wide.csv"),
  row.names = FALSE
)

# --------------------------------
# Generate a prediction grid.
# --------------------------------

# Vary feature values.
feature_values <- seq(
  min(data_wide[[feature]], na.rm = TRUE),
  max(data_wide[[feature]], na.rm = TRUE),
  length.out = n_grid
)

# Vary time-horizon values.
tend_values <- seq(
  0,
  max(data_wide[[duration_col]], na.rm = TRUE),
  length.out = n_grid
)

# Create prediction grid (age set to mean value).
newdata <- expand.grid(tend = tend_values, feature_value = feature_values)
colnames(newdata)[colnames(newdata) == "feature_value"] <- feature
newdata$age <- mean(data_wide$age, na.rm = TRUE)

# --------------------------------
# Output dir
# --------------------------------
output_dir <- file.path(model_dir, "model")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_array <- function(mat, rows, cols, file) {
  rownames(mat) <- round(rows, 3)
  colnames(mat) <- round(cols, 2)
  write.csv(as.data.frame(mat), file, row.names = TRUE, quote = FALSE)
}

# --------------------------------
# Fit model
# --------------------------------

model_formula <- as.formula(paste("ped_status ~ te(tend, ", feature, ") + age"))
model_pam <- gam(
  formula = model_formula,
  data    = data_long,
  family  = poisson(),
  method  = "REML",
  offset  = data_long$offset,
  gamma   = 1.5
)

# Compute cumulative hazard over time
pred <- predict(
  model_pam,
  newdata = newdata,
  type = "response",
  se.fit = TRUE
)

hazard <- matrix(pred$fit, nrow = n_grid, ncol = n_grid, byrow = TRUE)
hazard_se <- matrix(pred$se.fit, nrow = n_grid, ncol = n_grid, byrow = TRUE)

dt        <- diff(c(0, tend_values))  # time steps
cumhaz    <- t(apply(hazard, 1, function(x) cumsum(x * dt)))
cumhaz_se <- t(apply(hazard_se, 1, function(x) sqrt(cumsum((x * dt)^2))))


# select max follow-up time.
hazard_max_follow_up <- cumhaz[, ncol(cumhaz)]
hazard_max_se <- cumhaz_se[, ncol(cumhaz_se)]

# Reference for hr = mean feature value.
reference_value <- mean(data_long[[feature]], na.rm = TRUE)
reference_idx <- which.min(abs(feature_values - reference_value))

baseline <- hazard_max_follow_up[reference_idx]
baseline_se <- hazard_max_se[reference_idx]

# Hazard ratio and 95% CI.
hr <- hazard_max_follow_up / baseline
hr_se <- hazard_max_se / baseline

hr_lower <- hr * exp(-1.96 * hr_se / hr)
hr_upper <- hr * exp(1.96 * hr_se / hr)

save_array(
  hazard,
  feature_values,
  tend_values,
  file.path(output_dir, "hazard-array.csv")
)

save_array(
  cumhaz,
  feature_values,
  tend_values,
  file.path(output_dir, "cumulative_hazard-array.csv")
)

hr_df <- data.frame(
  feature = feature_values,
  hr = hr,
  hr_lower = hr_lower,
  hr_upper = hr_upper
)

write.csv(
  hr_df,
  file = file.path(output_dir, "hazard_ratio.csv"),
  row.names = FALSE
)

saveRDS(model_pam, file.path(output_dir, "model.rds"))
capture.output(
  summary(model_pam),
  file = file.path(output_dir, "model_summary.txt")
)
