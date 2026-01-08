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

ped_by       <- if (!is.null(arg_list[["--ped_by"]])) as.numeric(arg_list[["--ped_by"]]) else 1
n_grid       <- if (!is.null(arg_list[["--n_grid"]])) as.numeric(arg_list[["--n_grid"]]) else 100
n_bootstraps <- if (!is.null(arg_list[["--n_bootstraps"]])) as.numeric(arg_list[["--n_bootstraps"]]) else 0

# --------------------------------
# Read data
# --------------------------------
data_wide <- read.csv(data_path, stringsAsFactors = FALSE)
data_wide$id <- 1:nrow(data_wide)

# --------------------------------
# NHANES weight handling
# --------------------------------
weight_col <- "WTMEC2YR"

# Copy weight column for later use
data_wide$weight <- data_wide[[weight_col]]

# Drop subjects with NA weight
data_wide <- data_wide %>% filter(!is.na(weight))

# --------------------------------
# Reformat data from wide to long
# Include 'weight' in formula to carry through PED
# --------------------------------
surv_formula <- as.formula(
  paste0("Surv(", duration_col, ", ", event_col, ") ~ ", feature, " + age + weight")
)

data_long <- as_ped(
  data_wide,
  formula = surv_formula,
  cut = seq(0, max(data_wide[[duration_col]], na.rm = TRUE), by = ped_by)
)

# Normalize weight
data_long$w_norm <- data_long$weight / mean(data_long$weight, na.rm = TRUE)

# Save wide-form and long-form data
dir.create(file.path(model_dir, "data"), recursive = TRUE, showWarnings = FALSE)
write.csv(data_long, file.path(model_dir, "data", "data_long.csv"), row.names=FALSE)
write.csv(data_wide, file.path(model_dir, "data", "data_wide.csv"), row.names=FALSE)

# --------------------------------
# Generate prediction grid
# --------------------------------
feature_values <- seq(
  min(data_wide[[feature]], na.rm = TRUE),
  max(data_wide[[feature]], na.rm = TRUE),
  length.out = n_grid
)

tend_values <- seq(
  0,
  max(data_wide[[duration_col]], na.rm = TRUE),
  length.out = n_grid
)

newdata <- expand.grid(tend = tend_values, feature_value = feature_values)
colnames(newdata)[colnames(newdata) == "feature_value"] <- feature
newdata$age <- mean(data_wide$age, na.rm = TRUE)

# --------------------------------
# Output directory
# --------------------------------
output_dir <- file.path(model_dir, "model")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_array <- function(mat, rows, cols, file) {
  rownames(mat) <- round(rows, 3)
  colnames(mat) <- round(cols, 2)
  write.csv(as.data.frame(mat), file, row.names=TRUE, quote=FALSE)
}

# --------------------------------
# Fit PAMM model with weights
# --------------------------------
model_formula <- as.formula(paste("ped_status ~ te(tend, ", feature, ") + age"))

model_pam <- gam(
  formula = model_formula,
  data    = data_long,
  family  = poisson(),
  method  = "REML",
  offset  = data_long$offset,
  weights = data_long$w_norm,
  gamma   = 1.5
)

# --------------------------------
# Compute cumulative hazard over time
# --------------------------------
pred      <- predict(model_pam, newdata = newdata, type = "response", se.fit = TRUE)
hazard    <- matrix(pred$fit, nrow = n_grid, ncol = n_grid, byrow = TRUE)
hazard_SE <- matrix(pred$se.fit, nrow = n_grid, ncol = n_grid, byrow = TRUE)

dt        <- diff(c(0, tend_values))
cumhaz    <- t(apply(hazard, 1, function(x) cumsum(x * dt)))
cumhaz_SE <- t(apply(hazard_SE, 1, function(x) sqrt(cumsum((x * dt)^2))))

hazard_max_follow_up <- cumhaz[, ncol(cumhaz)]
hazard_max_SE        <- cumhaz_SE[, ncol(cumhaz_SE)]

# Weighted reference value for HR = 1
reference_value <- weighted.mean(data_long[[feature]], data_long$w_norm, na.rm = TRUE)
reference_idx <- which.min(abs(feature_values - reference_value))

baseline    <- hazard_max_follow_up[reference_idx]
baseline_SE <- hazard_max_SE[reference_idx]

# Hazard ratio and 95% CI
HR     <- hazard_max_follow_up / baseline
HR_SE  <- hazard_max_SE / baseline

HR_lower <- HR * exp(-1.96 * HR_SE / HR)
HR_upper <- HR * exp(1.96 * HR_SE / HR)

# --------------------------------
# Save outputs
# --------------------------------
save_array(hazard, feature_values, tend_values, file.path(output_dir, "hazard-array.csv"))
save_array(cumhaz, feature_values, tend_values, file.path(output_dir, "cumulative_hazard-array.csv"))

HR_df <- data.frame(
  feature  = feature_values,
  HR       = HR,
  HR_lower = HR_lower,
  HR_upper = HR_upper
)
write.csv(HR_df, file=file.path(output_dir, "hazard_ratio.csv"), row.names=FALSE)

saveRDS(model_pam, file.path(output_dir, "model.rds"))
capture.output(summary(model_pam), file = file.path(output_dir, "model_summary.txt"))
