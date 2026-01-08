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
# Reformat data from wide to long.
# --------------------------------
data_wide <- read.csv(data_path, stringsAsFactors = FALSE)
data_wide$id <- 1:nrow(data_wide)

surv_formula <- as.formula(paste0("Surv(", duration_col, ", ", event_col, ") ~ ", feature, " + age"))

data_long <- as_ped(
  data_wide,
  formula = surv_formula,
  cut = seq(0, max(data_wide[[duration_col]], na.rm = TRUE), by = ped_by)
)

# Save wide-form and long-form data.
dir.create(file.path(model_dir, "data"), recursive = TRUE, showWarnings = FALSE)
write.csv(data_long, file.path(model_dir, "data", "data_long.csv"), row.names=FALSE)
write.csv(data_wide, file.path(model_dir, "data", "data_wide.csv"), row.names=FALSE)

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
  max(data_wide[[duration_col]], na.rm=TRUE),
  length.out=n_grid
)

# Create prediction grid (age set to mean value; has no effect on cumulative hazard ratio).
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
  write.csv(as.data.frame(mat), file, row.names=TRUE, quote=FALSE)
}

# --------------------------------
# Fit model
# --------------------------------

# model_formula <- as.formula(paste("ped_status ~ te(tend, ", feature, ", bs=\"ps\") + age"))
model_formula <- as.formula(paste("ped_status ~ te(tend, ", feature, ") + age"))
model_pam <- gam(
  formula = model_formula,
  data    = data_long,
  family  = poisson(),
  method='REML',
  offset  = data_long$offset,
  gamma=1.5
)


# Compute cumulative hazard over time
pred      = predict(model_pam, newdata=newdata, type ="response", se.fit=TRUE)
hazard    = matrix(pred$fit, nrow=n_grid, ncol=n_grid, byrow=TRUE)
hazard_SE = matrix(pred$se.fit, nrow = n_grid, ncol = n_grid, byrow = TRUE)

dt        = diff(c(0, tend_values))  # time steps
cumhaz    = t(apply(hazard, 1, function(x) cumsum(x * dt)))
cumhaz_SE = t(apply(hazard_SE, 1, function(x) sqrt(cumsum((x * dt)^2))))


# Select max follow-up time.
hazard_max_follow_up = cumhaz[, ncol(cumhaz)]
hazard_max_SE = cumhaz_SE[, ncol(cumhaz_SE)]

# Reference for HR = mean feature value.
reference_value = mean(data_long[[feature]], na.rm = TRUE)
reference_idx = which.min(abs(feature_values - reference_value))

baseline = hazard_max_follow_up[reference_idx]
baseline_SE = hazard_max_SE[reference_idx]

# Hazard ratio and 95% CI.
HR = hazard_max_follow_up / baseline
HR_SE = hazard_max_SE / baseline

HR_lower = HR * exp(-1.96 * HR_SE / HR)
HR_upper = HR * exp(1.96 * HR_SE / HR)

save_array(hazard, feature_values, tend_values, file.path(output_dir, "hazard-array.csv"))
save_array(cumhaz, feature_values, tend_values, file.path(output_dir, "cumulative_hazard-array.csv"))

HR_df = data.frame(feature=feature_values, HR=HR, HR_lower=HR_lower, HR_upper=HR_upper)
write.csv(HR_df, file=file.path(output_dir, "hazard_ratio.csv"), row.names=FALSE)

saveRDS(model_pam, file.path(output_dir, "model.rds"))
capture.output(summary(model_pam), file = file.path(output_dir, "model_summary.txt"))



# # --------------------------------
# # Bootstrap
# # --------------------------------
# if (n_bootstraps > 0) {

#   boot_output_dir <- file.path(output_dir, "bootstrapped_estimates")
#   dir.create(boot_output_dir, recursive = TRUE, showWarnings = FALSE)

#   # Count existing CSVs and start from next index
#   existing <- list.files(boot_output_dir, pattern = "\\.csv$")
#   start_index <- length(existing) + 1

#   message("Bootstrapping...")
#   pb <- txtProgressBar(min = 0, max = n_bootstraps, style = 3)

#   for (i in 1:n_bootstraps) {
#     b <- start_index + i - 1
#     setTxtProgressBar(pb, i)
#     set.seed(1000 + b)

#     # Bootstrap sample
#     sample_ids <- sample(seq_len(nrow(data_wide)), replace = TRUE)
#     data_boot <- data_wide[sample_ids, ]
#     data_boot$id <- seq_len(nrow(data_boot))

#     data_boot_long <- as_ped(
#       data_boot, formula=surv_formula,
#       cut=seq(0, max(data_boot[[duration_col]], na.rm=TRUE), by=ped_by)
#     )

#     if (sum(data_boot[[event_col]]) < 5) next

#     model_b <- gam(
#       model_formula,
#       data = data_boot_long,
#       family = poisson(),
#       offset = data_boot_long$offset
#     )

#     hazard_b <- matrix(
#       predict(model_b, newdata=newdata, type="response"),
#       nrow=n_grid, ncol=n_grid, byrow=TRUE
#     )

#     dt <- diff(c(0, tend_values))
#     cumulative_hazard_b <- t(apply(hazard_b, 1, cumsum) * dt)

#     out_file <- file.path(
#       boot_output_dir,
#       paste0("cumulative_hazard-array__", b, ".csv")
#     )

#     save_array(cumulative_hazard_b, feature_values, tend_values, out_file)
#   }
# }
