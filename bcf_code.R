# Set working directory (adjust path as needed)
#setwd("/Users/zhouwenjun/Desktop/DATA 2020 HW/Final Project") 

## ── 1. Packages ───────────────────────────────────────────────────────────
# Ensure packages are installed and loaded
pkgs <- c("bcf", "tidyverse", "doParallel", "data.table", "stats")
inst <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(inst)) install.packages(inst)
invisible(lapply(pkgs, library, character.only = TRUE))

## ── 2. Data-generating function 
# This function will be called inside the simulation loop.
# NOTE: Includes hardcoded tau_x = 3 overwrite. True effect is always 3.
generate_data_with_Y <- function(n, homogeneous, linear) {
  p = 3; x_mat <- matrix(rnorm(n*p), nrow=n); x1 <- x_mat[,1]; x2 <- x_mat[,2]; x3 <- x_mat[,3]
  x4 <- rbinom(n, 1, 0.5); x5 <- sample(1:3, n, replace = TRUE)
  if (homogeneous) { tau_x_intended <- rep(3, n) } else { tau_x_intended <- 1 + 2 * x2 * x4 } # Calculated but overwritten
  if (linear) { g_x5 <- sapply(x5, function(x) if (x == 1) 2 else if (x == 2) -1 else -4); mu_x <- 1 + g_x5 + x1 * x3
  } else { g_x5 <- sapply(x5, function(x) if (x == 1) 2 else if (x == 2) -1 else -4); mu_x <- -6 + g_x5 + 6 * abs(x3 - 1) }
  # --- Add robustness checks for sd and sigma ---
  sd_mu_x <- sd(mu_x); if (is.na(sd_mu_x) || sd_mu_x < 1e-6) sd_mu_x <- 1
  pi_x <- 0.8 * pnorm(3 * mu_x / sd_mu_x - 0.5 * x1) + 0.05 + runif(n) / 10
  pi_x <- pmin(pmax(pi_x, 0), 1); T <- rbinom(n, 1, pi_x)
  tau_x <- rep(3, n); mu <- mu_x + tau_x * T # OVERWRITE tau_x=3 used here
  sigma_range <- range(mu_x + tau_x * pi_x); sigma <- diff(sigma_range) / 8; if (is.na(sigma) || sigma <= 1e-6) sigma <- 1
  Y <- mu + sigma * rnorm(n)
  # --- End robustness checks ---
  data <- data.frame(x1, x2, x3, x4, x5, mu_x, tau_x, pi_x, T, Y)
  data <- cbind(data, x_mat); return(data)
}


## ── 3. Define design settings and SUBSET to n=250 ───────────────────────
# Define all 8 settings using the provided list format
all_params <- list(
  list(homogeneous = TRUE, linear = TRUE, n = 250),  # 1 (n=250)
  list(homogeneous = TRUE, linear = FALSE, n = 250), # 2 (n=250)
  list(homogeneous = FALSE, linear = TRUE, n = 250), # 3 (n=250)
  list(homogeneous = FALSE, linear = FALSE, n = 250),# 4 (n=250)
  list(homogeneous = TRUE, linear = TRUE, n = 500),  # 5 (n=500)
  list(homogeneous = TRUE, linear = FALSE, n = 500), # 6 (n=500)
  list(homogeneous = FALSE, linear = TRUE, n = 500), # 7 (n=500)
  list(homogeneous = FALSE, linear = FALSE, n = 500) # 8 (n=500)
)
params_indices_to_run <- 1:4
params_subset <- all_params[params_indices_to_run]
param_names <- paste0("set_", params_indices_to_run)

## ── 4. Simulation Setup ──────────────────────────────────────────────────
N_REPETITIONS <- 50 # Number of repetitions per setting
all_run_metrics <- list() # List to store metrics from each run
dir.create("bcf_outputs", showWarnings = FALSE) # Create output directory for results
num_cores_to_use <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(num_cores_to_use)
registerDoParallel(cl)
cat(paste("Registered", getDoParWorkers(), "parallel workers.\n"))
analysis_base_seed <- 123 # Base seed for repetitions

## ── 5. Main Simulation Loop (Generate data within loop) ───────────────
for (rep in 1:N_REPETITIONS) { # Loop through the 50 repetitions
  
  # <<< Set a UNIQUE seed for THIS repetition's data generation & ANALYSIS >>>
  set.seed(analysis_base_seed + rep)
  
  cat(paste("Starting simulation repetition:", rep, "/", N_REPETITIONS, "\n"))
  
  # ---> Generate the 4 datasets (n=250 only) for THIS repetition <---
  current_data_sets <- vector("list", length(params_subset)) # length is 4
  names(current_data_sets) <- param_names
  generation_successful <- TRUE
  for (i in seq_along(params_subset)) {
    current_param_set <- params_subset[[i]]
    # *** Calling the data generation function directly ***
    current_data_sets[[i]] <- tryCatch({
      generate_data_with_Y(
        n = current_param_set$n,
        homogeneous = current_param_set$homogeneous,
        linear = current_param_set$linear
      )
    }, error = function(e) {
      warning("Data generation failed for setting ", param_names[i], " rep ", rep, ": ", e$message, call. = FALSE)
      generation_successful <- FALSE; return(NULL) })
    if (!generation_successful) break
  }
  # Skip analysis for this rep if data generation failed
  if (!generation_successful) { cat("Skipping analysis for repetition", rep, "\n"); next }
  
  # ---> Fit BCF to the newly generated data sets <---
  run_results <- foreach(
    ds = current_data_sets, # Use the generated datasets
    nm = names(current_data_sets),
    .packages = c("bcf", "SuperLearner", "dplyr", "tibble", "stats"),
    .errorhandling = "pass"
  ) %dopar% {
    
    # Check if ds is NULL (from potential generation error above)
    if (is.null(ds)) return(list(tag = nm, rep = rep, error = "Input data was NULL"))
    dat <- ds
    
    # --- Data Preparation ---
    prep_result <- tryCatch({
      if (!"x5" %in% names(dat)) stop("Column x5 not found")
      if (!is.factor(dat$x5)) dat$x5 <- factor(dat$x5) # Ensure factor
      formula_X <- ~ x1 + x2 + x3 + x4 + x5
      required_cols_X <- c("x1", "x2", "x3", "x4", "x5"); required_cols_other <- c("Y", "T", "tau_x")
      if (!all(c(required_cols_X, required_cols_other) %in% names(dat))) stop("Missing required columns")
      X  <- model.matrix(formula_X, data = dat)[, -1, drop = FALSE]
      list(X = X, y = dat$Y, z = dat$T, true_tau = dat$tau_x)
    }, error = function(e) { list(error = paste("Data Prep failed:", e$message)) })
    if ("error" %in% names(prep_result)) return(c(list(tag = nm, rep = rep), prep_result))
    X <- prep_result$X; y <- prep_result$y; z <- prep_result$z;
    true_cate_vec <- prep_result$true_tau # Vector of 3s
    true_ate_val  <- 3.0
    
    
    ps <- glm(z ~ ., family = binomial(), data = as.data.frame(X))$fitted.values
    ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)  
    
    # --- BCF Fit ---
    bcf_result <- tryCatch({
      fit <- bcf(y = y, z = z, x_control = X, x_moderate = X, pihat = ps,
                 nburn = 1000, nsim = 1000,
                 use_muscale = TRUE, use_tauscale = TRUE)
      list(fit = fit)
    }, error = function(e) { list(error = paste("BCF failed:", e$message)) })
    if ("error" %in% names(bcf_result)) return(c(list(tag = nm, rep = rep), bcf_result))
    fit <- bcf_result$fit
    
    # --- Calculate Metrics ---
    metrics_result <- tryCatch({
      ate_draws <- rowMeans(fit$tau); ate_hat <- mean(ate_draws)
      ate_ci <- stats::quantile(ate_draws, c(.025, .975), na.rm = TRUE)
      ate_len <- ate_ci[2] - ate_ci[1]; ate_coverage <- (ate_ci[1] <= true_ate_val & ate_ci[2] >= true_ate_val)
      ate_rmse <- sqrt((ate_hat - true_ate_val)^2)
      cate_hat <- colMeans(fit$tau); cate_rmse <- sqrt(mean((cate_hat - true_cate_vec)^2))
      cate_cis <- apply(fit$tau, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      if(!is.matrix(cate_cis) || nrow(cate_cis) != 2) stop("CATE CI calculation failed")
      cate_len_vec <- cate_cis[2, ] - cate_cis[1, ]; cate_cover_vec <- (cate_cis[1, ] <= true_cate_vec & cate_cis[2, ] >= true_cate_vec)
      mean_cate_len <- mean(cate_len_vec, na.rm = TRUE); mean_cate_coverage <- mean(cate_cover_vec, na.rm = TRUE)
      list(ate_hat=ate_hat, ate_coverage = ate_coverage, ate_len = ate_len, ate_rmse = ate_rmse,
           cate_rmse = cate_rmse, cate_coverage = mean_cate_coverage, cate_len = mean_cate_len)
    }, error = function(e) { list(error = paste("Metric calculation failed:", e$message)) })
    if ("error" %in% names(metrics_result)) return(c(list(tag = nm, rep = rep), metrics_result))
    
    # Combine results successfully
    return(c(list(tag = nm, sim_run = rep), metrics_result))
    
  } # End of foreach loop
  
  # Check for errors returned from foreach
  errors_this_rep <- Filter(function(x) !is.null(x) && "error" %in% names(x), run_results)
  if (length(errors_this_rep) > 0) {
    cat("Errors encountered in repetition", rep, ":\n")
    for(err_res in errors_this_rep){ cat("  Setting:", err_res$tag, "- Rep:", err_res$sim_run, "- Error:", err_res$error, "\n") }
  }
  
  # Store metrics from this repetition
  required_metric_names <- c("tag", "sim_run", "ate_coverage", "ate_len", "ate_rmse", "cate_rmse", "cate_coverage", "cate_len")
  valid_results <- Filter(function(x) !is.null(x) && !("error" %in% names(x)) && all(required_metric_names %in% names(x)), run_results)
  if (length(valid_results) > 0) {
    current_rep_metrics <- map_dfr(valid_results, ~ as_tibble(.x))
    all_run_metrics[[rep]] <- current_rep_metrics
  } else {
    if(length(errors_this_rep) == 0) cat("Warning: No valid results obtained for repetition", rep, "and no specific errors caught.\n")
  }
  
} # End of main simulation loop (rep)

# Stop the parallel cluster
stopCluster(cl)
registerDoSEQ()
cat("\n✓ Simulation loop finished.\n")

## ── 6. Aggregate, Summarize, and Save Metrics (for n=250 settings) ─────
if (length(all_run_metrics) > 0) {
  final_metrics <- bind_rows(all_run_metrics)
  
  # Convert the SUBSETTED params list to a data frame/tibble for joining
  # Need to handle the list structure correctly
  params_list_for_df <- all_params[params_indices_to_run]
  # Add the tag names before converting
  names(params_list_for_df) <- param_names
  params_df <- map_dfr(params_list_for_df, ~as_tibble(.x), .id = "tag")
  
  
  summary_metrics <- final_metrics %>% group_by(tag) %>%
    summarise(n_successful_runs = n(),
              ATE_RMSE_mean    = mean(ate_rmse, na.rm = TRUE), ATE_Coverage_mean= mean(ate_coverage, na.rm = TRUE), ATE_Len_mean = mean(ate_len, na.rm = TRUE),
              CATE_RMSE_mean   = mean(cate_rmse, na.rm = TRUE), CATE_Coverage_mean = mean(cate_coverage, na.rm = TRUE), CATE_Len_mean = mean(cate_len, na.rm = TRUE),
              .groups = 'drop' ) %>%
    left_join(params_df, by = "tag") %>% # Join with params info
    select(tag, homogeneous, linear, n, n_successful_runs,
           ATE_RMSE_mean, ATE_Coverage_mean, ATE_Len_mean,
           CATE_RMSE_mean, CATE_Coverage_mean, CATE_Len_mean) # Reorder cols
  
  # --- Create Table similar to Table 2 ---
  table_final <- summary_metrics %>%
    mutate(Effect_Type = ifelse(homogeneous, "Homogeneous", "Heterogeneous")) %>%
    select( Setting = tag, N = n, Effect_Type,
            ATE_rmse = ATE_RMSE_mean, ATE_cover = ATE_Coverage_mean, ATE_len = ATE_Len_mean,
            CATE_rmse = CATE_RMSE_mean, CATE_cover = CATE_Coverage_mean, CATE_len = CATE_Len_mean ) %>%
    arrange(Effect_Type, Setting) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  
  print("Summary Metrics Table (similar structure to Table 2):"); print(as.data.frame(table_final))
  
  # Save the raw metrics from all runs
  raw_output_filename <- file.path("bcf_outputs", paste0("metrics_n250_all_", N_REPETITIONS, "_reps_DGPinloop.csv"))
  fwrite(final_metrics, raw_output_filename)
  cat(paste("\n✓ Raw metrics (n=250, DGPinloop) from all replicates saved to:", raw_output_filename, "\n"))
  
  # Save the summary metrics (means)
  summary_output_filename <- file.path("bcf_outputs", paste0("metrics_n250_summary_T2format_", N_REPETITIONS, "_reps_DGPinloop.csv"))
  fwrite(summary_metrics, summary_output_filename) # Save the unformatted summary
  cat(paste("✓ Summary metrics (means) saved to:", summary_output_filename, "\n"))
} else { cat("No metrics were successfully collected across the simulation replicates.\n") }
cat("\n✓ Script finished.\n")