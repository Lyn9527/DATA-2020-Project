
setwd("/Users/zhouwenjun/Desktop/DATA 2020 HW/Final Project")  
pkgs <- c("bcf", "SuperLearner", "tidyverse", "doParallel", "data.table")
inst <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(inst)) install.packages(inst)
invisible(lapply(pkgs, library, character.only = TRUE))

# data-generating function
set.seed(42)

generate_data_with_Y <- function(n, homogeneous, linear) {
  p  <- 3
  xm <- matrix(rnorm(n * p), n); colnames(xm) <- paste0("x", 1:3)
  x1 <- xm[, 1]; x2 <- xm[, 2]; x3 <- xm[, 3]
  x4 <- rbinom(n, 1, 0.5)
  x5 <- sample(1:3, n, TRUE)
  
  ## treatment effect
  tau_x <- if (homogeneous) rep(3, n) else 1 + 2 * x2 * x5
  
  ## prognostic function
  g_x5 <- c(`1` = 2, `2` = -1, `3` = -4)[as.character(x5)]
  mu_x <- if (linear) 1 + g_x5 + x1 * x3 else -6 + g_x5 + 6 * abs(x3 - 1)
  
  ## propensity
  pi_x <- 0.8 * pnorm(3 * mu_x / sd(mu_x) - 0.5 * x1) + 0.05 + runif(n) / 10
  pi_x <- pmin(pmax(pi_x, 0), 1)
  
  ## treatment & outcome
  T  <- rbinom(n, 1, pi_x)
  mu <- mu_x + tau_x * T
  sigma <- diff(range(mu_x + tau_x * pi_x)) / 8
  Y  <- mu + sigma * rnorm(n)
  
  data.frame(Y, T, xm, x4, x5 = factor(x5),
             mu_x, tau_x_true = tau_x, pi_x_true = pi_x)
}

# 8 dataset designs (according to paper)
params <- tibble(
  homogeneous = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  linear      = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
  n           = c(250, 250, 250, 250, 500, 500, 500, 500),
  tag         = paste0("set_", 1:8)
)

# produce the data set
data_sets <- params %>%
  pmap(~ generate_data_with_Y(..3, ..1, ..2)) %>%
  setNames(params$tag)
print(head(data_sets$set_1))

# fit BCF to every data set
dir.create("bcf_outputs", showWarnings = FALSE)

cl <- makePSOCKcluster(max(1, parallel::detectCores() - 1))
registerDoParallel(cl)

results <- foreach(ds = data_sets, nm = names(data_sets),
                   .packages = c("bcf", "SuperLearner")) %dopar% {
                     
                     ## extract
                     dat <- ds
                     X   <- model.matrix(~ . -Y -T -mu_x -tau_x_true -pi_x_true, dat)[, -1]
                     y   <- dat$Y
                     z   <- dat$T
                     
                     ## propensity
                     ps <- SuperLearner(z, X, family = binomial(),
                                        SL.library = c("SL.glm","SL.gbm"))$SL.predict
                     
                     ## BCF
                     fit <- bcf(y, z,
                                x_control      = X, 
                                x_moderate     = X,    
                                pihat          = ps,      
                                nburn          = 50,      # Burn-in iterations
                                nsim           = 50,      # MCMC iterations to save
                                use_muscale = TRUE,      
                                use_tauscale = TRUE      
                     )
                     
                     ## summaries
                     cate_hat <- colMeans(fit$tau)
                     ate_hat  <- mean(cate_hat)
                     ci_ate   <- quantile(rowMeans(fit$tau), c(.025,.975))
                     
                     list(tag     = nm,
                          ate     = ate_hat,
                          ci_low  = ci_ate[1],
                          ci_high = ci_ate[2],
                          rmse    = sqrt(mean((cate_hat - dat$tau_x_true)^2)),
                          fit     = fit)       
                   }

stopCluster(cl)

# save metrics & full objects
metrics <- map_dfr(results, ~
                     tibble(tag   = .$tag,
                            ATE   = .$ate,
                            CI_lo = .$ci_low,
                            CI_hi = .$ci_high,
                            RMSE  = .$rmse))

fwrite(metrics, "bcf_outputs/metrics.csv")

for (r in results)
  saveRDS(r$fit, file = file.path("bcf_outputs",
                                  paste0(r$tag, "_bcf_fit.rds")))

cat("\n✓ finished – metrics written to  'bcf_outputs/metrics.csv'\n")
