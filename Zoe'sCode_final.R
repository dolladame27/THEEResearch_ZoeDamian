
#' ---
#' title: "Evolutionary Rescue Under Density-Dependent Regulation"
#' author: "Zoe Schuler"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---


#### MODEL PARAMETERS ####
N0  <- 10000   # Initial total population size (individuals at start of simulation)
p0  <- 0.01    # Initial mutant fraction (proportion of N0 that are mutants)
K   <- 50000   # Carrying capacity (environmental limit that scales density-dependent costs)

w_a <- 0.95    # Wildtype intrinsic fitness: per-capita offspring at low density (< 1 = declining)
w_A <- 1.05    # Mutant intrinsic fitness: per-capita offspring at low density (> 1 = growing)

f_a <- 1.0     # Wildtype density-sensitivity: how strongly fitness declines with crowding (wildtype density sensitivity is purely logistic)
f_A <- 1.25    # Mutant density-sensitivity: f_A > f_a means mutants pay higher density costs

tmax <- 300    # Maximum generations simulated (time horizon for rescue/extinction)


#### POPULATION INITIALIZATION ####
Na  <- numeric(tmax)           # wildtype counts per generation
N_A <- numeric(tmax)           # mutant counts per generation
Na[1]  <- N0 * (1 - p0)        # Initialize wildtypes: (1-p0) fraction of initial population
N_A[1] <- N0 * p0              # Initialize mutants: p0 fraction of initial population


#### DENSITY-DEPENDENT FITNESS FUNCTION ####
fitness <- function(Ntot, w, f) {
  # Logistic fitness function implements density dependence:
  # - w = intrinsic fitness at Ntot ≈ 0
  # - f scales how quickly fitness declines as Ntot approaches K
  # - At Ntot = K, fitness = 1 + (w-1)*(1-f), which can be < 1 for mutants when f_A > 1/w_A
  
  f_raw <- 1 + (w - 1) * (1 - f * Ntot / K)
  pmax(f_raw, 1e-9)            # fitness at 10^-9 to prevent negative/zero values
}


#### SIMULATION DYNAMICS ####
t_last <- 1                    # Track final generation reached

for (t in 2:tmax) {
  # 1. Calculate total population size from previous generation
  Ntot <- Na[t-1] + N_A[t-1]
  
  # 2. Determine genotype-specific effective fitness (density-dependent)
  wa_eff <- fitness(Ntot, w_a, f_a)  # Wildtype fitness this generation
  wA_eff <- fitness(Ntot, w_A, f_A)  # Mutant fitness this generation
  
  # 3. Stochastic reproduction via Poisson branching process
  # Each adult produces Poisson(w_eff) offspring, then dies (no overlap so discrete)
  Na[t]  <- rpois(1, Na[t-1]  * wa_eff)   # Realized wildtype offspring count
  N_A[t] <- rpois(1, N_A[t-1] * wA_eff)   # Realized mutant offspring count
  
  t_last <- t  # Update final generation count
  
  # 4. Termination condition: extinction if total population < 1 individual
  if ((Na[t] + N_A[t]) < 1) break
  
  # Print progress every 10 generations for monitoring
  if (t %% 10 == 0) {
    cat("Gen:", t, "Ntot =", Ntot, 
        "wa_eff =", round(wa_eff, 4), 
        "wA_eff =", round(wA_eff, 4), "\n")
  }
}


#### VISUALIZATION OF POPULATION TRAJECTORIES ####
time <- 1:t_last
plot(time, (Na + N_A)[time], type="l", lwd=2, col="black",
     ylab="Population size", xlab="Generation",
     main="Haploid Rescue with Differential Density Dependence")
lines(time, Na[time],  col="red",  lty=2)   # Dashed red: wildtype dynamics
lines(time, N_A[time], col="blue")          # Solid blue: mutant dynamics
legend("topright", legend=c("Total","Wildtype a","Mutant A"),
       col=c("black","red","blue"), lty=c(1,2,1), bty="n")


#### STAGE 1: RESCUE PROBABILITY vs. INITIAL MUTANT FRACTION ####

# --- EXPERIMENTAL DESIGN PARAMETERS ---
p0_grid <- exp(seq(log(1/N0), log(0.2), length.out = 20))  # 20 log-spaced p₀ values
nrep    <- 200                                            # Replicates per p₀ (increase to 500-1000 for smoother curves)
thresh_rescue <- 0.5                                       # Rescue threshold: N_tot ≥ 50% of K


# --- SINGLE SIMULATION ---
run_one_replicate <- function(p0, tmax = 300, rescue_thresh = 0.5*K) {
  # Initialize population
  Na  <- N0 * (1 - p0)    # Wildtype count
  N_A <- N0 * p0           # Mutant count
  
  for (t in 2:tmax) {
    Ntot <- Na + N_A
    
    # Termination checks
    if (Ntot < 1) return(FALSE)                     # Extinction
    if (Ntot >= rescue_thresh) return(TRUE)         # Rescue achieved
    
    # Density-dependent reproduction
    wa_eff <- fitness(Ntot, w_a, f_a)
    wA_eff <- fitness(Ntot, w_A, f_A)
    
    Na  <- rpois(1, Na  * wa_eff)
    N_A <- rpois(1, N_A * wA_eff)
  }
  
  # Final check at tmax
  return((Na + N_A) >= rescue_thresh)
}

# --- SWEEP ACROSS p₀ ---
rescue_stats <- do.call(rbind, lapply(p0_grid, function(p) {
  # Run nrep independent replicates and calculate rescue probability
  rescue_hits <- replicate(nrep, run_one_replicate(p))
  data.frame(p0 = p, rescue_prob = mean(rescue_hits), n_simulated = nrep)
}))

# --- IDENTIFY 90% THRESHOLD ---
P_target <- 0.9
p0_threshold <- min(rescue_stats$p0[rescue_stats$rescue_prob >= P_target])

# --- VISUALIZATION ---
plot(rescue_stats$p0, rescue_stats$rescue_prob, 
     log = "x", ylim = c(0, 1),
     type = "b", pch = 19, lwd = 2,
     xlab = expression(p[0] ~ "(initial mutant fraction, log scale)"),
     ylab = "Rescue probability",
     main = sprintf("Stage 1: Rescue Probability vs p₀ (nrep = %d)", nrep))
abline(v = p0_threshold, lty = 2, col = "red")
mtext(sprintf("p₀* ≈ %.3g for P ≥ %.2f", p0_threshold, P_target), 
      side = 3, line = 0.5, col = "red")

#### STAGE 2: SCALING OF p₀* WITH N₀ ####

# --- EXPERIMENTAL DESIGN PARAMETERS ---
N0_vec  <- round(exp(seq(log(200), log(min(0.8*K, 40000)), length.out = 8)))
# 8 N₀ values: 200, 500, 1,000, 2,500, 5,000, 10,000, 20,000, 40,000

p0_grid <- exp(seq(log(1 / max(N0_vec)), log(0.2), length.out = 25))
# Finer grid than Stage 1 (25 points) to improve threshold resolution

# --- SIMULATION (replicates Stage 1 mode) ---
estimate_threshold <- function(N0val) {
  sim_once <- function(p0) {
    Na  <- N0val * (1 - p0)
    N_A <- N0val * p0
    
    for (t in 2:tmax) {
      Ntot <- Na + N_A
      if (Ntot < 1) return(FALSE)                      # Extinction
      if (Ntot >= thresh_rescue * K) return(TRUE)      # Rescue threshold
      
      wa <- fitness(Ntot, w_a, f_a)
      wA <- fitness(Ntot, w_A, f_A)
      Na  <- rpois(1, Na  * wa)
      N_A <- rpois(1, N_A * wA)
    }
    (Na + N_A) >= thresh_rescue * K
  }
  
  # Rescue probability across p₀ grid
  probs <- sapply(p0_grid, function(p) mean(replicate(nrep, sim_once(p))))
  
  # Extract p₀*: smallest p₀ with P(rescue) ≥ 0.9
  i <- which(probs >= P_target)
  if (length(i)) min(p0_grid[i]) else NA_real_
}

# SWEEP ACROSS N₀ VALUES
p0_req <- sapply(N0_vec, estimate_threshold)

# --- STATISTICAL ANALYSIS ---
# Filter valid thresholds and compute absolute mutant counts
ok <- is.finite(p0_req) & p0_req > 0
scaling_data <- data.frame(
  N0 = N0_vec[ok],
  p0_req = p0_req[ok],
  N_mutant = N0_vec[ok] * p0_req[ok]
)

# Log-log regression: p₀* ∝ N₀^β
fit <- lm(log(p0_req[ok]) ~ log(N0_vec[ok]))
beta_coef <- coef(fit)[2]  # Scaling exponent β
r_squared <- summary(fit)$r.squared

# --- VISUALIZATION ---
plot(scaling_data$N0, scaling_data$p0_req, 
     log = "xy", type = "b", pch = 19, lwd = 2,
     xlab = expression(N[0] ~ "(initial population, log scale)"),
     ylab = expression(p[0]^"*" ~ "(threshold fraction, log scale)"),
     main = sprintf("Stage 2: Scaling β = %.3f (R² = %.3f)", beta_coef, r_squared))

# Add regression line
abline(fit, col = "red", lty = 2)

# Mark points where threshold wasn't reached
if (any(!ok)) {
  points(N0_vec[!ok], rep(min(p0_grid, na.rm = TRUE), sum(!ok)), pch = 1)
  mtext("open circles: p₀* not reached on grid", side = 3, line = 0.5, cex = 0.8)
}

# --- HEAD-COUNT RULE VERIFICATION ---
cat("Scaling exponent (β):", beta_coef, "\n")
print(scaling_data)  # Shows N_mutant ≈ constant across N₀ ranges

#### STAGE 2B:SIMPLE LOGISTIC, WITHOUT F ####
fitness <- function(Ntot, w) {
  # Simple logistic: per-capita fitness declines linearly with Ntot
  f_raw <- 1 + (w - 1) * (1 - Ntot / K)
  pmax(f_raw, 1e-9)
}
# --- EXPERIMENTAL DESIGN PARAMETERS ---
N0_vec  <- round(exp(seq(log(200), log(min(0.8*K, 40000)), length.out = 8)))
# 8 N₀ values: 200, 500, 1,000, 2,500, 5,000, 10,000, 20,000, 40,000

p0_grid <- exp(seq(log(1 / max(N0_vec)), log(0.2), length.out = 25))
# Finer grid than Stage 1 (25 points) to improve threshold resolution

# --- SIMULATION (replicates Stage 1 mode) ---
estimate_threshold <- function(N0val) {
  sim_once <- function(p0) {
    Na  <- N0val * (1 - p0)
    N_A <- N0val * p0
    
    for (t in 2:tmax) {
      Ntot <- Na + N_A
      if (Ntot < 1) return(FALSE)                      # Extinction
      if (Ntot >= thresh_rescue * K) return(TRUE)      # Rescue threshold
      
      wa <- fitness(Ntot, w_a)
      wA <- fitness(Ntot, w_A)
      Na  <- rpois(1, Na  * wa)
      N_A <- rpois(1, N_A * wA)
    }
    (Na + N_A) >= thresh_rescue * K
  }
  
  # Rescue probability across p₀ grid
  probs <- sapply(p0_grid, function(p) mean(replicate(nrep, sim_once(p))))
  
  # Extract p₀*: smallest p₀ with P(rescue) ≥ 0.9
  i <- which(probs >= P_target)
  if (length(i)) min(p0_grid[i]) else NA_real_
}

# SWEEP ACROSS N₀ VALUES
p0_req <- sapply(N0_vec, estimate_threshold)

# --- STATISTICAL ANALYSIS ---
# Filter valid thresholds and compute absolute mutant counts
ok <- is.finite(p0_req) & p0_req > 0
scaling_data <- data.frame(
  N0 = N0_vec[ok],
  p0_req = p0_req[ok],
  N_mutant = N0_vec[ok] * p0_req[ok]
)

# Log-log regression: p₀* ∝ N₀^β
fit <- lm(log(p0_req[ok]) ~ log(N0_vec[ok]))
beta_coef <- coef(fit)[2]  # Scaling exponent β
r_squared <- summary(fit)$r.squared

# --- VISUALIZATION ---
plot(scaling_data$N0, scaling_data$p0_req, 
     log = "xy", type = "b", pch = 19, lwd = 2,
     xlab = expression(N[0] ~ "(initial population, log scale)"),
     ylab = expression(p[0]^"*" ~ "(threshold fraction, log scale)"),
     main = sprintf("Stage 2: Scaling β = %.3f (R² = %.3f)", beta_coef, r_squared))

# Add regression line
abline(fit, col = "red", lty = 2)

# Mark points where threshold wasn't reached
if (any(!ok)) {
  points(N0_vec[!ok], rep(min(p0_grid, na.rm = TRUE), sum(!ok)), pch = 1)
  mtext("open circles: p₀* not reached on grid", side = 3, line = 0.5, cex = 0.8)
}

# --- HEAD-COUNT RULE VERIFICATION ---
cat("Scaling exponent (β):", beta_coef, "\n")
print(scaling_data)  # Shows N_mutant ≈ constant across N₀ ranges
#### STAGE 3: EFFECT OF f_A ON RESCUE AT DIFFERENT p₀ ####

# --- EXPERIMENTAL DESIGN PARAMETERS ---
p0_test_values <- c(0.005, 0.01, 0.02)  # Below, at, and above critical threshold
f_A_values <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)  #  density sensitivity scan
nrep_stage3 <- 200  # Replicates per (p₀, f_A) combination

# --- SIMULATION (vectorized for eff ---
run_rescue_trial <- function(p0_fixed, f_A_value) {
  Na  <- N0 * (1 - p0_fixed)
  N_A <- N0 * p0_fixed
  
  for (t in 2:tmax) {
    Ntot <- Na + N_A
    
    if (Ntot < 1) return(FALSE)
    if (Ntot >= thresh_rescue * K) return(TRUE)
    
    wa_eff <- fitness(Ntot, w_a, f_a)
    wA_eff <- fitness(Ntot, w_A, f_A_value)  # Vary crowding sensitivity
    
    Na  <- rpois(1, Na  * wa_eff)
    N_A <- rpois(1, N_A * wA_eff)
  }
  
  return((Na + N_A) >= thresh_rescue * K)
}

# --- SWEEP ---
# Results matrix: rows = f_A, cols = p₀
crowding_results <- matrix(0, nrow = length(f_A_values), ncol = length(p0_test_values))
dimnames(crowding_results) <- list(
  f_A = f_A_values,
  p0 = paste0("p0=", p0_test_values)
)

for (i in seq_along(f_A_values)) {
  for (j in seq_along(p0_test_values)) {
    cat(sprintf("Running f_A=%.2f, p0=%.3f (%d reps)...\n", 
                f_A_values[i], p0_test_values[j], nrep_stage3))
    
    rescue_hits <- replicate(nrep_stage3, 
                             run_rescue_trial(p0_test_values[j], f_A_values[i]))
    crowding_results[i, j] <- mean(rescue_hits)
  }
}

# --- VISUALIZATION (all curves together) ---
matplot(f_A_values, crowding_results, type = "b", pch = c(19, 21, 24), lty = 1:3,
        xlab = expression(f[A] ~ "(mutant crowding sensitivity)"),
        ylab = "Rescue probability",
        main = "Stage 3: Effect of f_A at Different p₀ Values",
        col = c("black", "blue", "red"))
legend("topright", legend = paste0("p₀ = ", p0_test_values),
       col = c("black", "blue", "red"), pch = c(19, 21, 24), lty = 1:3,
       title = "Initial mutant fraction")
abline(v = 1.25, col = "gray50", lty = 2)  # Mark default f_A

# --- STATISTICAL ANALYSIS ---
# Effect size of f_A at each p₀ (slope of rescue prob vs f_A)
effect_sizes <- apply(crowding_results, 2, function(col) {
  # Linear regression of rescue prob on f_A (excluding extremes)
  fit <- lm(col ~ f_A_values)
  coef(fit)[2]  # Slope: how much rescue drops per unit increase in f_A
})
names(effect_sizes) <- paste0("p0=", p0_test_values)

cat("\nEffect size of f_A (slope):\n")
print(effect_sizes)


#### STAGE 4: REBOUND METRICS ####

# --- EXPERIMENTAL DESIGN PARAMETERS ---
p0_subset <- c(0.001, 0.002, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1)
nrep_stage4 <- 200

# --- SIMULATION ENGINE (returns metrics, not full trajectories) ---
run_rebound_metrics <- function(p0, N0 = 10000, f_A = 2.0) {
  Na <- N0 * (1 - p0); N_A <- N0 * p0
  min_pop <- Na + N_A
  rebound_t <- NA
  
  for (t in 2:tmax) {
    Ntot <- Na + N_A
    if (Ntot < 1) return(list(time = NA, min = min_pop, outcome = "extinct"))
    if (Ntot >= thresh_rescue * K) {
      rebound_t <- t
      break
    }
    # Update with density dependence
    Na <- rpois(1, Na * fitness(Ntot, w_a, f_a))
    N_A <- rpois(1, N_A * fitness(Ntot, w_A, f_A))
    min_pop <- min(min_pop, Na + N_A)
  }
  
  list(time = rebound_t, min = min_pop, 
       outcome = if (!is.na(rebound_t)) "rescued" else "extinct")
}

# --- BATCH SWEEP (store only metrics) ---
rebound_summary <- do.call(rbind, lapply(p0_subset, function(p) {
  results <- replicate(nrep_stage4, run_rebound_metrics(p), simplify = FALSE)
  data.frame(
    p0 = p,
    rescue_prob = mean(sapply(results, `[[`, "outcome") == "rescued"),
    median_rebound = median(unlist(sapply(results, `[[`, "time")), na.rm = TRUE),
    median_min_pop = median(unlist(sapply(results, `[[`, "min"))),
    n_replicates = nrep_stage4
  )
}))

# --- PLOTTING ---
par(mfrow = c(1, 2))

# Plot Rebound time (rescued only)
valid <- !is.na(rebound_summary$median_rebound)
if (any(valid)) {
  plot(rebound_summary$p0[valid], rebound_summary$median_rebound[valid],
       log = "x", type = "b", pch = 19, lwd = 1.5,
       main = "Rebound Time", xlab = "p₀", ylab = "Generations")
}

# Plot: Minimum population (all data)
plot(rebound_summary$p0, rebound_summary$median_min_pop, 
     log = "xy", type = "b", pch = 19, lwd = 1.5,
     main = "Minimum Size", xlab = "p₀", ylab = "Individuals")


#  Population trajectories showing rebound patterns
p_trajectories <- ggplot() +
  lapply(names(example_results), function(p0_str) {
    result <- example_results[[p0_str]]
    p0_val <- as.numeric(p0_str)
    geom_line(data = data.frame(
      time = result$time,
      population = result$total_population,
      p0 = sprintf("p₀ = %.3f", p0_val)
    ), aes(x = time, y = population, color = p0), size = 1.2, alpha = 0.8)
  }) +
  geom_hline(yintercept = 25000, linetype = "dashed", alpha = 0.7, color = "red") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Population Dynamics: Single Population View",
    subtitle = "Total population trajectories showing rescue patterns",
    x = "Generation",
    y = "Total Population Size",
    color = "Initial Mutant Fraction",
    caption = "Red dashed line shows rescue threshold (50% of K = 50,000)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

print(p_trajectories)

