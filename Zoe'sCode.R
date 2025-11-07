#Question: What is the ideal fraction of mutants (p0) to introduce
# to rescue a declining population?
## PARAMETERS
N0  <- 10000   # initial total population size
p0  <- 0.002      # initial fraction of mutants
K   <- 50000   # carrying capacity

w_a <- 0.95      # intrinsic fitness of wildtype (below 1)
w_A <- 1.05      # intrinsic fitness of mutant (above 1)

f_a <- 1.0       # density-sensitivity of wildtype
f_A <- 1.25       # density-sensitivity of mutant

tmax <- 300      # number of generations


## INITIALIZATION
Na  <- numeric(tmax)           # wildtype counts
N_A <- numeric(tmax)           # mutant counts
Na[1]  <- N0 * (1 - p0)        # initialize wildtypes




N_A[1] <- N0 * p0              # initialize mutants


## LOGISTIC FITNESS FUNCTION (for both genotypes)
fitness <- function(Ntot, w, f) {
  f_raw <- 1 + (w - 1) * (1 - f * Ntot / K)  # logistic density dependence with f
  pmax(f_raw, 1e-9)                          # keep fitness positive
}


## SIMULATION
#set.seed(1)
t_last <- 1

for (t in 2:tmax) {
  Ntot <- Na[t-1] + N_A[t-1]   # total population size in previous generation
  
  # effective fitness for each genotype
  wa_eff <- fitness(Ntot, w_a, f_a)
  wA_eff <- fitness(Ntot, w_A, f_A)
  
  # stochastic reproduction
  Na[t]  <- rpois(1, Na[t-1]  * wa_eff)
  N_A[t] <- rpois(1, N_A[t-1] * wA_eff)
  
  t_last <- t
  
  # stop if population extinct
  if ((Na[t] + N_A[t]) < 1) break
  
  # print every 10 generations
  if (t %% 10 == 0) {
    cat("Gen:", t,
        "Ntot =", Ntot,
        "wa_eff =", round(wa_eff, 4),
        "wA_eff =", round(wA_eff, 4), "\n")
  }
}


## PLOT RESULTS
time <- 1:t_last
plot(time, (Na + N_A)[time], type="l", lwd=2, col="black",
     ylab="Population size", xlab="Generation",
     main="Haploid Rescue with Differential Density Dependence")
lines(time, Na[time],  col="red",  lty=2)
lines(time, N_A[time], col="blue")
legend("topright", legend=c("Total","Wildtype a","Mutant A"),
       col=c("black","red","blue"), lty=c(1,2,1), bty="n")

## rescue probability vs p0

# parameters
p0_grid <- exp(seq(log(1/N0), log(0.2), length.out = 20))  #log-spaced grid, for fine resolution at small p0
nrep    <- 200                  # replicates per p0
thresh_rescue <- 0.5             # "rescued if N >= thresh*K"
#set.seed(40)

# single replicate
.sim_once <- function(p0) {
  Na_i <- N0 * (1 - p0)
  NA_i <- N0 * p0
  for (t in 2:tmax) {
    Ntot <- Na_i + NA_i
    if (Ntot < 1) return(FALSE)                       # extinct
    if (Ntot >= thresh_rescue * K) return(TRUE)        # rescued
    
    wa_eff <- fitness(Ntot, w_a, f_a)
    wA_eff <- fitness(Ntot, w_A, f_A)
    
    Na_i <- rpois(1, Na_i * wa_eff)
    NA_i <- rpois(1, NA_i * wA_eff)
  }
  (Na_i + NA_i) >= thresh_rescue * K
}

# run sweep
res_stats <- do.call(rbind, lapply(p0_grid, function(p) {
  hits <- replicate(nrep, .sim_once(p))
  data.frame(p0 = p, rescue_prob = mean(hits))
}))
print(res_stats)

#plot
op <- par(no.readonly=TRUE); on.exit(par(op))
par(mar=c(4,4,1,1))
plot(res_stats$p0, res_stats$rescue_prob, log="x", ylim=c(0,1),
     type="b", pch=19, lwd=2,
     xlab=expression(p[0]~"(initial mutant fraction, log scale)"),
     ylab="Rescue probability")
P_target <- 0.9
idx <- which(res_stats$rescue_prob >= P_target)
if (length(idx) > 0) {
  p0_hit <- min(res_stats$p0[idx])
  abline(v = p0_hit, lty=2)
  mtext(sprintf("p0 ≈ %.3g reaches P≥%.2f", p0_hit, P_target),
        side=3, line=0.2)
}


## --- Minimal: N0 → p0* for 90% rescue  (fixed) -----------------------------
P_target <- 0.9

# set.seed(123)  # optional reproducibility

# N0 values to test
N0_vec  <- round(exp(seq(log(200), log(min(0.8*K, 40000)), length.out = 8)))

# p0 grid (log-spaced)
p0_grid <- exp(seq(log(1 / max(N0_vec)), log(0.2), length.out = 25))

# For each N0, find smallest p0 hitting the target
p0_req <- sapply(N0_vec, function(N0val) {
  sim_once <- function(p0) {
    Na  <- N0val * (1 - p0)
    N_A <- N0val * p0           
    for (t in 2:tmax) {
      Ntot <- Na + N_A
      if (Ntot < 1) return(FALSE)
      if (Ntot >= thresh_rescue * K) return(TRUE)
      wa <- fitness(Ntot, w_a, f_a)
      wA <- fitness(Ntot, w_A, f_A)
      Na  <- rpois(1, Na  * wa)
      N_A <- rpois(1, N_A * wA)
    }
    (Na + N_A) >= thresh_rescue * K
  }
  probs <- sapply(p0_grid, function(p) mean(replicate(nrep, sim_once(p))))
  i <- which(probs >= P_target)
  if (length(i)) min(p0_grid[i]) else NA_real_
})

# Plot (log–log). Filter out missing to avoid log() issues.
ok <- !is.na(p0_req) & p0_req > 0
plot(N0_vec[ok], p0_req[ok], log = "xy", type = "b", pch = 19, lwd = 2,
     xlab = expression(N[0]~"(initial population)"),
     ylab = expression(p[0]^"*"~"to reach P[rescue] >= 0.9"))

# Optionally show where target wasn't reached on the grid
if (any(!ok)) {
  points(N0_vec[!ok], rep(min(p0_grid), sum(!ok)), pch = 1)
  mtext("open circles: target not reached on grid", side = 3, line = 0.2, cex = 0.9)
}



cbind(N0 = N0_vec, p0_req = p0_req, N0_times_p0 = N0_vec * p0_req)
ok  <- is.finite(p0_req) & p0_req > 0
fit <- lm(log(p0_req[ok]) ~ log(N0_vec[ok]))
coef(fit)  # slope ≈ -1 confirms p0* ~ 1/N0




# Core Rebound Time Analysis
# Focus on single population dynamics and rebound patterns

library(ggplot2)
library(dplyr)
library(scales)

# Set seed for reproducibility  
#set.seed(42)

# Source the original fitness function
fitness <- function(Ntot, w, f, K = 50000) {
  f_raw <- 1 + (w - 1) * (1 - f * Ntot / K)
  return(max(f_raw, 1e-9))
}

# Core simulation function - focused on rebound analysis
run_rebound_simulation <- function(N0, p0, w_a = 0.95, w_A = 1.05, f_a = 1.0, f_A, 
                                   K = 50000, tmax = 300, thresh_rescue = 0.5) {
  
  Na <- numeric(tmax)
  N_A <- numeric(tmax)
  Ntot <- numeric(tmax)
  
  Na[1] <- N0 * (1 - p0)
  N_A[1] <- N0 * p0
  Ntot[1] <- Na[1] + N_A[1]
  
  # Track key metrics
  rebound_time <- NA
  max_population <- Ntot[1]
  min_population_before_rebound <- Ntot[1]
  time_to_min <- 1
  seen_minimum <- FALSE
  
  for (t in 2:tmax) {
    wa_eff <- fitness(Ntot[t-1], w_a, f_a, K)
    wA_eff <- fitness(Ntot[t-1], w_A, f_A, K)
    
    Na[t] <- rpois(1, Na[t-1] * wa_eff)
    N_A[t] <- rpois(1, N_A[t-1] * wA_eff)
    Ntot[t] <- Na[t] + N_A[t]
    
    max_population <- max(max_population, Ntot[t])
    
    # Find minimum population before rebound
    if (!seen_minimum && t > 5) {
      if (Ntot[t] < min_population_before_rebound) {
        min_population_before_rebound <- Ntot[t]
        time_to_min <- t
      } else if (Ntot[t] > min_population_before_rebound * 1.1) {
        seen_minimum <- TRUE
      }
    }
    
    # Check for rebound (reaching rescue threshold)
    if (is.na(rebound_time) && Ntot[t] >= thresh_rescue * K) {
      rebound_time <- t
      break
    }
    
    # Stop if extinct
    if (Ntot[t] < 1) {
      break
    }
  }
  
  # Calculate final metrics
  final_population <- Ntot[t]
  final_generations <- t
  outcome <- if (final_population >= thresh_rescue * K) "rescued" else "extinct"
  
  list(
    time = 1:final_generations,
    total_population = Ntot[1:final_generations],
    wildtype = Na[1:final_generations],
    mutant = N_A[1:final_generations],
    rebound_time = rebound_time,
    min_population = min_population_before_rebound,
    time_to_min = time_to_min,
    max_population = max_population,
    final_population = final_population,
    outcome = outcome,
    total_generations = final_generations
  )
}

# === MAIN REBOUND ANALYSIS ===

cat("=== CORE REBOUND TIME ANALYSIS ===\n\n")

# Test across different p0 values to see rebound patterns
p0_values <- c(0.001, 0.002, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1)
rebound_analysis <- data.frame(
  p0 = p0_values,
  rebound_time = NA,
  min_population = NA,
  time_to_min = NA,
  outcome = "",
  final_population = NA,
  total_generations = NA
)

for (i in seq_along(p0_values)) {
  cat(sprintf("Testing p0 = %.3f...\n", p0_values[i]))
  
  result <- run_rebound_simulation(
    N0 = 10000, p0 = p0_values[i], w_a = 0.95, w_A = 1.05,
    f_a = 1.0, f_A = 2.0, K = 50000, tmax = 300
  )
  
  rebound_analysis$rebound_time[i] <- result$rebound_time
  rebound_analysis$min_population[i] <- result$min_population
  rebound_analysis$time_to_min[i] <- result$time_to_min
  rebound_analysis$outcome[i] <- result$outcome
  rebound_analysis$final_population[i] <- result$final_population
  rebound_analysis$total_generations[i] <- result$total_generations
}

print("Rebound Analysis Results:")
print(rebound_analysis)

# === VISUALIZATION: SINGLE POPULATION VIEW ===

# 1. Population trajectories showing rebound patterns
example_results <- list()
selected_p0 <- c(0.005, 0.01, 0.02, 0.05)

for (p0 in selected_p0) {
  result <- run_rebound_simulation(
    N0 = 10000, p0 = p0, w_a = 0.95, w_A = 1.05,
    f_a = 1.0, f_A = 2.0, K = 50000, tmax = 300
  )
  example_results[[as.character(p0)]] <- result
}

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

# 2. Rebound time vs initial mutant fraction
p_rebound <- ggplot(rebound_analysis, aes(x = p0, y = rebound_time, color = outcome)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(size = 1.2, alpha = 0.7) +
  scale_x_continuous(trans = "log10", labels = scientific_format()) +
  scale_y_continuous(limits = c(0, 300)) +
  scale_color_manual(values = c("rescued" = "#2E8B57", "extinct" = "#CD5C5C")) +
  labs(
    title = "Rebound Time vs Initial Mutant Fraction (p₀)",
    subtitle = "Time to reach rescue threshold (50% of carrying capacity)",
    x = "Initial Mutant Fraction (p₀, log scale)",
    y = "Rebound Time (generations)",
    color = "Outcome",
    caption = "Clear threshold separates rescue from extinction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

print(p_rebound)

# 3. Minimum population analysis
p_minimum <- ggplot(rebound_analysis, aes(x = p0, y = min_population, color = outcome)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(size = 1.2, alpha = 0.7) +
  scale_x_continuous(trans = "log10", labels = scientific_format()) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("rescued" = "#2E8B57", "extinct" = "#CD5C5C")) +
  labs(
    title = "Minimum Population vs Initial Mutant Fraction",
    subtitle = "Lowest population reached before recovery begins",
    x = "Initial Mutant Fraction (p₀, log scale)",
    y = "Minimum Population Size",
    color = "Outcome",
    caption = "Higher p₀ leads to higher minimum population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

print(p_minimum)

# === STATISTICAL ANALYSIS ===

cat("\n=== STATISTICAL ANALYSIS ===\n")

# Correlation analysis
correlation_p0_rebound <- cor(rebound_analysis$p0, rebound_analysis$rebound_time, 
                              method = "spearman", use = "complete.obs")

cat(sprintf("Spearman correlation between p₀ and rebound time: %.4f\n", correlation_p0_rebound))

# Analysis of rescued vs extinct populations
rescued_cases <- rebound_analysis[rebound_analysis$outcome == "rescued", ]
extinct_cases <- rebound_analysis[rebound_analysis$outcome == "extinct", ]

if (nrow(rescued_cases) > 0) {
  cat(sprintf("\nRescued populations (n=%d):\n", nrow(rescued_cases)))
  cat(sprintf("  Average rebound time: %.1f generations\n", 
              mean(rescued_cases$rebound_time, na.rm = TRUE)))
  cat(sprintf("  Range: %.1f to %.1f generations\n", 
              min(rescued_cases$rebound_time, na.rm = TRUE), 
              max(rescued_cases$rebound_time, na.rm = TRUE)))
  cat(sprintf("  Average minimum population: %.0f individuals\n", 
              mean(rescued_cases$min_population)))
}

if (nrow(extinct_cases) > 0) {
  cat(sprintf("\nExtinct populations (n=%d):\n", nrow(extinct_cases)))
  cat(sprintf("  Average minimum population: %.0f individuals\n", 
              mean(extinct_cases$min_population)))
  cat(sprintf("  Average final population: %.0f individuals\n", 
              mean(extinct_cases$final_population)))
}

# Critical threshold analysis
if (nrow(rescued_cases) > 0 && nrow(extinct_cases) > 0) {
  critical_p0 <- max(extinct_cases$p0)
  cat(sprintf("\nCritical threshold: p₀ ≈ %.3f separates rescue from extinction\n", critical_p0))
}

# Regression analysis for rescued populations only
if (nrow(rescued_cases) > 5) {
  # Log-log regression for power law relationship
  log_p0 <- log(rescued_cases$p0)
  log_rebound <- log(rescued_cases$rebound_time)
  
  regression_result <- lm(log_rebound ~ log_p0)
  cat(sprintf("\nPower law relationship for rescued populations:\n"))
  cat(sprintf("  log(rebound_time) = %.2f + %.2f × log(p₀)\n", 
              coef(regression_result)[1], coef(regression_result)[2]))
  cat(sprintf("  R² = %.4f, p-value = %.4f\n", 
              summary(regression_result)$r.squared, 
              summary(regression_result)$coefficients[2,4]))
}

# === KEY INSIGHTS FOR PRESENTATION ===

cat("\n=== KEY INSIGHTS FOR PRESENTATION ===\n")

cat("

# === SAVE RESULTS ===

# Save analysis results
results_summary <- list(
  rebound_analysis = rebound_analysis,
  correlation = correlation_p0_rebound,
  critical_threshold = ifelse(exists("critical_p0"), critical_p0, NA),
  rescued_cases = rescued_cases,
  extinct_cases = extinct_cases
)

saveRDS(results_summary, "core_rebound_analysis_results.rds")
write.csv(rebound_analysis, "core_rebound_analysis_data.csv", row.names = FALSE)

# Save plots
ggsave("population_trajectories.png", p_trajectories, width = 10, height = 6, dpi = 300)
ggsave("rebound_time_analysis.png", p_rebound, width = 10, height = 6, dpi = 300)
ggsave("minimum_population.png", p_minimum, width = 10, height = 6, dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Core rebound analysis completed successfully.\n")
cat("Files saved:\n")
cat("  - core_rebound_analysis_results.rds\n")
cat("  - core_rebound_analysis_data.csv\n")
cat("  - population_trajectories.png\n")
cat("  - rebound_time_analysis.png\n")
cat("  - minimum_population.png\n")

cat(sprintf("\nFinal Answer: The ideal mutant fraction is p₀ ≈ %.3f-%.3f\n", 
            min(rescued_cases$p0), max(rescued_cases$p0)))

write.csv(rebound_analysis,
          file = "~/Desktop/core_rebound_analysis_data.csv",   # any path you like
          row.names = FALSE)


## === QUICK PATCH: sweep f_A without touching .sim_once =====================
library(ggplot2)

fA_grid  <- seq(0.5, 3, by = 0.1)
nrep     <- 50                # keep small while testing
p0_fixed <- 0.01

rescue_vec <- numeric(length(fA_grid))

for(i in seq_along(fA_grid)){
  cat("f_A =", fA_grid[i], "\n")
  f_A <<- fA_grid[i]          # GLOBAL overwrite (<<-)
  rescue_vec[i] <- mean(replicate(nrep, .sim_once(p0 = p0_fixed)))
}

rescue_by_fA <- data.frame(f_A = fA_grid, rescue_prob = rescue_vec)

ggplot(rescue_by_fA, aes(f_A, rescue_prob)) +
  geom_line(size = 1.2, colour = "#0072B2") +
  geom_point(size = 3, colour = "#0072B2") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Rescue probability vs mutant density-sensitivity f_A",
       x = expression(f[A]),
       y = "Rescue probability") +
  theme_minimal(base_size = 14)
## === p₀* vs f_A  (p₀ that gives 90 % rescue) =============================
library(ggplot2)

fA_grid   <- seq(0.5, 3, by = 0.1)   # density-sensitivity sweep
nrep      <- 200                     # replicates per (f_A, p₀) combo
P_target  <- 0.9                     # 90 % rescue
thresh_rescue <- 0.5                 # same as in your script
K         <- 50000                   # carrying capacity

# log-spaced p₀ search grid (wide enough)
p0_search <- exp(seq(log(1 / K), log(0.2), length.out = 25))

# helper: find smallest p₀ that reaches P_target for a given f_A
p0_star <- function(fA_val){
  f_A <<- fA_val                         
  probs <- vapply(p0_search, function(p) mean(replicate(nrep, .sim_once(p))), numeric(1))
  idx   <- which(probs >= P_target)
  if(length(idx)) return(min(p0_search[idx])) else return(NA_real_)
}

# build the curve
p0_req_df <- data.frame(
  f_A = fA_grid,
  p0_star = vapply(fA_grid, p0_star, numeric(1))
)

# plot
ggplot(p0_req_df, aes(f_A, p0_star)) +
  geom_line(size = 1.2, colour = "#0072B2") +
  geom_point(size = 3, colour = "#0072B2") +
  scale_y_continuous(trans = "log10", labels = scales::scientific_format()) +
  labs(title = "Minimal initial mutant fraction for 90 % rescue",
       subtitle = "p₀* vs mutant density-sensitivity f_A",
       x = expression(f[A]),
       y = expression(p[0]^"*" ~ "(90 % rescue, log scale)")) +
  theme_minimal(base_size = 14)


library(ggplot2)
df <- data.frame(N0 = c(1e3, 1e4, 4e4),
                 n_mut = c(21, 37, 48))        # your real 90 % points

ggplot(df, aes(x = n_mut, y = factor(N0), fill = "Mutants needed")) +
  geom_col(width = 0.6, colour = "white") +
  geom_text(aes(label = sprintf("%d", n_mut)), hjust = -0.3) +
  scale_x_continuous(limits = c(0, 60), name = "Mutant individuals (90 % rescue)") +
  scale_y_discrete(name = "Starting population N₀") +
  scale_fill_manual(values = "#0072B2") +
  ggtitle("Head-count, not percentage, predicts rescue") +
  theme_minimal(14) +
  theme(legend.position = "none")
ggsave("headcount.png", width = 5, height = 3, dpi = 300)

