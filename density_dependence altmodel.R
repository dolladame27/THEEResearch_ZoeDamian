## ---------- Lynx two-type model ----------

# One time step (generation)
step_lynx <- function(N_a, N_A, R_a = 0.9, R_A = 1.15,
                      s_a = 0.2, s_A = 0.2,
                      K = 500,
                      stochastic = FALSE) {
  
  Ntot <- N_a + N_A
  if (Ntot <= 0) return(c(a = 0, A = 0))
  
  pA <- N_A / Ntot              # frequency of A
  # frequency effects (centered at 0.5; change if you want a different midpoint)
  w_a <- 1 + s_a * (pA - 0.5)  # a has NEGATIVE freq dep. (rarer a -> higher pA -> higher w_a)
  w_A <- 1 + s_A * (pA - 0.5)  # A has POSITIVE  freq dep. (higher pA -> higher w_A)
  
  # keep multipliers non-negative
  w_a <- max(0, w_a)
  w_A <- max(0, w_A)
  
  # Beverton–Holt density regulation
  D <- 1 / (1 + Ntot / K)
  
  # per-type expected next counts
  lambda_a <- N_a * R_a * w_a * D
  lambda_A <- N_A * R_A * w_A * D
  
  if (!stochastic) {
    N_a1 <- lambda_a
    N_A1 <- lambda_A
  } else {
    # Poisson demographic noise
    N_a1 <- rpois(1, max(0, lambda_a))
    N_A1 <- rpois(1, max(0, lambda_A))
  }
  
  c(a = N_a1, A = N_A1)
}

# Simulate many generations
simulate_lynx <- function(Na0 = 100, NA0 = 10,
                          R_a = 0.9, R_A = 1.15,
                          s_a = 0.2, s_A = 0.2,
                          K = 500,
                          t_max = 300,
                          stochastic = FALSE) {
  out <- matrix(NA_real_, nrow = t_max + 1, ncol = 2,
                dimnames = list(0:t_max, c("a","A")))
  out[1, ] <- c(a = Na0, A = NA0)
  
  for (t in 1:t_max) {
    nxt <- step_lynx(out[t, "a"], out[t, "A"],
                     R_a = R_a, R_A = R_A,
                     s_a = s_a, s_A = s_A,
                     K = K, stochastic = stochastic)
    out[t + 1, ] <- nxt
    if (sum(nxt) <= 0) break  # extinct; remaining rows stay N_A
  }
  out
}

# Quick plot helper
plot_lynx <- function(out, K, main = "Lynx population dyN_Amics") {
  gens  <- as.numeric(rownames(out))
  valid <- complete.cases(out)
  gens  <- gens[valid]; out <- out[valid, , drop = FALSE]
  
  total <- rowSums(out)
  plot(gens, total, type = "l",
       xlab = "Generation", ylab = "Population size",
       ylim = c(0, max(total, K)),
       main = main)
  lines(gens, out[, "a"], col = "blue")
  lines(gens, out[, "A"], col = "red")
  abline(h = K, lty = 2)
  legend("topright",
         legend = c("Total", "a (wildtype)", "A (introduced)", "K"),
         lty = c(1,1,1,2), col = c("black","blue","red","black"),
         bty = "n")
}

## ---------- Examples ----------

# Deterministic run
out_det <- simulate_lynx(Na0 = 100, NA0 = 10,
                         R_a = 0.9, R_A = 1.15,
                         s_a = 0.2, s_A = 0.2,
                         K = 500, t_max = 300,
                         stochastic = FALSE)
plot_lynx(out_det, K = 500, main = "Deterministic")

# Stochastic run
set.seed(42)
out_sto <- simulate_lynx(Na0 = 100, NA0 = 10,
                         R_a = 0.9, R_A = 1.15,
                         s_a = 0.2, s_A = 0.2,
                         K = 500, t_max = 300,
                         stochastic = TRUE)
plot_lynx(out_sto, K = 500, main = "Stochastic")



