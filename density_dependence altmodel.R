# --- one generation ---
simulation <- function(N_A, N_a, decay_rate, sel_coeff, mut_rate,
                       R0 = 1.2, K = 1000) {
  N_tot <- N_A + N_a
  if (N_tot == 0) return(c(a = 0, A = 0))
  
  p <- N_A / N_tot
  
  # frequency dependence (your original form)
  w_a <- 1 + sel_coeff * (p - 1)
  w_A <- 1 + sel_coeff * (1 - p)
  
  # ---- CARRYING CAPACITY here ----
  density_factor <- R0 / (1 + N_tot / K)
  
  lambda_a <- N_a * w_a * density_factor
  lambda_A <- N_A * w_A * density_factor
  
  # guard against negatives
  lambda_a <- max(0, lambda_a)
  lambda_A <- max(0, lambda_A)
  
  c(a = rpois(1, lambda_a),
    A = rpois(1, lambda_A))
}
simulation_pop <- function(N_init_a, N_init_A, decay_rate, sel_coeff, mut_rate,
                           t_max, R0 = 1.2, K = 1000) {
  pop_new <- c(a = N_init_a, A = N_init_A)
  pop_vec <- matrix(pop_new, nrow = 1, dimnames = list("0", c("a","A")))
  start_total <- sum(pop_new)
  
  for (i in 1:t_max) {
    pop_new <- simulation(
      N_A = pop_new["A"],  # A first
      N_a = pop_new["a"],  # then a
      decay_rate = decay_rate,
      sel_coeff  = sel_coeff,
      mut_rate   = mut_rate,
      R0 = R0, K = K       # <- pass K (and R0)
    )
    pop_vec <- rbind(pop_vec, pop_new)
    rownames(pop_vec)[nrow(pop_vec)] <- as.character(i)
    
    total_now <- sum(pop_new)
    if (total_now >= 50 * start_total || total_now == 0) break
  }
  pop_vec
}
set.seed(1)
out <- simulation_pop(100, 10, decay_rate = 0, sel_coeff = 0.5, mut_rate = 0,
                      t_max = 500, R0 = 1.2, K = 300)
gens   <- 0:(nrow(out)-1)
total  <- rowSums(out)

plot(gens, total, type = "l",
     xlab = "Generation", ylab = "Population size",
     ylim = c(0, max(total, 1.05*300)))   # 300 = K; adjust if different

lines(gens, out[,"a"], col = "blue")
lines(gens, out[,"A"], col = "red")

abline(h = 300, lty = 2)                  # horizontal line at K
legend("topright",
       legend = c("Total", "a", "A", "K"),
       lty = c(1,1,1,2), col = c("black","blue","red","black"),
       bty = "n")

