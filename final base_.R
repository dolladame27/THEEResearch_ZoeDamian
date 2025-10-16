# Haploid population with differential density dependence
# Question: What is the ideal fraction of mutants (p0) to introduce
# to rescue a declining population?

## Parameters
N0   <- 10000     # initial total population size
p0   <- 0.01      # initial fraction of mutants
K    <- 50000     # carrying capacity
w_a  <- 0.95      # intrinsic fitness of wildtype (below 1)
w_A  <- 1.05      # intrinsic fitness of mutant (above 1)
tmax <- 300       # number of generations

## Initialization
Na  <- numeric(tmax)          # wildtype counts
N_A <- numeric(tmax)          # mutant counts
Na[1]  <- N0 * (1 - p0)       # initialize wildtypes
N_A[1] <- N0 * p0             # initialize mutants

## Logistic fitness function
fitness <- function(Ntot, w) {
  f_raw <- 1 + (w - 1) * (1 - Ntot / K)  # logistic density dependence
  pmax(f_raw, 1e-9)                      # ensure fitness stays positive
}

## Simulation
set.seed(1)
t_last <- 1

for (t in 2:tmax) {
  Ntot <- Na[t-1] + N_A[t-1]             # total population at previous gen
  
  # Effective density-dependent fitness each generation
  wa_eff <- fitness(Ntot, w_a)
  wA_eff <- fitness(Ntot, w_A)
  
  # Stochastic reproduction
  Na[t]  <- rpois(1, Na[t-1]  * wa_eff)
  N_A[t] <- rpois(1, N_A[t-1] * wA_eff)
  
  t_last <- t
  if ((Na[t] + N_A[t]) < 1) break        # stop if population extinct
  
  # Print every 10 generations
  if (t %% 10 == 0) {
    cat("Gen:", t, 
        "Ntot =", round(Ntot), 
        "wa_eff =", round(wa_eff, 4), 
        "wA_eff =", round(wA_eff, 4), "\n")
  }
}

## Plot results
time <- 1:t_last
plot(time, (Na + N_A)[time], type="l", lwd=2, col="black",
     ylab="Population size", xlab="Generation",
     main="Haploid Rescue with Logistic Density Dependence")
lines(time, Na[time],  col="red",  lty=2, lwd=1.5)
lines(time, N_A[time], col="blue", lwd=1.5)
legend("topright", legend=c("Total","Wildtype a","Mutant A"),
       col=c("black","red","blue"), lty=c(1,2,1), lwd=c(2,1.5,1.5), bty="n")


# p0 mainly changes time to rescue/extinction risk, not the final coexistence.
# ideal p0 is:
# (i) rescue probability, maximal
## If we repeat the stochastic simulation (e.g.100 replicates) for a given p0, some runs will go extinct, others will be rescued. The rescue probability is the percentage that survive.
# (ii) minimum N, maximal
## After introduction, the population plummets further until it’s low enough that the mutant, being more density-sensitive, can finally exploit the low-density environment and take over, rescuing the population.
# (iii) time to rebound, minimal
## Generation when growth resumes

# The ideal initial mutant fraction (p0) is the smallest proportion of mutants that reliably rescues the population. As p₀ increases, rescue success improves — the population is less likely to crash (higher minimum N) 
# and recovers faster (shorter rebound time). However, beyond a certain point, adding more mutants brings little additional benefit. The optimal p₀ is therefore at this plateau, where further increases no longer improve 
# rescue outcomes but would unnecessarily use more mutants or exceed carrying capacity.
# After a certain p₀, the mutant already guarantees rescue;
# further increases don’t change the outcome because the population is bounded by K and the mutant’s density-sensitive growth slows as soon as it becomes common.

