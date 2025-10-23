# Haploid population with differential density dependence
# Question:What is the ideal fraction of mutants to introduce, given variable density dependence,to rescue a declining population?
# Key Parameter(or Variable?): p0, initial fraction of mutants

## Parameters
N0  <- 100      # initial total population size
p0  <- 0.01      # initial fraction of mutants
K   <- 5000      # carrying capacity
w_a <- 0.95       # intrinsic fitness of wildtype (below 1)
w_A <- 1.05       # intrinsic fitness of mutant (above 1)
f_a <- 1.0        # wildtype density sensitivity
f_A <- 1.3        # mutant stronger negative freq. dependence
tmax <- 3000       # number of generations


## Initialization of first generation
Na <- numeric(tmax)    # wildtype counts
N_A <- numeric(tmax)   # mutant counts
Na[1] <- N0 * (1 - p0) # initialize wildtypes
N_A[1] <- N0 * p0      # initialize mutants


## Beverton–Holt crowding function
fitness <- function(Ntot, w, f) w / (1 + f * Ntot / K)  # always gives positive values (no negative growth)
# It smoothly saturates toward 0 as density increases.
# It’s biologically realistic, reflects limited resources, space, or competition.
# It allows different sensitivity between genotypes (f_a vs f_A).



## Simulation
set.seed(1)  # Each pass of the loop represents one generation 
# We start at t = 2 because the first generation (t = 1) was initialized already.
for(t in 2:tmax){
  Ntot <- Na[t-1] + N_A[t-1]   # We calculate the total population size in the previous generation.
  # This determines how crowded the population is, 
  # and so how much density dependence should reduce fitness.
  
  wa_eff <- fitness(Ntot, w_a, f_a)  # Each genotype gets its own densitydependence-ajusted growth rate for each generation.
  wA_eff <- fitness(Ntot, w_A, f_A)
  
  Na[t] <- rpois(1, Na[t-1] * wa_eff)  # Each genotype reproduces stochastically. 
  N_A[t] <- rpois(1, N_A[t-1] * wA_eff)# The expected number of offspring = previous number × realized fitness.
  # rpois() adds demographic stochasticity around that mean
  # example: 9000 × 0.8 = 7200, The realized number might be + - 7200
  
  if(Na[t] + N_A[t] < 1) break         # If total population drops below 1 (essentially extinct), stop the simulation.
}


## Plot population trajectories
time <- 1:tmax
plot(time, Na + N_A, type="l", lwd=2, col="black",               # total pop
     ylab="Population size", xlab="Generation",
     main="Haploid Rescue with Differential Density Dependence")
lines(time, Na, col="red", lty=2)                                # wildtype
lines(time, N_A, col="blue")                                     # mutant
legend("topright", legend=c("Total","Wildtype a","Mutant A"),    # legend
       col=c("black","red","blue"), lty=c(1,2,1), bty="n")

if(t %% 10 == 0) {   # print every 10 generations
  cat("Gen:", t, "Ntot =", Ntot, "wa_eff =", round(wa_eff,4), 
      "wA_eff =", round(wA_eff,4), "\n")
}

# p0 mainly changes time to rescue/extinction risk, not the final coexistence.
# ideal p0 is:
# (i) rescue probability
## If we repeat the stochastic simulation (e.g.100 replicates) for a given p0, some runs will go extinct, others will be rescued. The rescue probability is the percentage that survive.
# (ii) minimum N
## After introduction, the population plummets further until it’s low enough that the mutant, being more density-sensitive, can finally exploit the low-density environment and take over, rescuing the population.
# (iii) time to rebound
## Generation when growth resumes

# The ideal initial mutant fraction (p0) is the smallest proportion of mutants that reliably rescues the population. As p₀ increases, rescue success improves — the population is less likely to crash (higher minimum N) 
# and recovers faster (shorter rebound time). However, beyond a certain point, adding more mutants brings little additional benefit. The optimal p₀ is therefore at this plateau, where further increases no longer improve 
# rescue outcomes but would unnecessarily use more mutants or exceed carrying capacity.
# After a certain p₀, the mutant already guarantees rescue;
# further increases don’t change the outcome because the population is bounded by K and the mutant’s density-sensitive growth slows as soon as it becomes common.

fitness<- function(Ntot, w) 1 + w-w * Ntot / K

## Parameters
N0  <- 100     # initial total population size
p0  <- 0.1      # initial fraction of mutants
K   <- 500      # carrying capacity
w_a <- 0.95       # intrinsic fitness of wildtype (below 1)
w_A <- 1.05       # intrinsic fitness of mutant (above 1)
f_a <- 1.0        # wildtype density sensitivity
f_A <- 1.3        # mutant stronger negative freq. dependence
tmax <- 300       # number of generations


## Initialization of first generation
Na <- numeric(tmax)    # wildtype counts
N_A <- numeric(tmax)   # mutant counts
Na[1] <- N0 * (1 - p0) # initialize wildtypes
N_A[1] <- N0 * p0      # initialize mutants


## Beverton–Holt crowding function
fitness <- function(Ntot, w, f) 1+ w - w * Ntot / K  # always gives positive values (no negative growth)
# It smoothly saturates toward 0 as density increases.
# It’s biologically realistic, reflects limited resources, space, or competition.
# It allows different sensitivity between genotypes (f_a vs f_A).



## Simulation
set.seed(1)  # Each pass of the loop represents one generation of your haploid population. 
# We start at t = 2 because the first generation (t = 1) was initialized already.
for(t in 2:tmax){
  Ntot <- Na[t-1] + N_A[t-1]   # We calculate the total population size in the previous generation.
  # This determines how crowded the population is, 
  # and so how much density dependence should reduce fitness.
  
  #Each genotype gets its own densitydependence-ajusted growth rate for each generation.
 ## wa_eff <- fitness(Ntot, w_a, f_a)  
  wA_eff <- fitness(Ntot, w_A, f_A)
  wa_eff <- w_a  
  
  Na[t] <- rpois(1, Na[t-1] * wa_eff)  # Each genotype reproduces stochastically. 
  N_A[t] <- rpois(1, N_A[t-1] * wA_eff)# The expected number of offspring = previous number × realized fitness.
  # rpois() adds demographic stochasticity around that mean
  # example: 9000 × 0.8 = 7200, The realized number might be 7175, 7243, etc.
  
  if(Na[t] + N_A[t] < 1) break         # If total population drops below 1 (essentially extinct), stop the simulation.
}


## Plot population trajectories
time <- 1:tmax
plot(time, Na + N_A, type="l", lwd=2, col="black",               # total pop
     ylab="Population size", xlab="Generation",
     main="Haploid Rescue with Differential Density Dependence")
lines(time, Na, col="red", lty=2)                                # wildtype
lines(time, N_A, col="blue")                                     # mutant
legend("topright", legend=c("Total","Wildtype a","Mutant A"),    # legend
       col=c("black","red","blue"), lty=c(1,2,1), bty="n")

if(t %% 10 == 0) {   # print every 10 generations
  cat("Gen:", t, "Ntot =", Ntot, "wa_eff =", round(wa_eff,4), 
      "wA_eff =", round(wA_eff,4), "\n")
}

# p0 mainly changes time to rescue/extinction risk, not the final coexistence.
# ideal p0 is:
# (i) rescue probability
## If we repeat the stochastic simulation (e.g.100 replicates) for a given p0, some runs will go extinct, others will be rescued. The rescue probability is the percentage that survive.
# (ii) minimum N
## After introduction, the population plummets further until it’s low enough that the mutant, being more density-sensitive, can finally exploit the low-density environment and take over, rescuing the population.
# (iii) time to rebound
## Generation when growth resumes

# The ideal initial mutant fraction (p0) is the smallest proportion of mutants that reliably rescues the population. As p₀ increases, rescue success improves — the population is less likely to crash (higher minimum N) 
# and recovers faster (shorter rebound time). However, beyond a certain point, adding more mutants brings little additional benefit. The optimal p₀ is therefore at this plateau, where further increases no longer improve 
# rescue outcomes but would unnecessarily use more mutants or exceed carrying capacity.
# After a certain p₀, the mutant already guarantees rescue;
# further increases don’t change the outcome because the population is bounded by K and the mutant’s density-sensitive growth slows as soon as it becomes common.

print(simulate_gen(N0, p0, K, ))

