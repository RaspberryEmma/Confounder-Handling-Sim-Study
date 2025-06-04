# ****************************************
# Confounder Handling Simulation Study
#
# Random Number Generation
# Generates a CSV of random-number-generator seeds
# This ensures non-overlapping sequences between parallel jobs
#
# Emma Tarmey
#
# Started:          05/12/2024
# Most Recent Edit: 26/02/2025
# ****************************************

# clear R memory
rm(list=ls())

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# RNG State
set.seed(2024) # given from God
RNGkind(kind = "Mersenne-Twister")
RNGkind()
R.version

# Simulation parameters (upper limits where applicable)
n_run      <- 9
n_scen     <- 8
n_rep      <- 2000
n_obs      <- 10000
n_base_var <- 100  # number of variables with no prior cause (upper bound)

# Required samples for each individual job
n_samples_per_job <- as.integer( n_rep * n_obs * n_base_var )
print(n_samples_per_job)

# Required seeds
n_seeds <- n_run * n_scen
print(n_seeds)

# Check periodicity
RNG_period     <- (2^(19937) - 1)             # Mersenne-Twister period
min_req_period <- n_samples_per_job * n_seeds # Minimum period value needed to
                                              # guarantee zero overlap in simulation
message("\n\nChecking whether period of RNG algorithm is sufficient to produce non-overlapping sequences for all jobs")
print(RNG_period)
print(min_req_period)
message(ifelse((RNG_period > min_req_period), "YES", "NO"))

# Generate all seeds
seeds_df           <- data.frame(matrix(NA, nrow = n_seeds, ncol = 3)) 
colnames(seeds_df) <- c("simulation_run", "simulation_scenario", "seed")

for (i in 1:n_seeds) {
  # Extract next seed
  seed_sequence <- .Random.seed      # get current RNG state
  new_seed      <- seed_sequence[3]  # extract valid seed value
  
  message(paste("\n\nSeed", i))
  print(seed_sequence[c(1:10)])
  print(new_seed)
  
  seeds_df[i, "simulation_run"]      <- ceiling(i / n_scen)
  seeds_df[i, "simulation_scenario"] <- ifelse((i %% n_scen)==0, n_scen, (i %% n_scen)) 
  seeds_df[i, "seed"]                <- new_seed
  
  # generate and throw away 'n_samples_per_job'-many samples
  # sufficiently modifies RNG state
  set.seed(new_seed)
  runif(n_samples_per_job, min = 0, max = 1)
}

print(seeds_df)
write.csv(seeds_df, file = "../data/precomputed_RNG_seeds.csv", row.names = FALSE)


