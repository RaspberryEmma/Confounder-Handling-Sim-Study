# ****************************************
# Confounder Handling Simulation Study
#
# Snippets
# 
#
# Emma Tarmey
#
# Started:          10/02/2025
# Most Recent Edit: 14/03/2025
# ****************************************



# # datasets representative of our DAG
# generate_dataset <- function() {
#   # error terms
#   error_Z1 <- rnorm(n = n_obs, mean = 0, sd = 1)
#   error_Z2 <- rnorm(n = n_obs, mean = 0, sd = 1)
#   error_Z3 <- rnorm(n = n_obs, mean = 0, sd = 1)
#   error_Z4 <- rnorm(n = n_obs, mean = 0, sd = 1)
#   error_X  <- rnorm(n = n_obs, mean = 0, sd = 1)
#   error_Y  <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
#   
#   # shared prior U for all Z_i
#   prior_U <- rnorm(n = n_obs, mean = 0, sd = 1)
#   
#   # confounders Z
#   Z1 <- (alpha * prior_U) + (beta * error_Z1)
#   Z2 <- (alpha * prior_U) + (beta * error_Z2)
#   Z3 <- (alpha * prior_U) + (beta * error_Z3)
#   Z4 <- (alpha * prior_U) + (beta * error_Z4)
#   
#   # treatment variable X, outcome variable Y
#   X  <- (beta_X * Z1) + (beta_X * Z2) + (beta_X * Z3) + (beta_X * Z4) + error_X
#   Y  <- (causal * X)  + (beta_Y * Z1) + (beta_Y * Z2) + (beta_Y * Z3) + (beta_Y * Z4) + error_Y
#   
#   # combine into dataframe
#   dataset <- as.data.frame(cbind(Y, X, Z1, Z2, Z3, Z4))
#   
#   return (dataset)
# }



# # Prevalences of binary variables
# message("\n\nMeans of X, Y, Zi, equal to prevalence if binary")
# 
# message("\nX:")
# print(binary_X_prevalence)
# print(summary(dataset$X))
# print(var(dataset$X))
# 
# message("\nY:")
# print(binary_Y_prevalence)
# print(summary(dataset$Y))
# print(var(dataset$Y))
# 
# message("\nZ1:")
# print(binary_Z_prevalence)
# print(summary(dataset$Z1))
# print(var(dataset$Z1))
# 
# message("\n\nData:")
# print(head(dataset))




