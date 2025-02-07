# ****************************************
# Confounder Handling Simulation Study
#
# Simplified, less flexible refactor of all existing sim code
#
# Emma Tarmey
#
# Started:          06/02/2025
# Most Recent Edit: 07/02/2025
# ****************************************



# ----- Tech Setup ----- 

# clear R memory
rm(list=ls())

# check all external libraries
using<-function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs, require, character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need) > 0){ 
    install.packages(need)
    lapply(need, require, character.only=TRUE)
  }
}
using("dagitty", "dplyr", "ggcorrplot", "ggplot2", "glmnet",
      "igraph", "lars", "MASS", "matrixStats", "microbenchmark",
      "questionr", "sjmisc", "tidyr")

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}



# ----- Random Number Generation -----

# fix RNG seed based on current run and scenario
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>%
  filter(simulation_run      == 1) %>%
  filter(simulation_scenario == 1)
set.seed(seed$seed)



# ----- Parameters ------

n_simulation      <- "TEST" # see Table!
n_scenario        <- "TEST" # see Table!
n_rep             <- 100 # 1000
n_obs             <- 10000
num_total_conf    <- 4
num_meas_conf     <- 4
num_unmeas_conf   <- 0
Z_correlation     <- 0.2
target_r_sq_X     <- 0.2
target_r_sq_Y     <- 0.4
causal            <- 0.5
dissimilarity_rho <- 1.0
binary_X          <- FALSE
binary_Y          <- FALSE
binary_Z          <- FALSE



# ----- Helper Functions -----


# Round all numeric columns in a given data frame to a given number of digits
# Credit to: https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
round_df <- function(df, digits) {
  nums      <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return (df)
}


# Extract list of covariates who are selected in a given model
find_vars_in_model <- function(model = NULL) {
  vars <- names(model$coefficients)[vars != "(Intercept)"]
  vars <- vars[vars != "(Intercept)"]
  return (vars)
}


# Create a string corresponding to a regression of Y on a given set of covariates
make_model_formula <- function(vars_selected = NULL) {
  formula_string <- "Y ~ X"
  for (var in vars_selected) {
    formula_string <- paste(formula_string, " + ", var, sep = "")
  }
  return(formula_string)
}


# Create a string corresponding to a regression of X on a given set of covariates
make_X_model_formula <- function(vars_selected = NULL) {
  if ('X' %in% vars_selected) {
    vars_selected <- vars_selected[! vars_selected%in% c('X')]
  }
  
  formula_string <- "X ~ 0"
  for (var in vars_selected) {
    formula_string <- paste(formula_string, " + ", var, sep = "")
  }
  return(formula_string)
}


# The value for all beta coefficients used for generating X
beta_X_formula <- function(num_total_conf = NULL,
                           target_r_sq_X  = NULL,
                           Z_correlation  = NULL) {
  
  numerator   <- target_r_sq_X
  denominator <- num_total_conf * (1 - target_r_sq_X) * (1 + Z_correlation*(num_total_conf - 1))
  
  value <- sqrt(numerator / denominator)
  
  return (value)
}


# The proportion of variance in X which is explained by measured variables (Zs)
# i.e Not including unmeasured confounding or error terms
r_squared_X <- function(X_model     = NULL,
                        X_test_data = NULL) {
  R2 <- NaN
  
  # separate exposure X from Z's
  test_Zs <- X_test_data[, -1, drop = F]
  test_X  <- X_test_data[, 1]
  
  # generate predicted value vector
  pred_X <- predict( X_model, test_Zs ) %>% as.vector()
  
  # calculate residual and total sum of squares
  SSR <- sum((pred_X - test_X)^2)
  SST <- sum((test_X - mean(test_X))^2)
  R2  <- (1 - (SSR / SST))
  
  return (R2)
}


# The proportion of variance in Y which is explained by measured variables (X and Zs)
# i.e Not including unmeasured confounding or error terms
r_squared_Y <- function(model     = NULL,
                        test_data = NULL) {
  R2 <- NaN
  
  # separate outcome Y from other covariates
  test_X <- test_data[, -1, drop = F]
  test_y <- test_data[,  1]
  
  # generate predicted value vector for each model type
  pred_y <- predict( model, test_X ) %>% as.vector()
  
  SSR <- sum((pred_y - test_y)^2)
  SST <- sum((test_y - mean(test_y))^2)
  R2  <- (1 - (SSR / SST))
  
  return (R2)
}


determine_var_error_Y <- function(num_total_conf = NULL,
                                  beta_X         = NULL,
                                  causal         = NULL,
                                  Z_correlation  = NULL,
                                  target_r_sq_Y  = NULL) {
  
  LHS <- (1 - target_r_sq_Y)/target_r_sq_Y
  
  S   <- (causal^2 + 1) * ((num_total_conf * beta_X^2) + (num_total_conf * (num_total_conf - 1) * Z_correlation * beta_X^2))
  RHS <- S + causal^2
  
  return (LHS * RHS)
}


determine_cov_matrix <- function(num_total_conf = NULL,
                                 var_names      = NULL,
                                 beta_X         = NULL,
                                 beta_Y         = NULL,
                                 causal         = NULL,
                                 Z_correlation  = NULL,
                                 target_r_sq_X  = NULL,
                                 target_r_sq_Y  = NULL) {
  
  num_vars         <- length(var_names)
  num_of_cov_terms <- ((num_total_conf) * (num_total_conf - 1))
  new_coef         <- ((causal*beta_X) + beta_Y)^2
  
  analytic_cov           <- matrix(NaN, num_vars, num_vars)
  rownames(analytic_cov) <- var_names
  colnames(analytic_cov) <- var_names
  
  var_error_Y <- determine_var_error_Y(num_total_conf = num_total_conf,
                                       beta_X         = beta_X,
                                       causal         = causal,
                                       Z_correlation  = Z_correlation,
                                       target_r_sq_Y  = target_r_sq_Y)
  
  # variances
  for (i in 1:num_vars) {
    # variances here
    if (var_names[i] == 'X') {
      analytic_cov[i, i] <- (num_total_conf * beta_X^2) + 1 + (num_of_cov_terms * beta_X^2 * Z_correlation)
    }
    else if (var_names[i] == 'Y') {
      analytic_cov[i, i] <- (num_total_conf * new_coef) + (causal^2) + var_error_Y + (num_of_cov_terms * new_coef * Z_correlation)
    }
    else {
      analytic_cov[i, i] <- (sqrt(Z_correlation))^2 + (sqrt(1 - Z_correlation))^2
    }
  }
  
  # pairwise covariances
  for (i in 1:num_vars) {
    for (j in 1:num_vars) {
      if (i == j) {
        # diagonal, skip
      }
      else if (var_names[i] == 'X' & var_names[j] == 'Y') {
        analytic_cov[i, j] <- NaN
      }
      else if (var_names[i] == 'Y' & var_names[j] == 'X') {
        analytic_cov[i, j] <- NaN
      }
      else if (var_names[i] == 'X' & var_names[j] != 'Y') {
        analytic_cov[i, j] <- NaN
      }
      else if (var_names[i] == 'Y' & var_names[j] != 'X') {
        analytic_cov[i, j] <- NaN
      }
      else if (var_names[i] != 'X' & var_names[j] == 'Y') {
        analytic_cov[i, j] <- NaN
      }
      else if (var_names[i] != 'Y' & var_names[j] == 'X') {
        analytic_cov[i, j] <- NaN
      }
      else {
        analytic_cov[i, j] <- Z_correlation
      }
    }
  }
  
  return (analytic_cov)
}



# ----- Data generation mechanism -----


# coefficient values for DAG
alpha <- sqrt(Z_correlation)     # U on Z
beta  <- sqrt(1 - Z_correlation) # error_Z on Z
beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                          target_r_sq_X  = target_r_sq_X,
                          Z_correlation  = Z_correlation)
beta_Y  <- beta_X # Z on Y

# variance of the error term for Y
var_error_Y <- determine_var_error_Y(num_total_conf = num_total_conf,
                                     beta_X         = beta_X,
                                     causal         = causal,
                                     Z_correlation  = Z_correlation,
                                     target_r_sq_Y  = target_r_sq_Y)


# datasets representative of our DAG
generate_dataset <- function() {
  # error terms
  error_Z1 <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_Z2 <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_Z3 <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_Z4 <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_X  <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_Y  <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
  
  # shared prior U for all Z_i
  prior_U <- rnorm(n = n_obs, mean = 0, sd = 1)
  
  # confounders Z
  Z1 <- (alpha * prior_U) + (beta * error_Z1)
  Z2 <- (alpha * prior_U) + (beta * error_Z2)
  Z3 <- (alpha * prior_U) + (beta * error_Z3)
  Z4 <- (alpha * prior_U) + (beta * error_Z4)
  
  # treatment variable X, outcome variable Y
  X  <- (beta_X * Z1) + (beta_X * Z2) + (beta_X * Z3) + (beta_X * Z4) + error_X
  Y  <- (causal * X)  + (beta_Y * Z1) + (beta_Y * Z2) + (beta_Y * Z3) + (beta_Y * Z4) + error_Y
  
  # combine into dataframe
  dataset <- as.data.frame(cbind(Y, X, Z1, Z2, Z3, Z4))
  
  return (dataset)
}



# ----- Model fitting and metric measurement -----


model_methods <- c("linear", "linear_unadjusted",
                   "stepwise", "stepwise_X",
                   "two_step_lasso", "two_step_lasso_X")

results_methods <- c("pred_mse", "r_squared_X", "r_squared_Y")

# results_methods <- c("pred_mse", "r_squared_X", "r_squared_Y",
#                      "model_SE", "emp_SE",
#                      "data_odds_ratio", "model_odds_ratio",
#                      "causal_true_val", "causal_effect_est",
#                      "causal_effect_bias", "causal_effect_mcse",
#                      "avg_abs_param_bias", "coverage",
#                      "open_paths", "blocked_paths")


var_names <- c("Y", "X", "Z1", "Z2", "Z3", "Z4")

n_variables <- length(var_names)
n_methods   <- length(model_methods)
n_results   <- length(results_methods)


# track all results measures across all repetitions
results <- array(data     = NaN,
                 dim      = c(n_methods, n_results, n_rep),
                 dimnames = list(model_methods, results_methods, 1:n_rep))

# track all coefficient values across all repetitions
model_coefs <- array(data     = NaN,
                     dim      = c(n_methods, n_variables, n_rep),
                     dimnames = list(model_methods, var_names, 1:n_rep) )

# track covariate set selections across all repetitions
cov_selection <- array(data     = NaN,
                       dim      = c(n_methods, (n_variables - 1), n_rep),
                       dimnames = list(model_methods, var_names[-c(1)], 1:n_rep) )


# run simulation repetitions
for (repetition in 1:n_rep) {
  
  # track progress
  message( paste("\n\nRunning Simulation Run ", n_simulation,
                 ", Scenario ", n_scenario,
                 ", Iteration ", repetition, "/", n_rep, "\n", sep = "") )
  
  # generate data
  dataset   <- generate_dataset()
  
  # cut-up versions of the data as needed
  X_dataset <- subset(dataset, select=-c(Y))
  Z_dataset <- subset(dataset, select=-c(Y, X))
  Y_column  <- subset(dataset, select=c(Y))
  X_column  <- subset(dataset, select=c(X))
  
  for (method in model_methods) {
  
    # fit model
    model   <- NULL
    X_model <- NULL # required to make R2X well-defined
    
    if (method == "linear") {
      model   <- lm("Y ~ .", data = dataset)
      X_model <- lm("X ~ .", data = X_dataset)
    }
    
    else if (method == "linear_unadjusted") {
      model   <- lm("Y ~ X", data = dataset)
      X_model <- lm("X ~ 0", data = X_dataset)
    }
    
    else if (method == "stepwise") {
      stepwise_model      <- step(object    = lm("Y ~ .", data = dataset),            # all variable base
                                  direction = "both",                                 # stepwise, not fwd or bwd
                                  scope     = list(upper = "Y ~ .", lower = "Y ~ X"), # exposure X always included
                                  trace     = 0)                                      # suppress output
      
      vars_selected <- names(stepwise_model$coefficients)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      model   <- lm(model_formula,   data = dataset)
      X_model <- lm(X_model_formula, data = X_dataset)
    }
    
    else if (method == "stepwise_X") {
      stepwise_X_model   <- step(object    = lm("X ~ .", data = X_dataset),            # all variable base
                                 direction = "both",                                 # stepwise, not fwd or bwd
                                 scope     = list(upper = "X ~ .", lower = "X ~ 0"), # constant term
                                 trace     = 0)
      
      vars_selected <- names(stepwise_X_model$coefficients)
      vars_selected <- c('X', vars_selected)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      model   <- lm(model_formula,   data = dataset)
      X_model <- lm(X_model_formula, data = X_dataset)
    }
    
    else if (method == "two_step_lasso") {
      cv_lasso_model <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1)
      lambda         <- cv_lasso_model$lambda.min
      lasso_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1, lambda=lambda)
      
      vars_selected <- rownames(lasso_model$beta)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      model   <- lm(model_formula,   data = dataset)
      X_model <- lm(X_model_formula, data = X_dataset)
    }
    
    else if (method == "two_step_lasso_X") {
      cv_lasso_X_model <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1)
      lambda           <- cv_lasso_model$lambda.min
      lasso_model      <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1, lambda=lambda)
      
      vars_selected <- rownames(lasso_model$beta)
      vars_selected <- c('X', vars_selected)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      model   <- lm(model_formula,   data = dataset)
      X_model <- lm(X_model_formula, data = X_dataset)
    }
    
    # record metrics
    results[ method, "pred_mse", repetition]    <- mean(model$residuals^2)
    results[ method, "r_squared_X", repetition] <- r_squared_X(X_model = X_model, X_test_data = X_dataset)
    results[ method, "r_squared_Y", repetition] <- r_squared_Y(model = model, test_data = dataset)
  }
}


# ----- Process and save results -----

 # Coefficients fitted and error-variance fitted
coefs <- c(beta_X, beta_Y, causal)
names(coefs) <- c("beta_X", "beta_Y", "causal")

message("\n\nTrue Coefficients of DAG and Variance of error of Y")
print(coefs)
print(var_error_Y)

# Covariance Matrices
analytic_cov_matrix <- determine_cov_matrix(num_total_conf = num_total_conf,
                                            var_names      = var_names,
                                            beta_X         = beta_X,
                                            beta_Y         = beta_Y,
                                            causal         = causal,
                                            Z_correlation  = Z_correlation,
                                            target_r_sq_X  = target_r_sq_X,
                                            target_r_sq_Y  = target_r_sq_Y)
message("\n\nAnalytic Covariance:")
print(analytic_cov_matrix)

observed_cov_matrix <- round_df(as.data.frame(cov(dataset)), digits=3)
message("\n\nObserved Covariance:")
print(observed_cov_matrix)

# Results Table
final_results <- round_df(as.data.frame(apply(results, c(1,2), mean)), digits=3)
message("\n\nResults Table:")
print(final_results)
write.csv(final_results, "../data/final_results.csv")


