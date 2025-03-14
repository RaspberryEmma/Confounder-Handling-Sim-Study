# ****************************************
# Confounder Handling Simulation Study
#
# Simplified, less flexible refactor of all existing sim code
#
# Emma Tarmey
#
# Started:          06/02/2025
# Most Recent Edit: 14/03/2025
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
using("dplyr", "glmnet", "tidyr")

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}



# ----- Parameters ------

# Run
n_simulation      <- "TEST" # see Table!

n_obs             <- 1000 # 10000
n_rep             <- 2 # 2000
Z_correlation     <- 0.2
Z_subgroups       <- 4.0
target_r_sq_X     <- 0.1
target_r_sq_Y     <- 0.4
causal            <- 0.0 # 0.15 # 0.5

binary_X            <- TRUE
binary_X_prevalence <- 0.30
binary_Y            <- TRUE
binary_Y_prevalence <- 0.05 # 0.30
binary_Z            <- TRUE
binary_Z_prevalence <- 0.30

# Scenario
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  n_scenario      <- as.numeric(args[1])
  num_total_conf  <- as.numeric(args[2])
  num_meas_conf   <- as.numeric(args[3])
  num_unmeas_conf <- as.numeric(args[4])
  
} else {
  n_scenario      <- "TEST"
  num_total_conf  <- 16
  num_meas_conf   <- 16
  num_unmeas_conf <- 0
}



# ----- Random Number Generation -----

# fix RNG seed based on current run and scenario
# NB: Overwritten:
# filter(simulation_run      == n_simulation) %>%
# filter(simulation_scenario == n_scenario)
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>%
  filter(simulation_run      == 1) %>%
  filter(simulation_scenario == 1)
set.seed(seed$seed)



# ----- Helper Functions -----


# Inverse-logit transform function
# Used for simulating binary data
inverse_logit <- function(real_values = NULL) {
  probabilities <- (1)/(1 + exp(-1 * real_values))
  return (probabilities)
}


# Round all numeric columns in a given data frame to a given number of digits
# Credit to: https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
round_df <- function(df, digits) {
  nums      <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- signif(df[,nums], digits = digits)
  return (df)
}


# Map of "is.nan" for different object types
# Credit to: https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}


# Create a string corresponding to a regression of Y on a given set of covariates
make_model_formula <- function(vars_selected = NULL) {
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("Y ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "Y ~ 0"
  }
  
  return(formula_string)
}


# Create a string corresponding to a regression of X on a given set of covariates
make_X_model_formula <- function(vars_selected = NULL) {
  # remove X
  if ('X' %in% vars_selected) {
    vars_selected <- vars_selected[! vars_selected%in% c('X')]
  }
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("X ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "X ~ 0"
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


# The value for the beta coefficients used for generating X
# 4 different values corresponding to the 4 subgroups
beta_X_subgroups_formula <- function(beta_X = NULL) {
  beta_X_1 <- beta_X     # HH
  beta_X_2 <- beta_X     # HL
  beta_X_3 <- beta_X / 4 # LH
  beta_X_4 <- beta_X / 4 # LL
  
  beta_Xs <- c(beta_X_1, beta_X_2, beta_X_3, beta_X_4)
  
  return (beta_Xs)
}


# The value for the beta coefficients used for generating Y
# 4 different values corresponding to the 4 subgroups
beta_Y_subgroups_formula <- function(beta_X = NULL) {
  beta_Y_1 <- beta_X     # HH
  beta_Y_2 <- beta_X / 4 # HL
  beta_Y_3 <- beta_X     # LH
  beta_Y_4 <- beta_X / 4 # LL
  
  beta_Ys <- c(beta_Y_1, beta_Y_2, beta_Y_3, beta_Y_4)
  
  return (beta_Ys)
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


# Estimate the variance in Y before error
# Idea is that we fix this for given values of R2 and b
# We then get a compatible error term later on
determine_subgroup_var_Y <- function(num_total_conf = NULL,
                                     beta_X         = NULL,
                                     causal         = NULL,
                                     Z_correlation  = NULL,
                                     target_r_sq_Y  = NULL) {
  
  subgroup_size <- num_total_conf / 4
  beta_X_high   <- beta_X
  beta_X_low    <- beta_X / 4
  
  A <- (causal * beta_X_low)  + beta_X_high
  B <- (causal * beta_X_high) + beta_X_high
  C <- (causal * beta_X_high) + beta_X_low
  D <- (causal * beta_X_low)  + beta_X_low
  
  pairwise_combinations <- (A*B) + (A*C) + (A*D) + (B*C) + (B*D) + (C*D)
  
  LHS <- (A^2 + B^2 + C^2 + D^2) * (subgroup_size + (Z_correlation * subgroup_size * (subgroup_size - 1)))
  MID <- 2 * Z_correlation * subgroup_size^2 * pairwise_combinations
  RHS <- causal^2
  
  value <- LHS + MID + RHS
  
  return (value)
}


# Estimate the variance in the error term for Y
# Takes a given variance of Y and inverts the R2 formula
determine_subgroup_var_error_Y <- function(var_Y          = NULL,
                                           target_r_sq_Y  = NULL) {
  
  LHS   <- 1 - target_r_sq_Y
  RHS   <- var_Y / target_r_sq_Y
  value <- sqrt(LHS * RHS)
  
  return (value)
}


# Estimate covariance matrix of entire system
# Structurally same as above but with formulae swapped out
# Specifically this function gives the more-complicated "subgroups" formulae for variances
determine_subgroup_cov_matrix <- function(num_total_conf = NULL,
                                          var_names      = NULL,
                                          beta_X         = NULL,
                                          causal         = NULL,
                                          Z_correlation  = NULL,
                                          target_r_sq_X  = NULL,
                                          target_r_sq_Y  = NULL) {
  
  num_vars         <- length(var_names)
  num_of_cov_terms <- ((num_total_conf) * (num_total_conf - 1))
  subgroup_size    <- num_total_conf / 4
  
  analytic_cov           <- matrix(NaN, num_vars, num_vars)
  rownames(analytic_cov) <- var_names
  colnames(analytic_cov) <- var_names
  
  var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                    beta_X         = beta_X,
                                    causal         = causal,
                                    Z_correlation  = Z_correlation,
                                    target_r_sq_Y  = target_r_sq_Y)
  
  var_error_Y <- determine_subgroup_var_error_Y(var_Y          = var_Y,
                                                target_r_sq_Y  = target_r_sq_Y)
  
  real_var_Y <- var_Y + var_error_Y
  
  # variances
  # NB: i indexes over {1, ..., m, m+1, m+2} for variables {Y, X, Z1, ..., Zm}
  for (i in 1:num_vars) {
    # variances here
    if (var_names[i] == 'X') {
      analytic_cov[i, i] <- ((num_total_conf / 2) * (beta_X^2 / 16)) + ((num_total_conf / 2) * beta_X^2) + 1
    }
    else if (var_names[i] == 'Y') {
      analytic_cov[i, i] <- real_var_Y
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


# Helper function for recording coefficients fitted
# Accounts for the idea that some modelling methods will exclude variables
# i.e:  Variables correct order with NaNs filling in excluded values
fill_in_blanks <- function(coefs = NULL, labels = NULL) {
  for (label in labels) {
    # if a given variable doesn't exist, create as NaN
    if (!(label %in% names(coefs))) {
      coefs[label] <- NaN
    }
  }
  
  # assert variable ordering
  coefs <- coefs[order(factor(names(coefs), levels = labels))]
  
  return (coefs)
}



# ----- Data generation mechanism -----

# coefficient values for DAG
alpha <- sqrt(Z_correlation)     # U on Z
beta  <- sqrt(1 - Z_correlation) # error_Z on Z

beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                          target_r_sq_X  = target_r_sq_X,
                          Z_correlation  = Z_correlation)

beta_Xs <- beta_X_subgroups_formula(beta_X = beta_X)
beta_Ys <- beta_Y_subgroups_formula(beta_X = beta_X)

# variance of the error term for Y
var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                  beta_X         = beta_X,
                                  causal         = causal,
                                  Z_correlation  = Z_correlation,
                                  target_r_sq_Y  = target_r_sq_Y)

var_error_Y <- determine_subgroup_var_error_Y(var_Y          = var_Y,
                                              target_r_sq_Y  = target_r_sq_Y)


# datasets representative of our DAG
generate_dataset <- function() {
  num_of_batches       <- num_total_conf / 4
  num_censored_batches <- num_unmeas_conf / 4
  
  dataset <- data.frame(matrix(NaN, nrow = n_obs, ncol = length(var_names)))
  colnames(dataset) <- var_names
  
  # shared prior U for all Z_i
  prior_U <- rnorm(n = n_obs, mean = 0, sd = 1)
  
  X <- rep(0, length.out = n_obs)
  Y <- rep(0, length.out = n_obs)
  
  # add effect of confounders Z_i to X and Y in batches of 4
  if (!binary_Z) {
    for (i in 1:num_of_batches) {
      # error terms
      error_Z1 <- rnorm(n = n_obs, mean = 0, sd = 1)
      error_Z2 <- rnorm(n = n_obs, mean = 0, sd = 1)
      error_Z3 <- rnorm(n = n_obs, mean = 0, sd = 1)
      error_Z4 <- rnorm(n = n_obs, mean = 0, sd = 1)
      
      # confounders Z
      Z1 <- (alpha * prior_U) + (beta * error_Z1) # always X=L, Y=L
      Z2 <- (alpha * prior_U) + (beta * error_Z2) # always X=L, Y=H
      Z3 <- (alpha * prior_U) + (beta * error_Z3) # always X=H, Y=L
      Z4 <- (alpha * prior_U) + (beta * error_Z4) # always X=H, Y=H
      
      # add confounder effect on treatment variable X and outcome variable Y
      X  <- X + (beta_Xs[1] * Z1) + (beta_Xs[2] * Z2) + (beta_Xs[3] * Z3) + (beta_Xs[4] * Z4)
      Y  <- Y + (beta_Ys[1] * Z1) + (beta_Ys[2] * Z2) + (beta_Ys[3] * Z3) + (beta_Ys[4] * Z4)
      
      # record Zs
      # index formula places Zs such that all Zs of one sub-group are next to each other
      dataset[, (2 + (0 * num_of_batches) + i)] <- Z1
      dataset[, (2 + (1 * num_of_batches) + i)] <- Z2
      dataset[, (2 + (2 * num_of_batches) + i)] <- Z3
      dataset[, (2 + (3 * num_of_batches) + i)] <- Z4
    }
  }
  else {
    if (Z_correlation == 0.0) {
      # uncorrelated binary case
      for (i in 1:num_of_batches) {
        # independent error terms
        uniform_prior_U_Z1 <- runif(n = n_obs, min = 0, max = 1)
        uniform_prior_U_Z2 <- runif(n = n_obs, min = 0, max = 1)
        uniform_prior_U_Z3 <- runif(n = n_obs, min = 0, max = 1)
        uniform_prior_U_Z4 <- runif(n = n_obs, min = 0, max = 1)
        
        # confounders Z
        Z1 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U_Z1)) - 0.94
        Z2 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U_Z2)) - 0.94
        Z3 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U_Z3)) - 0.94
        Z4 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U_Z4)) - 0.94
        
        # add confounder effect on treatment variable X and outcome variable Y
        X  <- X + (beta_Xs[1] * Z1) + (beta_Xs[2] * Z2) + (beta_Xs[3] * Z3) + (beta_Xs[4] * Z4)
        Y  <- Y + (beta_Ys[1] * Z1) + (beta_Ys[2] * Z2) + (beta_Ys[3] * Z3) + (beta_Ys[4] * Z4)
        
        # record Zs
        # index formula places Zs such that all Zs of one sub-group are next to each other
        dataset[, (2 + (0 * num_of_batches) + i)] <- Z1
        dataset[, (2 + (1 * num_of_batches) + i)] <- Z2
        dataset[, (2 + (2 * num_of_batches) + i)] <- Z3
        dataset[, (2 + (3 * num_of_batches) + i)] <- Z4
      }
    }
    else {
      # correlated binary case
      uniform_prior_U <- runif(n = n_obs, min = 0, max = 1)
      
      for (i in 1:num_of_batches) {
        # confounders Z
        Z1 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U)) - 0.94
        Z2 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U)) - 0.94
        Z3 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U)) - 0.94
        Z4 <- 2.06 * rbinom(n = n_obs, size = 1, prob = inverse_logit(uniform_prior_U)) - 0.94
        
        # add confounder effect on treatment variable X and outcome variable Y
        X  <- X + (beta_Xs[1] * Z1) + (beta_Xs[2] * Z2) + (beta_Xs[3] * Z3) + (beta_Xs[4] * Z4)
        Y  <- Y + (beta_Ys[1] * Z1) + (beta_Ys[2] * Z2) + (beta_Ys[3] * Z3) + (beta_Ys[4] * Z4)
        
        # record Zs
        # index formula places Zs such that all Zs of one sub-group are next to each other
        dataset[, (2 + (0 * num_of_batches) + i)] <- Z1
        dataset[, (2 + (1 * num_of_batches) + i)] <- Z2
        dataset[, (2 + (2 * num_of_batches) + i)] <- Z3
        dataset[, (2 + (3 * num_of_batches) + i)] <- Z4
      }
    }
  }
  
  # add error term to X if not binary
  if (!binary_X) {
    error_X <- rnorm(n = n_obs, mean = 0, sd = 1)
    X <- X + error_X
  }
  
  # add error term to Y if not binary
  if (!binary_Y) {
    error_Y <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
    Y <- Y + error_Y
  }
  
  # binarize X if binary
  # NB: R2X = 0.6 here for binary
  if (binary_X) {
    # NB: intercept term of logit expression controls prevalence (mean) of binary var
    logit_prob_X  <- X - 1.40                                   # interpret existing values as logit(probability)
    prob_X        <- inverse_logit(logit_prob_X)                # apply inverse to obtain prob values
    binary_vals_X <- rbinom(n = n_obs, size = 1, prob = prob_X) # re-sample to obtain X
    X             <- binary_vals_X                              # write binary values over previous continuous values
  }
  
  # add causal effect (X on Y)
  # NB: causal = 0.15 for binary
  Y <- Y + (causal * X) 
  
  # binarize Y if binary
  if (binary_Y) {
    # NB: intercept term of logit expression controls prevalence (mean) of binary var
    # common Y: intercept = 1.45; rare Y: intercept = 3.71
    logit_prob_Y  <- Y - 1.45                                    # interpret existing values as logit(probability)
    prob_Y        <- inverse_logit(logit_prob_Y)                # apply inverse to obtain prob values
    binary_vals_Y <- rbinom(n = n_obs, size = 1, prob = prob_Y) # re-sample to obtain Y
    Y             <- binary_vals_Y                              # write binary values over previous continuous values
  }
  
  # record X and Y
  dataset[, 'Y'] <- Y
  dataset[, 'X'] <- X
  
  # censor covariates as appropriate
  if (num_censored_batches > 0) {
    for (i in 1:num_censored_batches) {
      # index formula censors Z's fairly between all subgroups
      drop <- paste('Z', c((0 * num_of_batches) + i,
                           (1 * num_of_batches) + i,
                           (2 * num_of_batches) + i,
                           (3 * num_of_batches) + i), sep='')
      
      # censor Zs
      dataset <- dataset[, !(names(dataset) %in% drop)]
    }
  }
  
  return (dataset)
}



# ----- Model fitting and metric measurement -----


model_methods <- c("linear", "linear_unadjusted",
                   "stepwise", "stepwise_X",
                   "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")

results_methods <- c("pred_mse", "model_SE", "emp_SE",
                     "r_squared_X", "r_squared_Y",
                     "causal_true_value", "causal_estimate", "causal_bias", "causal_coverage",
                     "open_paths", "blocked_paths")

var_names                         <- c("Y", "X", paste('Z', c(1:num_total_conf), sep=''))
var_names_except_Y                <- var_names[ !var_names == 'Y']
var_names_except_Y_with_intercept <- c("(Intercept)", var_names_except_Y)

n_variables <- length(var_names)
n_methods   <- length(model_methods)
n_results   <- length(results_methods)


# track all results measures across all repetitions
results <- array(data     = NaN,
                 dim      = c(n_methods, n_results, n_rep),
                 dimnames = list(model_methods, results_methods, 1:n_rep))

# track all coefficient values across all repetitions
model_coefs <- array(data     = NaN,
                     dim      = c(n_methods, (n_variables), n_rep),
                     dimnames = list(model_methods, var_names_except_Y_with_intercept, 1:n_rep) )

# track covariate set selections across all repetitions
cov_selection <- array(data     = NaN,
                       dim      = c(n_methods, (n_variables - 1), n_rep),
                       dimnames = list(model_methods, var_names_except_Y, 1:n_rep) )


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
      if (binary_Y) {
        model <- glm("Y ~ .", data = dataset, family = "binomial")
      }
      else {
        model <- lm("Y ~ .",  data = dataset)
      }
      if (binary_X) {
        X_model <- glm("X ~ .", data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm("X ~ .",  data = X_dataset)
      }
      
      vars_selected <- names(model$coefficients)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
    }
    
    else if (method == "linear_unadjusted") {
      if (binary_Y) {
        model <- glm("Y ~ X", data = dataset, family = "binomial")
      }
      else {
        model <- lm("Y ~ X", data = dataset)
      }
      if (binary_X) {
        X_model <- glm("X ~ 0", data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm("X ~ 0", data = X_dataset)
      }
      
      vars_selected <- names(model$coefficients)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
    }
    
    else if (method == "stepwise") {
      if (binary_Y) {
        stepwise_model <- step(object    = glm("Y ~ .", data = dataset, family = "binomial"), # all variable base
                               direction = "both",                                            # stepwise, not fwd or bwd
                               scope     = list(upper = "Y ~ .", lower = "Y ~ X"),            # exposure X always included
                               trace     = 0)                                                 # suppress output
      }
      else {
        stepwise_model <- step(object    = lm("Y ~ .", data = dataset),            # all variable base
                               direction = "both",                                 # stepwise, not fwd or bwd
                               scope     = list(upper = "Y ~ .", lower = "Y ~ X"), # exposure X always included
                               trace     = 0)
      }
      
      vars_selected <- names(stepwise_model$coefficients)
      vars_selected <- union(c('X'), vars_selected) # always select X
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      if (binary_Y) {
        model <- glm(model_formula, data = dataset, family = "binomial")
      }
      else {
        model <- lm(model_formula,  data = dataset)
      }
      if (binary_X) {
        X_model <- glm(X_model_formula, data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm(X_model_formula,  data = X_dataset)
      }
    }
    
    else if (method == "stepwise_X") {
      if (binary_X) {
        stepwise_X_model <- step(object    = glm("X ~ .", data = X_dataset, family = "binomial"), # all variable base
                                 direction = "both",                                              # stepwise, not fwd or bwd
                                 scope     = list(upper = "X ~ .", lower = "X ~ 0"),              # constant term
                                 trace     = 0)                                                   # suppress output
      }
      else {
        stepwise_X_model <- step(object    = lm("X ~ .", data = X_dataset),          # all variable base
                                 direction = "both",                                 # stepwise, not fwd or bwd
                                 scope     = list(upper = "X ~ .", lower = "X ~ 0"), # constant term
                                 trace     = 0)
      }
      
      vars_selected <- names(stepwise_X_model$coefficients)
      vars_selected <- union(c('X'), vars_selected) # always select X
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      if (binary_Y) {
        model <- glm(model_formula, data = dataset, family = "binomial")
      }
      else {
        model <- lm(model_formula,  data = dataset)
      }
      if (binary_X) {
        X_model <- glm(X_model_formula, data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm(X_model_formula,  data = X_dataset)
      }
    }
    
    else if (method == "two_step_lasso") {
      if (binary_Y) {
        cv_lasso_model <- cv.glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), family = 'binomial', alpha=1)
        lambda         <- cv_lasso_model$lambda.min
        lasso_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), family = 'binomial', alpha=1, lambda=lambda)
      }
      else {
        cv_lasso_model <- cv.glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1)
        lambda         <- cv_lasso_model$lambda.min
        lasso_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1, lambda=lambda)
      }
      
      lasso_coefs        <- as.vector(lasso_model$beta)
      names(lasso_coefs) <- rownames(lasso_model$beta)
      
      vars_selected <- names(lasso_coefs[lasso_coefs != 0.0])
      vars_selected <- union(c('X'), vars_selected) # always select X
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      if (binary_Y) {
        model <- glm(model_formula, data = dataset, family = "binomial")
      }
      else {
        model <- lm(model_formula,  data = dataset)
      }
      if (binary_X) {
        X_model <- glm(X_model_formula, data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm(X_model_formula,  data = X_dataset)
      }
    }
    
    else if (method == "two_step_lasso_X") {
      if (binary_X) {
        cv_lasso_X_model <- cv.glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), family = 'binomial', alpha=1)
        lambda           <- cv_lasso_X_model$lambda.min
        lasso_model      <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), family = 'binomial', alpha=1, lambda=lambda)
      }
      else {
        cv_lasso_X_model <- cv.glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1)
        lambda           <- cv_lasso_X_model$lambda.min
        lasso_model      <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1, lambda=lambda)
      }
      
      lasso_coefs        <- as.vector(lasso_model$beta)
      names(lasso_coefs) <- rownames(lasso_model$beta)
      
      vars_selected <- names(lasso_coefs[lasso_coefs != 0.0])
      vars_selected <- union(c('X'), vars_selected) # always select X
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      if (binary_Y) {
        model <- glm(model_formula, data = dataset, family = "binomial")
      }
      else {
        model <- lm(model_formula,  data = dataset)
      }
      if (binary_X) {
        X_model <- glm(X_model_formula, data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm(X_model_formula,  data = X_dataset)
      }
    }
    
    else if (method == "two_step_lasso_union") {
      # X model
      if (binary_X) {
        cv_lasso_X_model <- cv.glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), family = 'binomial', alpha=1)
        lambda_X         <- cv_lasso_X_model$lambda.min
        lasso_X_model    <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), family = 'binomial', alpha=1, lambda=lambda_X)
      }
      else {
        cv_lasso_X_model <- cv.glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1)
        lambda_X         <- cv_lasso_X_model$lambda.min
        lasso_X_model    <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1, lambda=lambda_X)
      }
      
      lasso_X_coefs        <- as.vector(lasso_X_model$beta)
      names(lasso_X_coefs) <- rownames(lasso_X_model$beta)
      
      X_vars_selected <- names(lasso_X_coefs[lasso_X_coefs != 0.0])
      X_vars_selected <- union(c('X'), X_vars_selected) # always select X
      X_vars_selected <- X_vars_selected[X_vars_selected != "(Intercept)"]
      
      
      # Y model
      if (binary_Y) {
        cv_lasso_Y_model <- cv.glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), family = 'binomial', alpha=1)
        lambda_Y         <- cv_lasso_Y_model$lambda.min
        lasso_Y_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), family = 'binomial', alpha=1, lambda=lambda_Y)
      }
      else {
        cv_lasso_Y_model <- cv.glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1)
        lambda_Y         <- cv_lasso_Y_model$lambda.min
        lasso_Y_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1, lambda=lambda_Y)
      }
      
      lasso_Y_coefs        <- as.vector(lasso_Y_model$beta)
      names(lasso_Y_coefs) <- rownames(lasso_Y_model$beta)
      
      Y_vars_selected <- names(lasso_Y_coefs[lasso_Y_coefs != 0.0])
      Y_vars_selected <- union(c('X'), Y_vars_selected) # always select X
      Y_vars_selected <- Y_vars_selected[Y_vars_selected != "(Intercept)"]
      
      
      # union model
      vars_selected   <- union(X_vars_selected, Y_vars_selected)
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
      
      if (binary_Y) {
        model <- glm(model_formula, data = dataset, family = "binomial")
      }
      else {
        model <- lm(model_formula,  data = dataset)
      }
      if (binary_X) {
        X_model <- glm(X_model_formula, data = X_dataset, family = "binomial")
      }
      else {
        X_model <- lm(X_model_formula,  data = X_dataset)
      }
    }
    
    # record model coefficients
    current_coefs <- fill_in_blanks(model$coefficients, var_names_except_Y_with_intercept)
    model_coefs[ method, , repetition] <- current_coefs
    
    # record covariate selection
    current_cov_selection        <- rep(0, times = length(var_names_except_Y))
    names(current_cov_selection) <- var_names_except_Y
    for (var in vars_selected) {
      current_cov_selection[var] <- 1
    }
    cov_selection[ method, , repetition] <- current_cov_selection
    
    
    # record results metrics
    results[ method, "pred_mse", repetition] <- mean(model$residuals^2)
    results[ method, "model_SE", repetition] <- (coef(summary(model))[, "Std. Error"])['X']
    results[ method, "emp_SE", repetition]   <- NaN # filled-in after
    
    results[ method, "r_squared_X", repetition] <- ifelse( binary_X, NaN, r_squared_X(X_model = X_model, X_test_data = X_dataset) )
    results[ method, "r_squared_Y", repetition] <- ifelse( binary_Y, NaN, r_squared_Y(model = model, test_data = dataset) )
    
    results[ method, "causal_true_value", repetition] <- causal
    results[ method, "causal_estimate", repetition]   <- current_coefs['X']
    results[ method, "causal_bias", repetition]       <- NaN # filled-in after
    
    within_CI <- 0.0
    CI        <- confint(model, 'X', level = 0.95)
    if ((causal > CI[1]) && (causal < CI[2])) {
      within_CI <- 1.0
    }
    results[ method, "causal_coverage", repetition]   <- within_CI
    
    results[ method, "open_paths", repetition]    <- num_total_conf
    results[ method, "blocked_paths", repetition] <- length(vars_selected[vars_selected != "X"])
  }
}

# Take mean across repetitions
final_results <- as.data.frame(apply(results, c(1,2), mean))

# fill-in other results
for (method in model_methods) {
  causal_effect_estimates               <- c(results[method, "causal_estimate", ])
  final_results[ method, "emp_SE"]      <- sd(causal_effect_estimates)
  final_results[ method, "causal_bias"] <- mean(causal_effect_estimates - causal)
}

# Round to 3 digits
final_results <- round_df(final_results, digits=3)

# Process coefficients (NB: we omit NaNs here for interpretability)
final_model_coefs <- as.data.frame(apply(model_coefs, c(1,2), function(x) mean(na.omit(x))))
final_model_coefs <- round_df(final_model_coefs, digits=3)

# Process cov selection
final_cov_selection <- as.data.frame(apply(cov_selection, c(1,2), mean))
final_cov_selection <- round_df(final_cov_selection, digits=3)



# ----- Present and save results -----

# Prevalences of binary variables
message("\n\nMeans of X, Y, Zi, equal to prevalence if binary")

message("\nX:")
print(binary_X_prevalence)
print(summary(dataset$X))
print(var(dataset$X))

message("\nY:")
print(binary_Y_prevalence)
print(summary(dataset$Y))
print(var(dataset$Y))

message("\nZ1:")
print(binary_Z_prevalence)
print(summary(dataset$Z1))
print(var(dataset$Z1))

print(head(dataset))
stop("dev")

# Coefficients fitted and error-variance fitted
message("\n\nTrue Coefficients of DAG and Variance of error of Y")
print("Coefficients of Z on X:")
print(beta_Xs)
print("Coefficients of Z on Y:")
print(beta_Ys)
print("Coefficient of X on Y:")
print(causal)
print(var_Y)
print(var_error_Y)
print(var_Y + var_error_Y)

message("\n\nObserved Coefficients")
print(final_model_coefs)

message("\n\nObserved Covariate Selection")
print(final_cov_selection)

# Covariance Matrices
analytic_cov_matrix <- determine_subgroup_cov_matrix(num_total_conf = num_total_conf,
                                                     var_names      = var_names,
                                                     beta_X         = beta_X,
                                                     causal         = causal,
                                                     Z_correlation  = Z_correlation,
                                                     target_r_sq_X  = target_r_sq_X,
                                                     target_r_sq_Y  = target_r_sq_Y)

message("\n\nNon-subgroup non-binary Analytic Covariance:")
print(analytic_cov_matrix)

observed_cov_matrix <- round_df(as.data.frame(cov(dataset)), digits=3)
message("\n\nObserved Covariance:")
print(observed_cov_matrix)

message("\n\nError Results:")
print(final_results[, c(1:3)])

message("\n\nObserved (Continuous) R2 Values:")
print(final_results[, c(4:5)])

message("\n\nCausal Effect Estimation:")
print(final_results[, c(6:9)])

message("\n\nBlocking Open Backdoor Pathways:")
print(final_results[, c(10:11)])


# Save to file
id_string <- paste("sim_run_", n_simulation, "_scenario_", n_scenario, sep='')
write.csv(final_results,       paste("../data/", id_string, "_results.csv", sep=''))
write.csv(final_model_coefs,   paste("../data/", id_string, "_model_coefs.csv", sep=''))
write.csv(final_cov_selection, paste("../data/", id_string, "_cov_selection.csv", sep=''))

write.csv(as.data.frame(analytic_cov_matrix), paste("../data/", id_string, "_analytic_cov_matrix.csv", sep=''))
write.csv(as.data.frame(observed_cov_matrix), paste("../data/", id_string, "_observed_cov_matrix.csv", sep=''))


