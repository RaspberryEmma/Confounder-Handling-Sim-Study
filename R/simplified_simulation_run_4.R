# ****************************************
# Confounder Handling Simulation Study
#
# Simplified, less flexible refactor of all existing sim code
#
# Emma Tarmey
#
# Started:          06/02/2025
# Most Recent Edit: 13/02/2025
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



# ----- Random Number Generation -----

# fix RNG seed based on current run and scenario
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>%
  filter(simulation_run      == 1) %>%
  filter(simulation_scenario == 1)
set.seed(seed$seed)



# ----- Parameters ------

# Run
n_simulation      <- 4 # see Table!

n_obs             <- 10000
n_rep             <- 2000
Z_correlation     <- 0.2
Z_subgroups       <- 2.0
target_r_sq_X     <- 0.4
target_r_sq_Y     <- 0.4
causal            <- 0.5

binary_X          <- FALSE
binary_Y          <- FALSE
binary_Z          <- FALSE

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



# ----- Helper Functions -----


# Round all numeric columns in a given data frame to a given number of digits
# Credit to: https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
round_df <- function(df, digits) {
  nums      <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return (df)
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


# Estrimate the variance in the error term for Y
# Guarantees our value of R2Y is respected for arbitrary value of causal effect
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


# Estimate covariance matrix of entire system
# Allows us to check in-one-go whether covariates are being generated properly
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


# Estimate covariance matrix of entire system
# Structurally same as above but with formulae swapped out
# Specifically this function gives the more-complicated "subgroups" formulae for variances
determine_subgroup_cov_matrix <- function(num_total_conf = NULL,
                                          var_names      = NULL,
                                          beta_Xs        = NULL,
                                          beta_Ys        = NULL,
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
  
  var_error_Y <- determine_subgroup_var_error_Y(num_total_conf = num_total_conf,
                                                beta_Xs        = beta_Xs,
                                                beta_Ys        = beta_Ys,
                                                causal         = causal,
                                                Z_correlation  = Z_correlation,
                                                target_r_sq_Y  = target_r_sq_Y)
  
  # variances
  for (i in 1:num_vars) {
    # variances here
    if (var_names[i] == 'X') {
      analytic_cov[i, i] <- NaN
    }
    else if (var_names[i] == 'Y') {
      analytic_cov[i, i] <- NaN
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
beta_Y  <- beta_X # Z on Y

# variance of the error term for Y
var_error_Y <- determine_var_error_Y(num_total_conf = num_total_conf,
                                     beta_X         = beta_X,
                                     causal         = causal,
                                     Z_correlation  = Z_correlation,
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
    X  <- X + (beta_X * Z1) + (beta_X * Z2) + (beta_X * Z3) + (beta_X * Z4)
    Y  <- Y + (beta_Y * Z1) + (beta_Y * Z2) + (beta_Y * Z3) + (beta_Y * Z4)
    
    # record Zs
    # index formula places Zs such that all Zs of one sub-group are next to each other
    dataset[, (2 + (0 * num_of_batches) + i)] <- Z1
    dataset[, (2 + (1 * num_of_batches) + i)] <- Z2
    dataset[, (2 + (2 * num_of_batches) + i)] <- Z3
    dataset[, (2 + (3 * num_of_batches) + i)] <- Z4
    
    # TESTING
    # print( (2 + (0 * num_of_batches) + i) )
    # print( (2 + (1 * num_of_batches) + i) )
    # print( (2 + (2 * num_of_batches) + i) )
    # print( (2 + (3 * num_of_batches) + i) )
  }
  
  # add error terms to X and Y
  error_X <- rnorm(n = n_obs, mean = 0, sd = 1)
  error_Y <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
  X <- X + error_X
  Y <- Y + error_Y
  
  # add causal effect (X on Y)
  Y <- Y + (causal * X)
  
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
                   "two_step_lasso", "two_step_lasso_X")

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
      model   <- lm("Y ~ .", data = dataset)
      X_model <- lm("X ~ .", data = X_dataset)
      
      vars_selected <- names(model$coefficients)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
    }
    
    else if (method == "linear_unadjusted") {
      model   <- lm("Y ~ X", data = dataset)
      X_model <- lm("X ~ 0", data = X_dataset)
      
      vars_selected <- names(model$coefficients)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
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
      stepwise_X_model   <- step(object    = lm("X ~ .", data = X_dataset),          # all variable base
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
    
    results[ method, "r_squared_X", repetition] <- r_squared_X(X_model = X_model, X_test_data = X_dataset)
    results[ method, "r_squared_Y", repetition] <- r_squared_Y(model = model, test_data = dataset)
    
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
  causal_effect_estimates            <- c(results[method, "causal_estimate", ])
  
  final_results[ method, "emp_SE"]      <- sd(causal_effect_estimates)
  final_results[ method, "causal_bias"] <- mean(causal_effect_estimates - causal)
}

# Round to 3 digits
final_results <- round_df(final_results, digits=4)

# Process coefficients
final_model_coefs <- as.data.frame(apply(model_coefs, c(1,2), mean))
final_model_coefs <- round_df(final_model_coefs, digits=4)

# Process cov selection
final_cov_selection <- as.data.frame(apply(cov_selection, c(1,2), mean))
final_cov_selection <- round_df(final_cov_selection, digits=4)



# ----- Present and save results -----

 # Coefficients fitted and error-variance fitted
coefs <- c(beta_X, beta_Y, causal)
names(coefs) <- c("beta_X", "beta_Y", "causal")

message("\n\nTrue Coefficients of DAG and Variance of error of Y")
print(coefs)
print(var_error_Y)

message("\n\nObserved Coefficients")
print(final_model_coefs)

message("\n\nObserved Covariate Selection")
print(final_cov_selection)

# Covariance Matrices
analytic_cov_matrix <- determine_cov_matrix(num_total_conf = num_total_conf,
                                            var_names      = var_names,
                                            beta_X         = beta_X,
                                            beta_Y         = beta_Y,
                                            causal         = causal,
                                            Z_correlation  = Z_correlation,
                                            target_r_sq_X  = target_r_sq_X,
                                            target_r_sq_Y  = target_r_sq_Y)

# analytic_subgroup_cov_matrix <- determine_subgroup_cov_matrix(num_total_conf = num_total_conf,
#                                                               var_names      = var_names,
#                                                               beta_Xs        = c(beta_X, beta_X, beta_X, beta_X),
#                                                               beta_Ys        = c(beta_X, beta_X, beta_X, beta_X),
#                                                               causal         = causal,
#                                                               Z_correlation  = Z_correlation,
#                                                               target_r_sq_X  = target_r_sq_X,
#                                                               target_r_sq_Y  = target_r_sq_Y)

message("\n\nNon-subgroup Analytic Covariance:")
print(analytic_cov_matrix)

# message("\n\n(TBC) Subgroup Analytic Covariance:")
# print(analytic_subgroup_cov_matrix)

observed_cov_matrix <- round_df(as.data.frame(cov(dataset)), digits=3)
message("\n\nObserved Covariance:")
print(observed_cov_matrix)

message("\n\nError Results:")
print(final_results[, c(1:3)])

message("\n\nOberserved R2 Values:")
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


