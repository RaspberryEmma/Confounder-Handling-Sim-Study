# ****************************************
# Confounder Handling Simulation Study
#
# Simplified, less flexible refactor of all existing sim code
#
# Emma Tarmey
#
# Started:          06/02/2025
# Most Recent Edit: 04/06/2025
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
using("dplyr", "glmnet", "speedglm", "tidyr")

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}




# ----- Parameters ------

# Run
n_simulation      <- 4 # see Table!

n_obs             <- 10000
n_rep             <- 2000
Z_correlation     <- 0.2
Z_subgroups       <- 4.0
target_r_sq_X     <- 0.6
target_r_sq_Y     <- 0.4
causal            <- 0.15

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
  include_step    <- as.numeric(args[5])
  
} else {
  n_scenario      <- 1
  num_total_conf  <- 32
  num_meas_conf   <- 32
  num_unmeas_conf <- 0
  include_step    <- 0
}



# ----- Random Number Generation -----

# fix RNG seed based on current run and scenario
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>%
  filter(simulation_run      == n_simulation) %>%
  filter(simulation_scenario == n_scenario)
set.seed(seed$seed)



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
  # first item does not require a '+' sign
  formula_string <- paste("Y ~ ", vars_selected[1], sep = "")
  
  for (var in vars_selected[-1]) {
    formula_string <- paste(formula_string, " + ", var, sep = "")
  }
  return(formula_string)
}


# Create a string corresponding to a regression of X on a given set of covariates
make_X_model_formula <- function(vars_selected = NULL) {
  # remove X
  if ('X' %in% vars_selected) {
    vars_selected <- vars_selected[! vars_selected%in% c('X')]
  }
  
  # first item does not require a '+' sign
  formula_string <- paste("X ~ ", vars_selected[1], sep = "")
  
  for (var in vars_selected[-1]) {
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



# Determine the mean squared error of the model on a new dataset
# Generates new test set internally, in principle works for all basic model types
# lm, glm, speedglm etc
glm_mse <- function(model = NULL) {
  test_data        <- generate_dataset()
  predicted_values <- predict(model, type = "response", newdata = test_data)
  return( mean((test_data$Y - predicted_values)^2) )
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

if (include_step == 1) {
  model_methods <- c("stepwise")
} else if (include_step == 2) {
  model_methods <- c("stepwise_X")
} else {
  model_methods <- c("linear", "linear_unadjusted", "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")
}

results_methods <- c("pred_mse", "model_SE", "emp_SE",
                     "r_squared_X", "r_squared_Y",
                     "causal_true_value", "causal_estimate", "causal_bias", "causal_coverage",
                     "open_paths", "blocked_paths", "convergence_rate")

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


# track progress through simulation
previous_time <- as.numeric(Sys.time())*1000
current_time  <- as.numeric(Sys.time())*1000
delta_time    <- current_time - previous_time
remaining_rep <- n_rep + 1

# run simulation repetitions
for (repetition in 1:n_rep) {
  
  remaining_rep   <- remaining_rep - 1
  previous_time   <- current_time
  current_time    <- as.numeric(Sys.time())*1000
  delta_time      <- current_time - previous_time
  remaining_time  <- remaining_rep * delta_time
  expected_finish <- signif(((current_time + remaining_time)/1000), digits = 13) %>% as.POSIXct()
  
  # track progress
  message( paste("\n\nRunning Simulation Run ", n_simulation,
                 ", Scenario ", n_scenario,
                 ", Iteration ", repetition, "/", n_rep,
                 "\nLast iteration computed in ", signif(delta_time/1000, 4), " seconds",
                 "\nExpected finish time ", expected_finish, " for current confounder set",
                 "\n", sep = "") )
  
  # generate data
  dataset   <- generate_dataset()
  
  # cut-up versions of the data as needed
  X_dataset <- subset(dataset, select=-c(Y))
  Z_dataset <- subset(dataset, select=-c(Y, X))
  Y_column  <- subset(dataset, select=c(Y))
  X_column  <- subset(dataset, select=c(X))
  
  for (method in model_methods) {
    message(paste0("Fitting model for method ", method))
    method_start_time <- as.numeric(Sys.time())*1000
    
    # fit model
    model   <- NULL
    X_model <- NULL # required to make R2X well-defined
    
    if (method == "linear") {
      vars_selected   <- colnames(X_dataset) # all
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
    }
    
    else if (method == "linear_unadjusted") {
      vars_selected   <- c("X") # none except X
      model_formula   <- "Y ~ X"
      X_model_formula <- "X ~ 0"
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
    }
    
    else if (method == "two_step_lasso") {
      cv_lasso_model <- cv.glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1)
      lambda         <- cv_lasso_model$lambda.min
      lasso_model    <- glmnet(x = data.matrix(X_dataset), y = data.matrix(Y_column), alpha=1, lambda=lambda)
      
      lasso_coefs        <- as.vector(lasso_model$beta)
      names(lasso_coefs) <- rownames(lasso_model$beta)
      
      vars_selected <- names(lasso_coefs[lasso_coefs != 0.0])
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
    }
    
    else if (method == "two_step_lasso_X") {
      cv_lasso_X_model <- cv.glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1)
      lambda           <- cv_lasso_X_model$lambda.min
      lasso_model      <- glmnet(x = data.matrix(Z_dataset), y = data.matrix(X_column), alpha=1, lambda=lambda)
      
      lasso_coefs        <- as.vector(lasso_model$beta)
      names(lasso_coefs) <- rownames(lasso_model$beta)
      
      vars_selected <- names(lasso_coefs[lasso_coefs != 0.0])
      vars_selected <- c('X', vars_selected)
      vars_selected <- vars_selected[vars_selected != "(Intercept)"]
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
    }
    
    else if (method == "two_step_lasso_union") {
      # extract covariate selections from prev models (named binary vectors)
      X_vars_selected <- cov_selection[ "two_step_lasso_X", , repetition]
      Y_vars_selected <- cov_selection[   "two_step_lasso", , repetition]
      
      # convert to vectors of variable names
      X_vars_selected <- names(X_vars_selected[X_vars_selected == 1])
      Y_vars_selected <- names(Y_vars_selected[Y_vars_selected == 1])
      
      # union
      vars_selected <- union(X_vars_selected, Y_vars_selected)
      
      model_formula   <- make_model_formula(vars_selected = vars_selected)
      X_model_formula <- make_X_model_formula(vars_selected = vars_selected)
    }
    
    # fit models
    model   <- lm(model_formula,   data = dataset)
    X_model <- lm(X_model_formula, data = X_dataset)
    
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
    results[ method, "pred_mse", repetition] <- glm_mse(model = model)
    if (binary_Y) {
      # logistic speedglm object
      results[ method, "model_SE", repetition] <- (coef(summary(model)))['X', 'Std. Error']
    }
    else {
      # linear lm object
      results[ method, "model_SE", repetition] <- (coef(summary(model))[, "Std. Error"])['X']
    }
    results[ method, "emp_SE", repetition]   <- NaN # filled-in after
    
    results[ method, "r_squared_X", repetition] <- r_squared_X(X_model = X_model, X_test_data = X_dataset)
    results[ method, "r_squared_Y", repetition] <- r_squared_Y(model = model, test_data = dataset)
    
    results[ method, "causal_true_value", repetition] <- causal
    results[ method, "causal_estimate", repetition]   <- current_coefs['X']
    results[ method, "causal_bias", repetition]       <- NaN # filled-in after
    
    if (binary_Y && (method != "linear_unadjusted")) {
      # speedglm object
      results[ method, "convergence_rate", repetition] <- as.integer(model$convergence)
    }
    else {
      # lm object
      results[ method, "convergence_rate", repetition] <- NaN
    }
    
    within_CI <- 0.0
    CI        <- confint(model, 'X', level = 0.95)
    if ((!is.na(CI[1])) && (!is.na(CI[2]))) {
      if ((causal > CI[1]) && (causal < CI[2])) {
        within_CI <- 1.0
      }
    }
    results[ method, "causal_coverage", repetition]   <- within_CI
    
    results[ method, "open_paths", repetition]    <- num_total_conf
    results[ method, "blocked_paths", repetition] <- length(vars_selected[vars_selected != "X"])
    
    method_end_time <- as.numeric(Sys.time())*1000
    method_elapsed_time <- (method_end_time - method_start_time)/1000
    message(paste0("Method ", method, " ran in ", signif(method_elapsed_time, digits=4), " seconds\n"))
    
  } # methods loop
} # repetitions loop

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

# # Covariance Matrices
analytic_cov_matrix <- determine_subgroup_cov_matrix(num_total_conf = num_total_conf,
                                                     var_names      = var_names,
                                                     beta_X         = beta_X,
                                                     causal         = causal,
                                                     Z_correlation  = Z_correlation,
                                                     target_r_sq_X  = target_r_sq_X,
                                                     target_r_sq_Y  = target_r_sq_Y)

message("\n\nNon-subgroup Analytic Covariance:")
print(analytic_cov_matrix)

observed_cov_matrix <- round_df(as.data.frame(cov(dataset)), digits=3)
message("\n\nObserved Covariance:")
print(observed_cov_matrix)

message("\n\nError Results:")
print(final_results[, c(1:3)])

message("\n\nObserved R2 Values:")
print(final_results[, c(4:5)])

message("\n\nCausal Effect Estimation:")
print(final_results[, c(6:9)])

message("\n\nBlocking Open Backdoor Pathways:")
print(final_results[, c(10:11)])


# Save to file
id_string <- paste("sim_run_", n_simulation, "_scenario_", n_scenario, "_step_", include_step, sep='')
write.csv(dataset,             paste("../data/", id_string, "_dataset.csv", sep=''))
write.csv(final_results,       paste("../data/", id_string, "_results.csv", sep=''))
write.csv(final_model_coefs,   paste("../data/", id_string, "_model_coefs.csv", sep=''))
write.csv(final_cov_selection, paste("../data/", id_string, "_cov_selection.csv", sep=''))

write.csv(as.data.frame(analytic_cov_matrix), paste("../data/", id_string, "_analytic_cov_matrix.csv", sep=''))
write.csv(as.data.frame(observed_cov_matrix), paste("../data/", id_string, "_observed_cov_matrix.csv", sep=''))

