# A package for randomly sampling dirichlet random variables.
# install.packages('DirichletReg')
# needed for rdirichlet() fxn
library(DirichletReg)
library(tidyverse)
library(Matrix)
# Set parameters (n, p, q) ------------------------------------------------

n <- 100
p <- 15
q <- 1


# X is for each sample
# single Binary covariate
# Almost even split of case/control
X <- c(rep(0, .4 * n), rep(1, .6 * n))
# Set betas (need one for each p)
beta <- c(rep(-1, 5), rep(2, 10))


# Calculate alpha ---------------------------------------------------------
source(here::here("R", "helpers.R"))
source(here::here("R", "dirichlet_functions.R"))
source(here::here("R", "gee_functions.R"))

# True eta has no intercept included
eta <- get_eta(X, beta, n, p)
# log(alpha) = Xbeta
alpha <- exp(eta)
alpha
# Calculate the "true" variance and correlation
alpha0 <- rowSums(matrix(alpha, nrow = n, ncol = p, byrow = T))
alpha0



true_var <- get_dirichlet_var(alpha, n, p)
true_var

true_cor <- get_dirichlet_cor(alpha, n, p)
true_cor

true_A <- Diagonal(n*p)
diag(true_A) <- sqrt(1/get_dirichlet_var(alpha,n,p))
true_R_inv <- get_R_inv(alpha, omega = 0 , rho = 1, D = diag(n), n, p)

V_inv <- true_A %*% true_R_inv %*% true_A

# Generate Dirichlet RV ---------------------------------------------------

# Simulate from Dirichlet distribution given alpha, n, p
# Need to draw n different samples, since each sample is dirichlet distributed
# Then combine alphas to one vector.

simulate_dirichlet_y <- function(alpha, n, p, seed = 1225){
  set.seed(seed)
  y_i <- list()
  mat <- matrix(alpha, nrow = n, ncol = p, byrow = T)
  for (i in 1:n) {
    alphai <- as.vector(mat[i, ])
    y_i[[i]] <- rdirichlet(1, alphai)
  }
  ys <- unlist(y_i)
  return(ys)
}

ys <- simulate_dirichlet_y(alpha, n, p)


# True GEE eqn values -----------------------------------------------------
X <- model.matrix(~x, data.frame(x = X))
ymean <- alpha / rep(alpha0, each =15)
partials <- calculate_partials(alpha, alpha0, n, p, X)
true_G <- partials %*% V_inv %*% as.matrix(ys - ymean)
true_H <- -partials %*% V_inv %*% t(partials)

true_update <- Matrix::solve(true_H, true_G)
true_update

# Calculate new beta value based on
# beta+ = beta + gamma H^-1 G
beta_new <- beta - gamma * as.vector(update)


# Run model ---------------------------------------------------------------

true_beta <- c(rep(c(0, -1), 5), rep(c(0,2),10))

# Try simulation!
# Fake distance matrix that is the pxp identity
source(here::here("R", "main_gee_dir_phy.R"))
model_output <- dm_cor_gee(
  Y = ys,
  X = X,
  sample_id = rep(1:n, each = p),
  ASV_id = rep(letters[1:15], n),
  distance_matrix = diag(p),
  max_iter = 20,
  tol = .0001,
  gamma = .8,
  lambda = .0001,
  save_beta = T,
  intercept = T,
  only_dir_cor = T, 
  fixed = T,
  A = true_A, 
  R_inv = true_R_inv
 # init_beta = true_beta 
)

matrix(model_output$betas[[20]], nrow = p, ncol = q+1, byrow = T)

descr <- "Simulation, only dirichlet"



# Diagnostics -------------------------------------------------------------
# currently wont work since not outputting omega values... 
source(here::here("R", "diagnostics_plots.R"))
plots <- fun_diagnostics(model_output, descr)
plots
