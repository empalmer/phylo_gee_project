# A package for randomly sampling dirichlet random variables.
# install.packages('DirichletReg')
# needed for rdirichlet() fxn
library(DirichletReg)


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


# Run model ---------------------------------------------------------------

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
  gamma = 1,
  lambda = .0001,
  save_beta = T,
  intercept = T,
  only_dir_cor = T, 
  fixed = T,
  A = true_A, 
  R_inv = true_R_inv
)

model_output$betas
old_betas_plus <- model_output$betas
descr <- "Simulation, only dirichlet"


# Diagnostics -------------------------------------------------------------

source(here::here("R", "diagnostics_plots.R"))
plots <- fun_diagnostics(model_output, descr)
plots
