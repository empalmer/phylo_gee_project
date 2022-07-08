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
# This is one x for each sample

X <- c(rep(0, .4 * n), rep(1, .6 * n))

X <- model.matrix(~x, data.frame(x = X))
X <- as.list(as.data.frame(t(X)))
X <- map(X, ~kronecker(diag(p),t(.x))) %>% 
  reduce(rbind)
X


# Set betas (need one for each p)
#beta <- c(rep(-1, 5), rep(2, 10))
# All positive 
beta <- rep(3, 15)


# Calculate alpha ---------------------------------------------------------
source(here::here("R", "helpers.R"))
source(here::here("R", "dirichlet_functions.R"))
source(here::here("R", "gee_functions.R"))

# True eta has no intercept included
eta <- as.vector(kronecker(X, diag(p)) %*% beta)
# log(alpha) = X*beta
# alpha = e^eta
alpha <- exp(eta)
alpha
# Calculate the "true" variance and correlation
alpha0 <- rowSums(matrix(alpha, nrow = n, ncol = p, byrow = T))
alpha0



true_var <- get_dirichlet_var(alpha, n, p)
true_cor <- get_dirichlet_cor(alpha, n, p)
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



# how close is the sample variance to true variance?
var(ys[1:15])
true_var[1]


# True GEE eqn values -----------------------------------------------------

X <- model.matrix(~x, data.frame(x = X))
X <- as.list(as.data.frame(t(X)))
map(X, ~kronecker(diag(2),t(.x))) %>% 
  reduce(rbind)

ymean <- alpha / rep(alpha0, each =15)
partials <- calculate_partials(alpha, alpha0, n, p, X)


true_resid <- diag(true_A %*% Diagonal(x = ys - ymean))
true_phi <- 1 / (as.numeric(sum(true_resid^2) * (1 / (n*p - (p*q - 1)))))

V_inv <- true_phi * true_A %*% true_R_inv %*% true_A

true_G <- partials %*% V_inv %*% as.matrix(ys - ymean)
sum(as.vector(true_G))
true_H <- -partials %*% V_inv %*% t(partials)

true_update <- Matrix::solve(true_H, true_G)
true_update

# Calculate new beta value based on
# beta+ = beta + gamma H^-1 G
beta_new <- beta -  as.vector(true_update)
beta_new

beta <- beta_new

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


# Stability/function simulation -------------------------------------------------------------
source(here::here("R", "helpers.R"))
source(here::here("R", "dirichlet_functions.R"))
source(here::here("R", "gee_functions.R"))

n <- 5000
p <- 15
q <- 2
lambda <- 10000

# Set betas (need one for each p)
# beta <- rep(3, 15)
beta <- c(rep(-4, 5), rep(2, 10))
# add true 0 intercept
beta <- as.vector(t(matrix(c( rep(0, p),beta), nrow = p)))
# If intercept: 
x <- c(rep(0, .4 * n), rep(1, .6 * n))
# change to 1 if including intercept
x <- model.matrix(~ 1 +  x)
X <- as.list(as.data.frame(t(x)))
X <- map(X, ~kronecker(diag(p),t(.x))) %>% 
  reduce(rbind)
# All positive 
#beta <- rep(3, 15)
# True eta has no intercept included
eta <- as.vector( X  %*% beta)
# log(alpha) = X*beta
# alpha = e^eta
alpha <- exp(eta)
# Calculate the "true" variance and correlation
alpha0 <- rowSums(matrix(alpha, nrow = n, ncol = p, byrow = T))
ys <- simulate_dirichlet_y(alpha, n, p)

# if initilizing beta as 0 to make algorithm figure it out
beta <- rep(0, p*q)

update_list <- list()
g_list <- list()
beta_list <- list()
diff_list <- list()
phi_list <- list()
n_rep <- 20
for(i in 1:n_rep){
  beta_list[[i]] <- beta
  eta <- as.vector(X %*% beta)
  alpha <- exp(eta)
  
  # Calculate the "true" variance and correlation
  alpha0 <- rowSums(matrix(alpha, nrow = n, ncol = p, byrow = T))
  
  var <- get_dirichlet_var(alpha, n, p)
  cor <- get_dirichlet_cor(alpha, n, p)
  A <- Diagonal(n*p)
  diag(A) <- sqrt(1/get_dirichlet_var(alpha,n,p))
  R_inv <- get_R_inv(alpha, omega = 0 , rho = 1, D = diag(n), n, p)
  
  ymean <- alpha / rep(alpha0, each = p)
  partials <- calculate_partials(alpha, alpha0, n, p, x)
  
  resid <- diag(A %*% Diagonal(x = ys - ymean))
  
  phi <- 1 / (as.numeric(sum(resid^2) * (1 / (n*p - (p*q - 1)))))
  phi_list[[i]] <- phi
  V_inv <- phi * A %*% R_inv %*% A
  
  G <- partials %*% V_inv %*% as.matrix(ys - ymean) +
        2 * lambda * tcrossprod(rep(1, p*q)) %*% beta
  #  2 * lambda * p * q *sum(beta)
  g_list[[i]] <- sum(as.vector(G))
  H <- -partials %*% V_inv %*% t(partials) +
        2 * lambda * tcrossprod(rep(1, p*q)) 
  #  2 * lambda * p * q 
  
  update <- Matrix::solve(H, G)
  update_list[[i]] <- update
  
  # Calculate new beta value based on
  # beta+ = beta + gamma H^-1 G
  beta_new <- beta -  as.vector(update)
  
  
  diff_list[[i]] <- sum(abs(beta_new - beta))
  
  beta <- beta_new
  print(i)
}

# Plot beta values across iterations. Are they close to true values? 
# plot changes of betas between iterations
num_beta <- p*q
betas <- as.data.frame(beta_list) 
colnames(betas) <- 1:n_rep
betas %>% 
  mutate(type = rep(c('x'), num_beta), 
         id = factor(1:(num_beta)), 
         type = factor(rep(c("int","x"), p))) %>% 
  pivot_longer(1:n_rep, names_to = "iteration") %>% 
  ggplot(aes(x = as.numeric(iteration), y = value, group = id, linetype = type)) +
  geom_line() + 
  labs(x = "iteration", y = "beta") + 
  ggtitle("beta -4,2, lambda = 10000, add penalty, init 0, n=5000")

beta_list[[n_rep]]



# plot GEE values, stable? Close to 0? 
data.frame(iter = 1:n_rep, 
           g = unlist(g_list)) %>%
  ggplot(aes(x = iter, y = g)) +
  geom_line()

# plot differences of summed abs diff 
data.frame(iter = 1:n_rep, 
           diff = unlist(diff_list)) %>%
  ggplot(aes(x = iter, y = diff)) +
  geom_line()

# plot phi, around 1? 
data.frame(iter = 1:n_rep, 
           phi = unlist(phi_list)) %>%
  ggplot(aes(x = iter, y = phi)) +
  geom_line()


