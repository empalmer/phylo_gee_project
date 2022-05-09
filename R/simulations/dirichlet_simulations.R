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
X <- c(rep(0, .4*n), rep(1, .6*n))
# Set betas (need one for each p)
beta <- c(rep(-1, 5), rep(2, 10))


# Calculate alpha ---------------------------------------------------------
source(here::here("R", "helpers.R"))
source(here::here("R", "dirichlet_functions.R"))

# True eta has no intercept included 
eta <- get_eta(X, beta, n, p)
# log(alpha) = Xbeta
alpha <- exp(eta)


# Generate Dirichlet RV ---------------------------------------------------

# Simulate from Dirichlet distribution given alpha, n, p
# Need to draw n different samples, since each sample is dirichlet distributed
# Then combine alphas to one vector. 

y_i <- list()
mat <- matrix(alpha, nrow = p, ncol = n, byrow = T)
mat[1,]
for (i in 1:1000) {
  alphai <- mat[i, ]
  y_i[[i]] <- rdirichlet(1, alphai)
}

ys <- unlist(y_i)