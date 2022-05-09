# A package for randomly sampling dirichlet random variables. 
# install.packages('DirichletReg')
library(DirichletReg)


# Set up simulation -------------------------------------------------------

# Get a distance matrix from phyloseq package 
data("GlobalPatterns")
#devtools::install_github("david-barnett/microViz@0.9.1")
library(microViz)
GlobalPatterns
GP = ps_filter(GlobalPatterns, SampleType == "Feces")
gp1  = transform_sample_counts(GP, function(x) x / sum(x) )
filtered <-  filter_taxa(gp1, function(x) mean(x) > 1e-4, TRUE)
filtered

library(adephylo)
D <- as.matrix(distTips(phy_tree(filtered), method = "patristic"))
dim(D)
set.seed(506)
i <- sample(1:329, 15)

D_sim <- D[i,i]
dim(D_sim)

# set parameters
omega <- 1
rho <- 30

n <- 1000
p <- 15
q <- 1

# include intercept? 
# X is for each sample
X <- c(rep(0, 400), rep(1,600))
# beta is for each taxa 
beta <- c(rep(-1, 5), rep(2,10))


# source helper functions 
source(here::here("R","dm_cor_gee_clean.R"))

eta <- get_eta(X, beta , n, p)
alpha <- exp(eta)
alpha
summary(alpha)
alpha0 <- rowSums(matrix(alpha, nrow = p, byrow = T))

y_mean <- alpha/alpha0  


library(Matrix)
# dirichlet correlation (block diagonal of different blocks for each sample)
cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*D_sim))
R <- bdiag(Rs)

# 1/sqrt(Var_D)
A <- Diagonal(n*p)
diag(A) <- sqrt(var_dirichlet(alpha,n,p))
# sqrt(Var_D) * R * (sqrt(Var_D))
V <- A %*% R %*% A


# Generate 1 sample with the mean and variance structure 
library(MASS)
ys <- mvrnorm(1, y_mean, V )
ys


# Try dirichlet now. 
y_i <- list()
mat <- matrix(alpha, nrow = n, ncol = p, byrow = T)
for(i in 1:1000){
  alphai <- mat[i,]
  y_i[[i]] <- rdirichlet(1, alphai)
}

ys <- unlist(y_i)

# Standardizing: 
# How much standardizing is needed? 



# Set any negative ys to 0. 
mean(ys < 0)
ys_pos <- ifelse(ys < 0, 0, ys)


# Do col sums need to add to 1? 
# What do colsums look like? 
rowSums(matrix(ys, nrow = n, ncol =p, byrow = T))



# Shouldnt these be close to 1?  
# Standardize: 
ys_one <- data.frame(id = rep(1:100, each = 15), 
                     ys_pos) %>%
  group_by(id) %>%
  mutate(ysum = sum(ys_pos)) %>% 
  mutate(ys_one = ys_pos/ ysum) %>% 
  pull(ys_one)
ys_one

rowSums(matrix(ys_one, nrow = n, ncol = p, byrow = T))

# Check if y has the same moments now? how to do this
# since cant get alpha from y? right? 






# Run model ---------------------------------------------------------------


# Try simulation! 
source(here::here("R","dm_cor_gee_clean.R"))
model_output <- dm_cor_gee(Y = ys, X = X,
                           sample_id = rep(1:n, each = 15), 
                           ASV_id = rep(letters[1:15], n),
                           distance_matrix = D_sim, 
                           max_iter = 20, tol = .00001, 
                           gamma = 1, lambda = .0001, 
                           save_beta = T, intercept = T)

model_output$betas
descr <- "Simulation, only dirichlet"



# Diagnostics -------------------------------------------------------------


plots <- fun_diagnostics(model_output, descr)
plots
