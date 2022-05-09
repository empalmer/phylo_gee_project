# A package for randomly sampling dirichlet random variables. 
# install.packages('DirichletReg')
library(DirichletReg)


n <- 1000
p <- 15
q <- 1

# include intercept? no, should be default to make in model matrix. 
# X is for each sample
X <- c(rep(0, .4*n), rep(1, .6*n))
# beta is for each taxa
beta <- c(rep(-1, 5), rep(2, 10))



# Generate Dirichlet RV ---------------------------------------------------

# Simulate from Dirichlet distribution given alpha, n, p
# Need to draw n different samples, since each sample is dirichlet distributed
# Then combine alphas to one vector. 

y_i <- list()
mat <- matrix(alpha, nrow = n, ncol = p, byrow = T)
for (i in 1:1000) {
  alphai <- mat[i, ]
  y_i[[i]] <- rdirichlet(1, alphai)
}

ys <- unlist(y_i)