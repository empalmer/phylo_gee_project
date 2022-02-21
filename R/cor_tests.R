
# testing dirichlet correlation 

alpha <- c(1,2,3)


alpha0 <- sum(alpha)

num <- - sqrt(alpha %*% t(alpha)) 
denom <- sqrt((alpha0-alpha) %*% t((alpha0 -alpha)))
t_cor <- num/denom

solve(t_cor)

A <- matrix(0,nrow = 3, ncol = 3)
diag(A) <- sqrt(1/var_dirichlet(alpha,1,3))
get_dirichlet_cor(alpha, 1, 3)

V <- A %*% t_cor %*% A
solve(V)
