#' Phylogeny-aware Dirichlet GEE 
#' 
#' 
#' @param Y should be of length n*p
#' @param X should be dimension n*q
#' @param id should be of length n*p
#' @param distance_matrix An object of type dist 
#'
#' @return
#' @export
#'
#' @examples
dm_cor_gee2 <- function(Y, X, id, distance_matrix, intercept = T){
  require(tidyverse)
  require(Matrix)
  require(MASS)
  

  # Set up indeces 
  # Number of subjects i = 1,...,n
  n <- length(unique(id))
  # Number of ASVs j = 1,...,p
  p <- length(Y)/n
  # Convert into correct model matrix format 
  # Number of covariates: k = 1,...,q
  q <- 1 + intercept
  
  # Save a design matrix form for GEE portion
  X_vecs <- X
  # Add intercept, so a dataframe of the Xs plus a col of 1
  X_intercept <- data.frame(int = rep(1, n), 
                            X = X_vecs)
  # convert from dataframe to matrix
  X_mat <- as.matrix(X_intercept, nrow = n)
  
  # This makes the kronecker product 
  X_big_i <- list()
  for( i in 1:n){
    X_big_i[[i]] <- matrix(
      bdiag(
        rep(list(matrix(X_mat[i,], nrow =1)),p)), nrow =p)
  }
  X <- reduce(X_big_i, rbind)
  
  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 

  beta_matrix <- matrix(0, nrow = q, ncol = p)
  
  # fix this: 
  # The Dirichlet link function is the log
  y_means <- colMeans(matrix(Y, ncol = p))
  # Initialize the beta intercept term as ybar/ sum(ybar)
  # But sum = 1
  beta_matrix[1,] <- log(y_means/sum(y_means))

  # convert matrix to vector
  beta <- as.vector(beta_matrix)
  
  # Initialize the A matrix which is a diagonal matrix of the inverse of the square root of the variance
  # Here it is (nxp)x(nxp) 
  # This is using the Matrix package diagonal function
  A <- Diagonal(n*p)
  # Set up distance matrix
  # Switch to take a dist object as an argument 
  d_jk <- distance_matrix[upper.tri(distance_matrix)]
  
  # Since d_jk doesnt have an i index we need to repeat it for each sample
  d_ijk <- rep(d_jk,n)
  
  
  
  # Set up storage to keep track of parameter values in each iteration 
  eta_list <- list()
  alpha_list <- list()
  omega_list <- list()
  rho_list <- list()
  beta_list <- list()
  beta_diffs <- list()
  
  # Initial values for rho, omega, phi 
  rho <- 1
  omega <- 0.5
  phi <- 1
  
  # Initial value for R, R_inverse 
  eta <- as.vector(X %*% beta) 
  alpha <- exp(eta) 
  alpha0 <- colSums(matrix(alpha, nrow = p))
  cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
  Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))
  # Use Moore-Penrose generalized inverse
  R_invs <- map(Rs, ginv)
  R_inv <- bdiag(R_invs)
  
  # Main loop
  count <- 0
  while(count < 5){
    count <- count + 1
    print(paste0("Iteration: ", count))
    
    browser()
    
    # Update beta loop 
    # Depends on "fixed" values of rho, omega and phi, 
    # Which are used to make R_inv 
    beta <- update_beta(Y, X_mat, beta, R_inv, phi, n_iter = 1)
  
    # Update R inverse by updating omega and rho 
    temp_res <- update_R_phi <- function(Y, X, beta)
    R_inv <- temp_res$R_inv
    phi <- temp_res$phi
    
    
    # Save estimates when things start working 
    # eta_list[[count]] <- eta
    # alpha_list[[count]] <- alpha
    # omega_list[[count]] <- omega
    # rho_list[[count]] <- rho
    # beta_list[[count]] <- beta
    # beta_diffs[[count]] <- diffs
  }
  return(list(etas = list(eta_list), 
              alphas = list(alpha_list), 
              omegas = omega_list, 
              rhos = rho_list, 
              betas = list(beta_list), 
              beta_convergences = beta_diffs,
              num_iter = count))
}




#' Update rho, omega, and R inv 
#'
#' @param X 
#' @param beta 
#' @param Y 
#'
#' @return
#' @export
#'
#' @examples
update_phi_R_inv <- function(X, beta, Y){
  eta <- as.vector(X %*% beta) 
  alpha <- exp(eta) 
  alpha0 <- colSums(matrix(alpha, nrow = p))
  mu <-  alpha / rep(alpha0, each = p)
  diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
  
  ### Now update omega, rho, 
  # We use the residuals as the responses for the NLS regression
  resid <- diag(A %*% Diagonal(x = Y - mu))
  # Setup calculated residuals for nls fxn
  # Need to "flatten" so each residual in the residual matrix 
  # These are squared residuals 
  cross_resids <- data.frame(id,resid) %>% 
    group_by(id) %>% 
    group_map(~tcrossprod(.x$resid))
  
  
  cross_resid_vec <- unlist(map(cross_resids, ~.x[upper.tri(.x)]))
  
  # resids_sparse <- as.matrix(bdiag(resids))
  # flat_resids <- resids_sparse[upper.tri(resids_sparse)]
  # flat_resids <- flat_resids[flat_resids != 0]
  # Overdispersion 
  # This is actually phi inverse? According to GEE paper
  # sum of squared residuals 
  phi <- as.numeric(sum(resid^2)*(1/(n-(q-1))))
  
  print(paste0("phi = ", phi))
  
  browser()
  # Get dirichlet correlation part
  cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
  cor_dirichlet_vec <- unlist(map(cor_dirichlet_list, ~.x[upper.tri(.x)]))
  
  
  # cor_dirichlet <- as.matrix(bdiag(cor_dirichlet_list))
  # flat_cor_dir <- cor_dirichlet[upper.tri(cor_dirichlet)]
  # flat_cor_dir <- flat_cor_dir[flat_cor_dir != 0]
  
  estimates <- nls_optim(cross_resid_vec/ phi, cor_dirichlet_vec, d_ijk)
  
  omega <- estimates[1]
  rho <- estimates[2]
  
  print(paste0("omega  = ", omega, "rho = ", rho))
  
  # use current values of omega and rho to create
  # present iteration of combined working correlation matrix. 
  
  Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))
  
  # invert
  # Use Moore-Penrose generalized inverse
  R_invs <- map(Rs, ginv)
  #R_invs <- map(Rs, solve)
  R_inv <- bdiag(R_invs)
  
  
  return(list(omega = omega, 
              rho = rho, 
              R_inv = R_inv))
}



#' Update beta 
#'
#' @param Y 
#' @param X_mat 
#' @param beta 
#' @param R_inv 
#' @param phi 
#' @param n_iter 
#'
#' @return
#' @export
#'
#' @examples
update_beta <- function(Y, X_mat, beta, R_inv, phi, n_iter = 1){
  beta.new <- beta
  diffs <- numeric(1)
  for(s in 1:n_iter){
    print(paste0("Beta iteration ", s))
    beta.old <- beta.new
    eta <- as.vector(X %*% beta.new) 
    alpha <- exp(eta) 
    alpha0 <- colSums(matrix(alpha, nrow = p))
    mu <-  alpha / rep(alpha0, each = p)
    
    diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
    # Make the partials for each block.
    partiali <- list()
    for(i in 1:n){
      alphai0 <- alpha0[i]
      xi <- X_mat[i,]
      alphai <- alpha[((i-1)*p + 1):(i*p)]
      
      # this xi is the vector xi. 
      partiali[[i]] <- (1/alphai0)^2*kronecker((alphai0*diag(alphai) - 
                                                  alphai %*% t(alphai)),xi)
    }
    # NOT bdiag. Instead a long cbinded matrix of pqxnp
    partials <- purrr::reduce(partiali, cbind)
    
    
    # Save V inv 
    # reminder A is A^{-1/2}, and A^t = A (A is diagonal)
    # V_inv <- 1/phi * A %*% R_inv %*% A
    V_inv <- phi * A %*% R_inv %*% A
    
    hess <- partials %*% V_inv %*% t(partials) + diag(rep(.001, q*p))
    # GEE estimating equations/ gradient 
    esteq <- partials %*% V_inv %*% as.matrix(Y - mu)
    
    update <- solve(hess, esteq)
    beta.new <- beta.new - as.vector(update)
    
    # Save to see speed of convergence
    diffs[s] <- sum((beta.new - beta.old)^2)
    print(paste0("beta = ",beta.new))
  }
  beta <- beta.new
  return(beta)
}


#' Function for calculating Dirichlet variance 
#'
#' @param alpha a vector of n*p alpha values 
#' @param n a n*p vector of what sample each alpha is
#' @param p 
#'
#' @return Matrix of variances for each sample+ASV
#' @export
#'
#' @examples
var_dirichlet <- function(alpha,n,p){
  # TODO fix when allow unequal cluster size
  id <- rep(1:n, each = p)
  data.frame(alpha, id) %>% 
    group_by(id) %>% 
    mutate(alpha0 = sum(alpha)) %>%
    ungroup() %>% 
    mutate(v = (alpha*(alpha0 - alpha))/(alpha0^2*(alpha + 1))) %>% 
    pull(v)
}


#' Calculate the formula for Dirichlet correlation
#' 
#' Mostly used as a helper function 
#'
#' @param alpha 
#' @param n 
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
get_dirichlet_cor <- function(alpha,n,p){
  id_c <- rep(1:n, each = p)
  #v <- var_dirichlet(alpha,n,p)
  
  # Make correlation block matrices for each subject
  cor_list <- list()
  for(i in 1:n){
    # or change this to be a matrix. 
    # TODO make more efficient
    alpha_i <- data.frame(alpha, id_c) %>% 
      filter(id_c == i) %>% 
      pull(alpha)
    # Dirichlet correlation is: 
    # Cor(i,j) = - 
    alpha0 <- sum(alpha_i)
    num <- - sqrt(alpha_i %*% t(alpha_i)) 
    denom <- sqrt((alpha0-alpha_i) %*% t((alpha0 -alpha_i)))
    cor_i <- num/denom
    
    # dirichlet covariance defined differently for diagonal, so 
    # to simplify coding, set correlation diagonal equal to 1. 
    diag(cor_i) <- 1
    
    cor_list[[i]] <- cor_i
  }
  
  return(cor_list)
}


#' Non-linear least squares Optimization 
#' 
#' Solve for omega and rho using non-linear least squares optimization 
#' Use the equation omega R_dir + (1 - omega) e^-2 rho D 
#'
#' @param resid_vec 
#' @param cor_vec 
#' @param D_vec 
#'
#' @return A list of parameters omega and rho 
#' @export
#'
#' @examples
nls_optim <- function(resid_vec, cor_vec, D_vec){
  par <- list(omega = .5, 
              rho = 1)
  fun <- function(par){
    sum((resid_vec - (par[1]*cor_vec + (1 - par[1])*exp(-2*par[2]*D_vec)))^2)
  }  
  
  estimates <- optim(par, fun,  lower=c(0,0), upper= c(1, Inf),method="L-BFGS-B")
  return(estimates$par)
}

