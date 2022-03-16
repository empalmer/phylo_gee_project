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
dm_cor_gee <- function(Y, X, sample_id, ASV_id,
                       distance_matrix, intercept = T){
  require(tidyverse)
  require(Matrix)
  require(MASS)
  
  # Set up indeces 
  # Number of subjects i = 1,...,n
  n <- length(unique(sample_id))
  # Number of ASVs j = 1,...,p
  p <- length(Y)/n
  # Number of covariates: k = 1,...,q
  q <- 1 + intercept
  
  
  # Add intercept, so a dataframe of the Xs plus a col of 1
  if(intercept){
    X <- model.matrix(~x, data.frame(x = X))
  }


  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 
  beta_matrix <- matrix(0, nrow = q, ncol = p)
  # The Dirichlet link function is the log
  y_means <- colMeans(matrix(Y, ncol = p))
  # Initialize the beta intercept term as ybar/ sum(ybar) = 1
  beta_matrix[1,] <- log(y_means/sum(y_means))
  # convert matrix to vector
  beta <- as.vector(beta_matrix)
  
  
  # Set up distance matrix
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
  phi_list <- list()
  
  # Initial values for rho, omega, phi 
  rho <- 1
  omega <- 0.5
  phi <- 1
  
  # Initialize R_inv 
  R_inv <- get_R_inv(X, beta, omega, rho, D, n, p)
  
  # Main loop
  count <- 0
  diff <- 100
  while( diff > .1 & count < 25){
    count <- count + 1
    print(paste0("Iteration: ", count))
    
    # Update beta loop 
    # Depends on "fixed" values of rho, omega and phi, 
    # Which are used to make R_inv 
    beta.old <- beta
    beta <- update_beta(Y = Y, X = X, beta = beta, R_inv = R_inv,
                        phi = phi, n_iter = 1, n=n, p=p, q=q, ASV_id)
  
    # Update R inverse by updating omega and rho 
    temp_res <- update_phi_R_inv(Y = Y, X = X, id = sample_id,
                                 distance_matrix = distance_matrix,
                                 d_ijk = d_ijk,
                                 beta = beta, n = n, p = p, q = q)
    R_inv <- temp_res$R_inv
    phi <- temp_res$phi
    
    
    diff <- sum(abs(beta.old - beta))
    print(paste0("Difference = ", diff))
    # Save estimates when things start working 
    # eta_list[[count]] <- eta
    # alpha_list[[count]] <- alpha
     omega_list[[count]] <- temp_res$omega
     rho_list[[count]] <- temp_res$rho
     beta_list[[count]] <- beta
     beta_diffs[[count]] <- diff
     phi_list[[count]] <- phi
    
  }
  return(list(betas = list(beta_list),
              omegas = unlist(omega_list), 
              rhos = unlist(rho_list), 
              phis = unlist(phi_list),
              beta_convergences = unlist(beta_diffs),
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
update_phi_R_inv <- function(Y, X, id, distance_matrix, d_ijk, beta, n, p, q){
  
  eta <- get_eta(X, beta, n, p)
  alpha <- exp(eta) 
  alpha0 <- colSums(matrix(alpha, nrow = p))
  mu <-  alpha / rep(alpha0, each = p)
  A <- Diagonal(n*p)
  diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
  
  if(any(is.infinite(diag(A))) | any(is.nan(diag(A)))) {
    stop("Infinite values due to infinite valued alpha or A")
  }
  
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
  
  # Overdispersion 
  # This is actually phi inverse. According to GEE paper
  # sum of squared residuals 
  phi <- as.numeric(sum(resid^2)*(1/(n-(q-1))))
  print(paste0("phi = ", phi))
  
  # Get dirichlet correlation part
  cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
  cor_dirichlet_vec <- unlist(map(cor_dirichlet_list, ~.x[upper.tri(.x)]))
  
  # standardize cross prod residuals by dividing by phi. 
  # CHECK THIS
  estimates <- nls_optim(cross_resid_vec*phi, cor_dirichlet_vec, d_ijk)
  
  omega <- estimates[1]
  rho <- estimates[2]
  print(paste0("omega  = ", omega, " , rho = ", rho))
  
  # use current values of omega and rho to create
  # present iteration of combined working correlation matrix. 
  Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))
  
  # invert using Moore-Penrose generalized inverse (solve will not work)
  R_invs <- map(Rs, ginv)
  R_inv <- bdiag(R_invs)
  
  return(list(phi = phi, 
              R_inv = R_inv, 
              omega = omega, 
              rho = rho))
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
update_beta <- function(Y, X, beta, R_inv, phi, n_iter = 1, n, p, q, ASV_id){
  beta.new <- beta
  diffs <- numeric(1)
  A <- Diagonal(n*p)
  for(s in 1:n_iter){
    print(paste0("Beta iteration ", s))
    # Used for convergence criteria 
    beta.old <- beta.new
    
    eta <- get_eta(X, beta, n, p)
    # eta <- as.vector(X %*% beta.new) 
    alpha <- exp(eta) 
    alpha0 <- colSums(matrix(alpha, nrow = p))
    mu <-  alpha / rep(alpha0, each = p)
    
    # A is A^{-1/2}
    diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
    
    # Make the partials for each block.
    partiali <- list()
    for(i in 1:n){
      alphai0 <- alpha0[i]
      xi <- X[i,]
      alphai <- alpha[((i-1)*p + 1):(i*p)]
      
      # this xi is the vector xi. 
      partiali[[i]] <- (1/alphai0)^2*kronecker((alphai0*diag(alphai) - 
                                                  alphai %*% t(alphai)),xi)
    }
    # dimension of partials is pq * np 
    partials <- purrr::reduce(partiali, cbind)
    
    # Save V inv 
    # reminder A is A^{-1/2}, and A^t = A (A is diagonal)
    # V_inv <- 1/phi * A %*% R_inv %*% A
    # Did i mess up the phis? 
    V_inv <- (1/phi) * A %*% R_inv %*% A
    

    # Since we have more parameters than samples 
    # likely hessian matrix will be singular 
    # add a small diagonal lambda to hessian matrix. 
    hess <- partials %*% V_inv %*% t(partials) + diag(rep(.001, q*p))
    # GEE estimating equations/ gradient 
    esteq <- partials %*% V_inv %*% as.matrix(Y - mu)
    
    update <- solve(hess, esteq)
    beta.new <- beta.new - as.vector(update)
    
    # Save to see speed of convergence
    diffs[s] <- sum((beta.new - beta.old)^2)

    print(paste0("beta ",  rep(ASV_id, each = q), " = ",beta.new))
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



get_eta <- function(X, beta, n, p){
  X_list <- as.list(data.frame(t(X)))
  beta_mat <- matrix(beta, nrow = p, byrow = F)
  unname(unlist(map( X_list, ~ as.numeric(tcrossprod(.x, beta_mat)))))
}

get_R_inv <- function(X, beta, omega, rho, D, n, p){
  # Initialize eta 
  eta <- get_eta(X, beta, n, p)
  alpha <- exp(eta) 
  alpha0 <- colSums(matrix(alpha, nrow = p))
  
  
  cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
  Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*D))
  # Use Moore-Penrose generalized inverse
  R_invs <- map(Rs, ginv)
  R_inv <- bdiag(R_invs)
}