#' Phylogeny-aware Dirichlet GEE 
#' 
#' 
#' @param Y should be of length n*p
#' @param X should be dimension n*q
#' @param distance_matrix An object of type dist 
#' @param sample_id 
#' @param ASV_id 
#' @param intercept 
#' @param max_iter 
#' @param tol 
#' @param gamma 
#' @param lambda 
#' @param save_beta 
#'
#' @return
#' @export
#'
#' @examples
dm_cor_gee <- function(Y, X, sample_id, ASV_id,
                       distance_matrix, intercept = T, max_iter = 100, 
                       tol, gamma = 1, lambda = .01, save_beta = F){
  start.time <- Sys.time()
  require(tidyverse)
  require(Matrix)
  require(MASS)
  

# Source functions --------------------------------------------------------
  # Has get_eta function
  source(here::here("helpers.R"))
  # has get_dirichlet_var and get_dirichlet-cor
  source(here::here("dirichlet_functions.R"))
  # had update_beta
  source(here::here("update_beta.R"))
  # has calculate_equations, get_R_inv, calculate_partials
  source(here::here("gee_functions.R"))
  # has update_phi_rho_omega, nls_optim
  source(here::here("update_phi_rho_omega.R"))
  
# Initialize ----------------------------------------------------------

  function_call <- match.call()
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
  #beta_matrix[1,] <- log(y_means/sum(y_means))
  # convert matrix to vector
  beta <- as.vector(beta_matrix)
  
  # Set up distance matrix
  d_jk <- distance_matrix[upper.tri(distance_matrix)]
  # Since d_jk doesnt have an i index we need to repeat it for each sample
  d_ijk <- rep(d_jk,n)
  
  # Set up storage to keep track of parameter values in each iteration 
  res <- list()

# Main loop ---------------------------------------------------------------
  count <- 0
  diff <- 100
  while( diff > tol & count < max_iter){
    count <- count + 1
    print(paste0("Iteration: ", count))
    
# Step 1: R, w, rho ---------------------------------------------------------
    # wRrho_res <- update_phi_rho_omega(Y = Y, X = X, id = sample_id,
    #                              distance_matrix = distance_matrix,
    #                              d_ijk = d_ijk,
    #                              beta = beta, n = n, p = p, q = q)
    # phi <- wRrho_res$phi
    # rho <- wRrho_res$rho
    # omega <- wRrho_res$omega
    omega <- 1
    rho <- 5
    phi <- 1

# Step 2: Beta  -----------------------------------------------------
    # Depends on "fixed" values of rho, omega and phi, 
    # Which are used to make R_inv 
    beta.old <- beta
    vals <- update_beta(Y = Y, X = X, beta = beta, R_inv = R_inv,
                        phi = phi, n_iter = 1, n=n, p=p, q=q, ASV_id,
                        rho, omega, D = distance_matrix, gamma = gamma, 
                        lambda = lambda)
    beta <- vals$beta
    
# Convergence criteria and save results----------------------------------------
    diff <- sum((beta.old - beta)^2)
    print(paste0("Difference = ", diff))
    
    # Save results that are updated each loop
    res$omega[count] <- omega
    res$rho[count] <- rho
    res$diff[count] <- diff
    res$phi[count] <- phi
    res$G[count] <- vals$G
    res$G_new[count] <- vals$G_new
    if(save_beta){
      res$betas[count] <- list(beta)
    }
  }
  
  # Save results that are updated after last iteration.
  res$num_iter <- count
  res$time <- Sys.time() - start.time
  res$function_call <- function_call
  #res$resids <- wRrho_res$st_resid
  if(!save_beta){
    res$beta <- beta
  }
  
  return(res)
}


