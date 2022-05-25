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
                       tol, gamma = 1, lambda = .01, save_beta = F, 
                       only_dir_cor = F, fixed = F, A, R_inv, 
                       init_beta,...) {
  start.time <- Sys.time()
  require(tidyverse)
  require(Matrix)
  require(MASS)

  
  # Source functions --------------------------------------------------------
  # Has get_eta function

  source(here::here("R","helpers.R"))
  # has get_dirichlet_var and get_dirichlet-cor
  source(here::here("R","dirichlet_functions.R"))
  # had update_beta
  source(here::here("R","update_beta.R"))
  # has calculate_equations, get_R_inv, calculate_partials
  source(here::here("R","gee_functions.R"))
  # has update_phi_rho_omega, nls_optim
  source(here::here("R","update_phi_rho_omega.R"))

  # Initialize ----------------------------------------------------------

  function_call <- match.call()
  # Set up indeces
  # Number of subjects i = 1,...,n
  n <- length(unique(sample_id))
  # Number of ASVs j = 1,...,p
  p <- length(Y) / n
  # Number of covariates: k = 1,...,q
  q <- 1 + intercept

  # Add intercept, so a dataframe of the Xs plus a col of 1
  if (intercept) {
    X <- model.matrix(~x, data.frame(x = X))
  }

  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0.
  #beta_matrix <- matrix(0, nrow = q, ncol = p)
  # The Dirichlet link function is the log
  # y_means <- colMeans(matrix(Y, ncol = p))
  # Initialize the beta intercept term as ybar/ sum(ybar) = 1
  # beta_matrix[1,] <- log(y_means/sum(y_means))
  # convert matrix to vector
  #beta <- as.vector(beta_matrix)
  beta <- init_beta
  
  # Set up distance matrix
  d_jk <- distance_matrix[upper.tri(distance_matrix)]
  # Since d_jk doesnt have an i index we need to repeat it for each sample
  d_ijk <- rep(d_jk, n)

  # Set up storage to keep track of parameter values in each iteration
  results <- list()

  # Main loop ---------------------------------------------------------------
  count <- 0
  diff <- 100
  while (diff > tol & count < max_iter) {
    count <- count + 1
    print(paste0("Iteration: ", count))


    # Step 1: R, w, rho ---------------------------------------------------------
    if(!only_dir_cor){
      wRrho_res <- update_phi_rho_omega(Y = Y, X = X, id = sample_id,
                                        distance_matrix = distance_matrix,
                                        d_ijk = d_ijk,
                                        beta = beta, n = n, p = p, q = q)
      phi <- wRrho_res$phi
      rho <- wRrho_res$rho
      omega <- wRrho_res$omega
    } else{
      omega <- 1
      rho <- 5
      phi <- 1
    }


    # Step 2: Beta  -----------------------------------------------------
    # Depends on "fixed" values of rho, omega and phi,
    # Which are used to make R_inv
    beta.old <- beta
    
    if(fixed){
      beta_step <- update_beta(
        Y = Y, X = X, beta = beta,
        ASV_id = ASV_id, n_iter = 1,
        n = n, p = p, q = q, 
        rho, omega,
        D = distance_matrix, 
        gamma = gamma, lambda = lambda, A = A, R_inv = R_inv, fixed = T
      )
    }else{
      beta_step <- update_beta(
        Y = Y, X = X, beta = beta,
        ASV_id = ASV_id, n_iter = 1,
        n = n, p = p, q = q, 
        rho, omega,
        D = distance_matrix, 
        gamma = gamma, lambda = lambda, fixed = F
      )
    }
       
    
    beta <- beta_step$beta

    # Convergence criteria and save results----------------------------------------
    diff <- sum((beta.old - beta)^2)
    print(paste0("Difference = ", diff))

    # Save results that are updated each loop
    results$omega[count] <- omega
    results$rho[count] <- rho
    results$diff[count] <- diff
    results$phi[count] <- phi
    results$G[count] <- beta_step$G
    results$G_new[count] <- beta_step$G_new
    if (save_beta) {
      results$betas[count] <- list(beta)
    }
  }

  # Save results that are updated after last iteration.
  results$num_iter <- count
  results$time <- Sys.time() - start.time
  results$function_call <- function_call
  # res$resids <- wRrho_res$st_resid
  if (!save_beta) {
    results$beta <- beta
  }

  return(results)
}
