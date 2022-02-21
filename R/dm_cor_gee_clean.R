#' Phylogeny-aware Dirichlet GEE 
#' 
#' 
#' @param Y should be of length n*p
#' @param X should be dimension n*q
#' @param id should be of length n*p
#' @param distance_matrix 
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
  # Add intercept
  X_intercept <- data.frame(int = rep(1, n), 
                            X = X_vecs)
  X_mat <- as.matrix(X_intercept, nrow = n)
  
  X_big_i <- list()
  for( i in 1:n){
    X_big_i[[i]] <- matrix(
      bdiag(
        rep(list(matrix(X_mat[i,], nrow =1)),p)), nrow =p)
  }
  X <- reduce(X_big_i, rbind)
  
  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 
  beta_matrix <- matrix(0, nrow = q, ncol = p)
  # The Dirichlet link function is the log
  y_means <- matrix(Y, ncol = p)
  # is this link or inv_link? 
  beta_matrix[1,] <- log(colMeans(y_means))
  beta <- as.vector(beta_matrix)
  # Initialize the A matrix which is a diagonal matrix of the inverse of the square root of the variance
  # Here it is (nxp)x(nxp) 
  # This is using the Matrix package diagonal function...
  A <- Diagonal(n*p)
  
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
  
  
  # Main loop
  count <- 0
  while(count < 12){
    count <- count + 1
    print(paste0("Iteration: ", count))
    # eta is g(mu) = g(alpha) the link between mean response and covariates
    # eta = log(alpha)
    eta <- as.vector(X %*% beta)
    # Use the dirichlet link function that links mu to eta. (links alpha to beta)
    # First iteration will have all eta 0, alpha 1. 
    alpha <- exp(eta)
    alpha0 <- colSums(matrix(alpha, nrow = p))
    mu <- alpha / rep(alpha0, each = p)
    
    # Calculate variance from values of alpha 
    # since variance is function of alpha 
    diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
    resid <- diag(A %*% Diagonal(x = Y - mu))
    # do we need this?
    phi <- as.numeric(crossprod(resid,resid)*(1/(n-(q-1))))

    # Get dirichlet correlation part
    cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
    cor_dirichlet <- as.matrix(bdiag(cor_dirichlet_list))
    diag(cor_dirichlet) <- 1
    flat_cor_dir <- cor_dirichlet[upper.tri(cor_dirichlet)]
    flat_cor_dir <- flat_cor_dir[flat_cor_dir != 0]
    
    # Setup calculated residuals for nls fxn
    resids <- data.frame(id,resid) %>% 
      group_by(id) %>% 
      group_map(~tcrossprod(.x$resid))
    resids_sparse <- as.matrix(bdiag(resids))
    flat_resids <- resids_sparse[upper.tri(resids_sparse)]
    flat_resids <- flat_resids[flat_resids != 0]
    browser()
    # try to standardize with phi
    nls_data_frame <- data.frame(y = flat_resids/sqrt(phi),
                                 r = flat_cor_dir,
                                 d = d_ijk)  
    nls_mod <- nls(y ~ w*r + (1 - w)*exp(-2*rho*d), 
                   data = nls_data_frame, 
                   start = list(w = .5, 
                                rho = 1))
    res_nls <- summary(nls_mod)
    
    # save estimated omega and rho
    omega <- res_nls$coefficients[1,1]
    rho <- res_nls$coefficients[2,1]
    rho <- ifelse(rho < 0 , 0, rho)
    
    # use current values of omega and rho to create
    # present iteration of combined working correlation matrix. 
    
    Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))

    # invert
    # Use Moore-Penrose generalized inverse
    R_invs <- map(Rs, ginv)
    #R_invs <- map(Rs, solve)
    R_inv <- bdiag(R_invs)
    
    
    # Step for updating beta 
    beta.new <- beta
    diffs <- numeric(1)
    for(s in 1:1){
      print(paste0("Beta iteration ", s))
      beta.old <- beta.new
      eta <- as.vector(X%*%beta.new) 
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
      
      # Save V inv so doesnt need to be 
      # reminder A is A^{-1/2}, and A^t = A (its diagonal)
      V_inv <- A %*% R_inv %*% A
  
      hess <- partials %*% V_inv %*% t(partials)
      # GEE estimating equations/ gradient 
      esteq <- partials %*% V_inv %*% as.matrix(Y - mu)
      
      update <- solve(hess, esteq)
      # Should there be a minus here? 

      beta.new <- beta.new - as.vector(update)
      # Save to see speed of convergence
      diffs[s] <- sum((beta.new - beta.old)^2)
      print(paste0("beta",beta.new))
      
    }
    beta <- beta.new
    
    eta_list[[count]] <- eta
    alpha_list[[count]] <- alpha
    omega_list[[count]] <- omega
    rho_list[[count]] <- rho
    beta_list[[count]] <- beta
    beta_diffs[[count]] <- diffs
  }
  return(list(etas = list(eta_list), 
              alphas = list(alpha_list), 
              omegas = omega_list, 
              rhos = rho_list, 
              betas = list(beta_list), 
              beta_convergences = beta_diffs,
              num_iter = count))
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
  
  # Make covariance block matrices for each subject
  cor_list <- list()
  for(i in 1:n){
    # or change this to be a matrix. 
    # TODO make more efficient
    alpha_i <- data.frame(alpha, id_c) %>% 
      filter(id_c == i) %>% 
      pull(alpha)
    # Dirichlet covariance is: 
    # Cor(i,j) = - 
    alpha0 <- sum(alpha_i)
    num <- - sqrt(alpha_i %*% t(alpha_i)) 
    denom <- sqrt((alpha0-alpha_i) %*% t((alpha0 -alpha_i)))
    cor_i <- num/denom
    cor_list[[i]] <- cor_i
  }
  
  return(cor_list)
}


