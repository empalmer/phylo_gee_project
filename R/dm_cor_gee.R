
# switch and pass in Y and X separetly. 
# Y should be dimension n*p
# X should be dimension n*q
# id should be dimension n*p


dm_cor_gee <- function(Y, X, id, distance_matrix){
  require(tidyverse)
  require(Matrix)
  
  # For dirichlet distribution use log link (g(mu))
  # Save inverse link function 
  link_fun <- function(x){log(x)}
  inv_link_fun <- function(x){exp(x)}
  
  # Set up indeces 
  # Number of subjects i = 1,...,n
  n <- length(unique(id))
  # Number of ASVs j = 1,...,p
  # p <- n
  # TODO change when we allow unequal cluster size 
  p <- length(Y)/n
  # Convert into correct model matrix format 
  # X compact is the model matrix of just the given data
  # Not yet in the correct format for GEE 
  # X_compact <- model.matrix(formula,model_data)
  # Number of covariates: k = 1,...,q
  # minus one to accouting for intercept
  # TODO change when q > 1
  # q <- dim(X)[2] 
  q <- 1
  
  # Start setting up the design matrix for GEE
  Ip <- diag(p)
  # kronecker product with the pxp identity matrix
  # Assume there is an intercept column.
  X_no_intercept <- kronecker(X, Ip)
  intercept <- rep(1,n*p)
  X <- cbind(intercept, X_no_intercept)
  # soo big already??
  
  
  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 
  beta_index <- c("intercept", paste0("otu",rep(1:p, each = q),",cov",rep(1:q,p)))
  beta <- rep(0, p*q + 1)
  # The Dirichlet link function is the log
   beta[1] <- link_fun(mean(Y))
  
  # Initialize the A matrix which is a diagonal matrix of the inverse of the square root of the variance
  # Here it is (nxp)x(nxp) 
  # This is using the Matrix package diagonal function...
  A <- Diagonal(n*p)
  
  dg_inv_d_eta <- Diagonal(n*p)
  
  # Is there an overdispersion parameter phi for dirichlet distribution? - zhao paper seems like it does? 
  
  # Main fisher scoring loop
  count <- 0
  while(count < 30){
    count <- count + 1
    
    # eta is g(mu) = g(alpha) the link between mean response and covariates
    # eta = log(alpha)
    eta <- as.vector(X %*% beta)
    # Use the dirichlet link function that links mu to eta. (links alpha to beta)
    #alpha is equivalent to mu in prev notation
    # (Include code of what the function link is)
    alpha <- inv_link_fun(eta)

    # Calculate variance from values of alpha 
    # since variance is function of alpha 
    diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
    
    # is this right? 
    # TODO update if using weights. 

    resid <- diag(A %*% Diagonal(x = Y - alpha))
    # q or q+1
    # do we need this?
    phi <- crossprod(resid,resid)*(1/(n-q))
    
    # Get dirichlet correlation part
    cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
    cor_dirichlet <- bdiag(cor_dirichlet_list)
    diag(cor_dirichlet) <- 1
    flat_cor_dir <- cor_dirichlet[upper.tri(cor_dirichlet)]
    flat_cor_dir <- flat_cor_dir[flat_cor_dir != 0]
    
    # get 
    d_jk <- distance_matrix[upper.tri(distance_matrix)]
    # Since d_jk doesnt have an i index we need to repeat it for each sample
    d_ijk <- rep(d_jk,n)
    d2_ijk <- 2*d_ijk

    #initialize mult_e column
    resids <- list()
    for(i in 1:n){
      # get matrix of residuals for just one sample
      resid_i <- data.frame(resid, id) %>%
        filter(id == i) %>%
        pull(resid)
      # matrix of e_ij*e_ik
      resids[[i]] <- resid_i %*% t(resid_i)
    }
    resids_sparse <- bdiag(resids)
    flat_resids <- resids_sparse[upper.tri(resids_sparse)]
    flat_resids <- flat_resids[flat_resids != 0]
    
    nls_data_frame <- data.frame(y = flat_resids,
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
    
    # now unnecessary 
    # make block matrix of 
    # diag_D <- bdiag(rep(list(distance_matrix),n))
    # R_etm <- exp(-2*rho*diag_D)
    # 
    # use current values of omega and rho to create
    # present iteration of combined working correlation matrix. 
    
    Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))
    
    R_invs <- map(Rs, solve)
    
    
    R_inv_block <- bdiag(R_invs)
    
    
    # Step for updating beta 
    # limit convergence steps for speedier algorithm? 
    
    beta.new <- beta
    for(i in 1:10){
      eta <- as.vector(X%*%beta.new) 
      alpha <- inv_link_fun(eta) #mu <- InvLink(eta)
      diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
      # derivative inverse link
      # for dirichlet distribution the derivative of exp(x) is exp(x)
      diag(dg_inv_d_eta) <- exp(eta)
      
      
      hess <- crossprod(A %*% dg_inv_d_eta %*% X,
                        R_inv_block %*% A %*% dg_inv_d_eta %*% X)
      esteq <- crossprod(A %*% dg_inv_d_eta %*% X,
                         R_inv_block %*% A %*% as.matrix(Y - alpha))
      update <- solve(hess, esteq)
      beta.new <- beta.new + as.vector(update)

    }
  }
}


# alphas is a vector of n*p alpha values 
# id is a n*p vector of what sample each alpha is
# returns the matrix of variances for all 
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


get_dirichlet_cor <- function(alpha,n,p){
  id_c <- rep(1:n, each = p)
  v <- var_dirichlet(alpha,n,p)
  
  # Make covariance block matrices for each subject
  cor_list <- list()
  for(i in 1:n){
    alpha_i <- data.frame(alpha, id_c) %>% 
      filter(id_c == i) %>% 
      pull(alpha)
    
    cross <- alpha_i %*% t(alpha_i)
    alpha0 <- sum(alpha_i)
    denom <- -alpha0^2*(alpha0 + 1)
    
    v_i <- data.frame(v,id_c) %>% 
      filter(id == i) %>% 
      pull(v)
    V_prod <- v_i %*% t(v_i)
    V_denom <- sqrt(V_prod)
    
    cor_i <- cross/(denom *V_denom)
    cor_list[[i]] <- cor_i
  }

  return(cor_list)
}


