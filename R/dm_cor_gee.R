
# Can load practice data from https://github.com/lichen-lab/glmmTree in file glmmTree_data.R
source('~/Desktop/ResearchFall2021/R/glmmTree_data.R')
# Filtered to have 50 samples and 48 OTUs so computation goes quickly 

# switch and pass in Y and X separetly. 
# Y should be dimension n*p
# X should be dimension n*q
# id should be dimension n*p

#dm_cor_gee <- function(formula, data, id){
dm_cor_gee <- function(Y, X, id, distance_matrix){
  require(tidyverse)
  require(Matrix)
  require(brms) # For Dirichlet helper functions
  browser()
  # Save formula and data in correct format
  # old version
  # model_data <- model.frame(formula, data, na.action=na.pass)
  # id <- dplyr::pull(data, !!enquo(id))
  
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
  # X_no_intercept <- kronecker(X_compact[,-1], Ip)
  X_no_intercept <- kronecker(X, Ip)
  intercept <- rep(1,n*p)
  X <- cbind(intercept, X_no_intercept)
  # soo big already??
  
  # Y <- model.response(model_data)
  # Set up functions for dirichlet distribution
  dirichlet <- brms::dirichlet()
  
  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 
  beta_index <- c("intercept", paste0("otu",rep(1:p, each = q),",cov",rep(1:q,p)))
  beta <- rep(0, p*q + 1)
  # The dirichlet mean function is the logit. 
  beta[1] <- dirichlet$linkfun(mean(Y))
  
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
    # eta = log(alpha/(1-alpha))
    eta <- as.vector(X %*% beta)
    # Use the dirichlet link function that links mu to eta. (links alpha to beta)
    #alpha is equivalent to mu in prev notation
    # (Include code of what the function link is)
    alpha <- dirichlet$linkinv(eta)

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

    # test <- data.frame(resid,id) %>% 
    #   group_split() %>% 
    #   map(~.x$resid%*% t(.x$resid)) %>% 
    #   map(~.x[upper.tri(.x)])
    # View(test[[1]])
    #initialize mult_e column
    resids <- list()
    for(i in 1:n){
      # get matrix of residuals for just one sample
      resid_i <- data.frame(resid, id) %>%
        filter(id == i) %>%
        pull(resid)
      # matrix of e_ij*e_ik
      resids[[i]] <- resid_i %*% t(resid_i)
      # mult_e <- c(mult_e,resid_mat[upper.tri(resid_mat)])
    }
    resids_sparse <- bdiag(resids)
    flat_resids <- resids_sparse[upper.tri(resids_sparse)]
    flat_resids <- flat_resids[flat_resids != 0]
    
    # combos <- expand_grid(j = 1:p, k = 1:p)
    # indeces <- map_dfr(1:n, ~cbind(.x*combos,combos))
    # colnames(indeces) <- c("jn","kn","j","k")
    # 
    # big_resid_mat <- resid %*% t(resid)
    # 
    # mult_e <- big_resid_mat[indeces$jn, indeces$kn]
    # cor_dir <- cor_dirichlet[indeces$jn, indeces$kn]
     

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
      alpha <- dirichlet$linkinv(eta) #mu <- InvLink(eta)
      diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
      # derivative inverse link
      diag(dg_inv_d_eta) <- exp(eta)/(1 + exp(eta))^2
      
      
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
  # combine into block matrix
  # cor_mat <- bdiag(cor_list)
  # diag(cor_mat) <- 1
  return(cor_list)
}


