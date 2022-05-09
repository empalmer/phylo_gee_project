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
dm_cor_gee <- function(Y, X, id, distance_matrix, intercept = T){
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
  # minus one to accounting for intercept
  # TODO change when q > 1
  # q <- dim(X)[2] 
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
  # Start setting up the design matrix for GEE
  #Ip <- diag(p)
  # For use in the beta loop
  #Ip_n <- bdiag(rep(list(Ip),n))
  # kronecker product with the pxp identity matrix
  # Assume there is an intercept column.
  # X_no_intercept <- kronecker(X, Ip)
  # change when including intercept
  # X <- X_no_intercept
  #intercept <- rep(1,n*p)
  #X <- cbind(intercept, X_no_intercept)
  # Double check this
  #X <- kronecker(X_mat, Ip)
  # so big already??
  # Change to be x_i version? - not how matrix GEE setup is 
  # Big X is needed for eta <- X%>% beta in each loop
  # Later check if there is a speed/memory tradeoff here
  
  # Initialize beta column. Intercept beta0 is the mean of the Y, and the rest are 0. 
  # Just a helper index title to keep track 
  #beta_index <- c("intercept", paste0("otu",rep(1:p, each = q),
  #beta_index <- paste0("otu",rep(1:p, each = q),",cov",rep(1:q,p))
  #beta <- rep(0, p*q)
  beta_matrix <- matrix(0, nrow = q, ncol = p)
  #beta <- rep(0, p*q)
  # Fix this - 
  # The Dirichlet link function is the log
  # Set the intercept term
  # Is there just one? Or one for each OTU? 
  # 1/19 initialize the intercept term for each j
  # which is mean of the ys? 

  y_means <- matrix(Y, ncol = p)
  # is this link or inv_link? 
  beta_matrix[1,] <- link_fun(colMeans(y_means))
  beta <- as.vector(beta_matrix)
  # beta[1] <- link_fun(mean(Y))
  # Initialize the A matrix which is a diagonal matrix of the inverse of the square root of the variance
  # Here it is (nxp)x(nxp) 
  # This is using the Matrix package diagonal function...
  A <- Diagonal(n*p)
  
  # Is this the dimension when there are more covariates? 
  #partials <- Matrix(nrow = n*p, ncol = n*p)
  # Initialize the partials matrix 
  #partials <- Matrix(nrow = n*p*q, ncol = n*p)
  # dg_inv_d_eta <- Diagonal(n*p)
  
  # Is there an overdispersion parameter phi for dirichlet distribution? - zhao paper seems like it does? 
  # 1/19 Check histogram of residuals to see if overdispersion is necessary - 
  # maybe not since it is proportions not counts 
  
  # get 
  d_jk <- distance_matrix[upper.tri(distance_matrix)]
  # Since d_jk doesnt have an i index we need to repeat it for each sample
  d_ijk <- rep(d_jk,n)
  #d2_ijk <- 2*d_ijk
  
  # Set up storage to keep track of parameter values in each iteration 
  eta_list <- list()
  alpha_list <- list()
  omega_list <- list()
  rho_list <- list()
  beta_list <- list()
  beta_diffs <- list()
  
  
  # Main fisher scoring loop
  count <- 0
  while(count < 12){
    count <- count + 1
    print(paste0("Iteration: ", count))
    # eta is g(mu) = g(alpha) the link between mean response and covariates
    # eta = log(alpha)
    # do we make this a data frame too? 
    #eta <- map_dbl(1:p, ~ X_mat[.x,] %*% beta)
    eta <- as.vector(X %*% beta)
    # Use the dirichlet link function that links mu to eta. (links alpha to beta)
    #alpha is equivalent to mu in prev notation
    # (Include code of what the function link is)
    # First iteration will have all eta 0, alpha 1. 
    alpha <- inv_link_fun(eta)
    alpha0 <- colSums(matrix(alpha, nrow = p))
    mu <- alpha / rep(alpha0, each = p)

    # Calculate variance from values of alpha 
    # since variance is function of alpha 
    diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
    
    # is this right? 
    # TODO update if using weights. 

    # CHANGED from Y - alpha to Y - mu 
    resid <- diag(A %*% Diagonal(x = Y - mu))
    # q or q+1
    # do we need this?
    phi <- crossprod(resid,resid)*(1/(n-(q-1)))
    
    # Get dirichlet correlation part
    cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
    cor_dirichlet <- bdiag(cor_dirichlet_list)
    diag(cor_dirichlet) <- 1
    # somehow this line doesnt work??? 
    flat_cor_dir <- cor_dirichlet[upper.tri(cor_dirichlet)]
    flat_cor_dir <- flat_cor_dir[flat_cor_dir != 0]
    
    #initialize mult_e column
    #resids <- list()
  
    # new cross prod version 
    resids <- data.frame(id,resid) %>% 
      group_by(id) %>% 
      group_map(~tcrossprod(.x$resid))
    
    # for(i in 1:n){
    #   # get matrix of residuals for just one sample
    #   resid_i <- data.frame(resid, id) %>%
    #     filter(id == i) %>%
    #     pull(resid)
    #   # matrix of e_ij*e_ik
    #   resids[[i]] <- resid_i %*% t(resid_i)
    # }
    resids_sparse <- bdiag(resids)
    # find a more efficient way to do this. 
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
    rho <- ifelse(rho < 0 , 0, rho)
    
    # now unnecessary 
    # make block matrix of 
    # diag_D <- bdiag(rep(list(distance_matrix),n))
    # R_etm <- exp(-2*rho*diag_D)
    # 
    # use current values of omega and rho to create
    # present iteration of combined working correlation matrix. 
    
    Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*distance_matrix))
    for(i in 1:n){
      print(i)
      solve(Rs[[i]])
    }
    
    R_invs <- map(Rs, solve)
    
    
    R_inv <- bdiag(R_invs)
    
    
    # Step for updating beta 
    # TODO limit convergence steps for speedier algorithm? 
    
    beta.new <- beta
    diffs <- numeric(1)
    for(s in 1:1){
      print(paste0("Beta iteration ", s))
      beta.old <- beta.new
      eta <- as.vector(X%*%beta.new) 
      alpha <- inv_link_fun(eta) 
      alpha0 <- colSums(matrix(alpha, nrow = p))
      mu <-  alpha / rep(alpha0, each = p)
      
  
      diag(A) <- sqrt(1/var_dirichlet(alpha,n,p))
      #X <- X[,-1]
      # Make the partials for each block.
      partiali <- list()
      for(i in 1:n){
        alphai0 <- alpha0[i]
        #Xi <- X[((k-1)*p + 1):((k*p)),]
        # included intercept here FIXME
        
        # FIX to change X to a df for generality. 
        # breaks dimensions
        #xi <- c(1,X_vec_df[i])
        #xi <- X_vec_df[i]
        xi <- X_mat[i,]
        alphai <- alpha[((i-1)*p + 1):(i*p)]
        
        # this xi is the vector xi. 
        partiali[[i]] <- (1/alphai0)^2*kronecker((alphai0*diag(alphai) - 
                                                  alphai %*% t(alphai)),xi)
        #partiali[[i]] <- (1/alphai0)*(kronecker(diag(p), xi))%*%diag(alphai) - (1/alphai0)^2*kronecker(alphai %*% t(alphai), xi)

      }
      # NOT bdiag. Instead a long cbinded matrix of pqxnp
      partials <- purrr::reduce(partiali, cbind)
      #partials <- bdiag(partiali)
      
      # now no longer a diagonal matrix. 
      # derivative inverse link
      # for dirichlet distribution the derivative of exp(x) is exp(x)
      # diag(dg_inv_d_eta) <- exp(eta)
      
      # Save V inv so doesnt need to be 
      # reminder A is A^{-1/2}, and A^t = A (its diagonal)
      V_inv <- A %*% R_inv %*% A
      
      # not sure why crossprod is used here - hides the transposes 
      # maybe faster - yes? 
      #browser()
      # A = A^{-1/2}, simplified for notation here 
      #hess <- t(X) %*% t(partials) %*% A %*% R_inv %*% A %*% partials %*% X
      #hess <- crossprod(A %*% partials %*% X,
      #                   R_inv %*% A %*% partials %*% X)
      # Change to be crossprod() eventually so it runs faster?
      hess <- partials %*% V_inv %*% t(partials)
      #esteq <- crossprod(A %*% partials %*% X,
      #                    R_inv %*% A %*% as.matrix(Y - mu))
      # GEE estimating equations/ gradient 
      esteq <- partials %*% V_inv %*% as.matrix(Y - mu)
      
      update <- solve(hess, esteq)
      # Should there be a minus here? 
      # Try it ? 
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
    # or use crossprod function? might speed up 
    # cross <- alpha_i %*% t(alpha_i)
    # alpha0 <- sum(alpha_i)
    # denom <- -alpha0^2*(alpha0 + 1)
    # 
    # v_i <- data.frame(v,id_c) %>% 
    #   filter(id_c == i) %>% 
    #   pull(v)
    # V_prod <- v_i %*% t(v_i)
    # V_denom <- sqrt(V_prod)
    # 
    # cor_i <- cross/(denom *V_denom)
    cor_list[[i]] <- cor_i
  }

  return(cor_list)
}


