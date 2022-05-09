#' Update rho, omega, and R inv 
#'
#' @param X 
#' @param beta 
#' @param Y 
#' @param id 
#' @param distance_matrix 
#' @param d_ijk 
#' @param n 
#' @param p 
#' @param q 
#'
#' @return phi , omega , rho, standardized residuals
#' @export
#'
#' @examples
update_phi_rho_omega <- function(Y, X, id, distance_matrix, d_ijk,
                                 beta, n, p, q){
  # Setup calculations
  eta <- get_eta(X, beta, n, p)
  alpha <- exp(eta) 
  alpha0 <- colSums(matrix(alpha, nrow = p, byrow = T))
  mu <-  alpha / rep(alpha0, each = p)
  A <- Diagonal(n*p)
  diag(A) <- sqrt(1/get_dirichlet_var(alpha,n,p))
  
  if(any(is.infinite(diag(A))) | any(is.nan(diag(A)))) {
    stop("Stop condition met: Infinite values from alpha or A")
  }
  
  ### Now update omega, rho, 
  # We use the residuals as the responses for the NLS regression
  resid <- diag(A %*% Diagonal(x = Y - mu))
  
  # Overdispersion 
  # phi = 1/(sum sum r^2)/(N-p)
  # sum of squared residuals 
  # Note this is the inverse of some definitions of phy 
  # Using the original GEE paper definition
  phi <- 1/(as.numeric(sum(resid^2)*(1/(n*p-(p*q-1)))))
  
  # Setup calculated residuals for nls fxn
  # Need to "flatten" so each residual in the residual matrix 
  # These are squared residuals 
  cross_resids <- data.frame(id,resid) %>% 
    group_by(id) %>% 
    group_map(~tcrossprod(.x$resid))
  cross_resid_vec <- unlist(map(cross_resids, ~.x[upper.tri(.x)]))
  
  # Get dirichlet correlation part
  cor_dirichlet_list <- get_dirichlet_cor(alpha,n,p)
  cor_dirichlet_vec <- unlist(map(cor_dirichlet_list, ~.x[upper.tri(.x)]))
  
  # standardize cross prod residuals by multiplying by phi. 
  # since both are multiplied by sqrt phy. 
  st_resid <- sqrt(phi)*resid
  st_cross_resid <- phi*cross_resid_vec
  estimates <- nls_optim(st_cross_resid, cor_dirichlet_vec, d_ijk)
  
  omega <- estimates[1]
  rho <- estimates[2]
  print(paste0("phi = ", phi,", omega  = ", omega, " , rho = ", rho))
  
  return(list(phi = phi, 
              omega = omega, 
              rho = rho, 
              st_resid = st_resid))
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
  # This function is the sum of squares of 
  fun <- function(par){
    sum((resid_vec - (par[1]*cor_vec + (1 - par[1])*exp(-2*par[2]*D_vec)))^2)
  }  
  
  estimates <- optim(par, fun,  lower=c(0,0), upper= c(1, Inf),
                     method="L-BFGS-B")
  
  return(estimates$par)
}

