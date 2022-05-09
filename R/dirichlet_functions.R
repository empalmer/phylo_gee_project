#' Function for calculating Dirichlet variance 
#' 
#' Helper function
#'
#' @param alpha a vector of n*p alpha values 
#' @param n a n*p vector of what sample each alpha is
#' @param p 
#'
#' @return Matrix of variances for each sample+ASV
#' @export
#' 
#' @examples
get_dirichlet_var <- function(alpha,n,p){
  id <- rep(1:n, each = p)
  data.frame(alpha, id) %>% 
    group_by(id) %>% 
    mutate(alpha0 = sum(alpha)) %>%
    ungroup() %>% 
    mutate(v = (alpha*(alpha0 - alpha))/(alpha0^2*(alpha0 + 1))) %>% 
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
    num <- -sqrt(alpha_i %*% t(alpha_i)) 
    denom <- sqrt((alpha0-alpha_i) %*% t((alpha0 -alpha_i)))
    cor_i <- num/denom
    
    # dirichlet covariance defined differently for diagonal, so 
    # to simplify coding, set correlation diagonal equal to 1. 
    diag(cor_i) <- 1
    
    cor_list[[i]] <- cor_i
  }
  
  return(cor_list)
}
