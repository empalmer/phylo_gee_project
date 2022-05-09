#' Calculate eta 
#' 
#' Get eta based on values of X and beta (and n and p). eta = X^t beta. Helper function to calculate eta based on the same X value for each beta for different samples. 
#'
#' @param X 
#' @param beta 
#' @param n 
#' @param p 
#'
#' @return vector of eta values
#' @export
#'
#' @examples
get_eta <- function(X, beta, n, p){
  X_list <- as.list(data.frame(t(X)))
  beta_mat <- matrix(beta, nrow = p, byrow = F)
  unname(unlist(map( X_list, ~ as.numeric(tcrossprod(.x, beta_mat)))))
}