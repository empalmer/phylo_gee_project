
# Helper fxn
get_eta <- function(X, beta, n, p){
  X_list <- as.list(data.frame(t(X)))
  beta_mat <- matrix(beta, nrow = p, byrow = F)
  unname(unlist(map( X_list, ~ as.numeric(tcrossprod(.x, beta_mat)))))
}