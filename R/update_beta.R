#' Update beta
#'
#' Inner function for updating beta each iteration using current values of rho and omega. Calculates phi and alpha.
#'
#' @param Y
#' @param beta
#' @param ASV_id
#' @param phi
#' @param n_iter
#' @param X
#' @param n
#' @param p
#' @param q
#' @param rho
#' @param omega
#' @param D
#' @param gamma
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
update_beta <- function(Y, X, beta, ASV_id, phi, n_iter = 1,
                        n, p, q, rho, omega, D, gamma, lambda, fixed, A, R_inv) {
  # beta.new <- beta
  # G_init <-  calculate_equations(beta, n=n, p=p,q, Y=Y, X=X, hess = F, omega, rho, D,phi)$G

  # Save dummy G_new to start loop (LS)
  # G_new <- 0
  # G_prev <- 0
  count <- 0
  # Commented out condition for start of line search
  # while((count < 1 | sum(abs(G_new)) > sum(abs(G_prev))) & count < 5){
  # This is a loop but it is always just 1 iteration
  for (s in 1:n_iter) {
    count <- count + 1
    # beta.old <- beta.new
    # Calculate Hessian and GEE values based on prev beta value

    if(!fixed){
      eqns <- get_gee_equations(beta, n, p, q, Y, X, hess = T, omega, rho, D, lambda, fixed, ...)
    } else{
      eqns <- get_gee_equations(beta, n, p, q, Y, X, hess = T, omega, rho,
                                D, lambda, fixed = fixed, A=A, R_inv = R_inv)
    }
    
    
    hess <- eqns$H
    esteq <- eqns$G
    G_prev <- esteq

    update <- Matrix::solve(hess, esteq)
    # Calculate new beta value based on
    # beta+ = beta + gamma H^-1 G
    beta_new <- beta - gamma * as.vector(update)

    # Checks that GEE eqn is reduced
    G_new <- get_gee_equations(beta_new, n, p, q, Y, X, hess = F, omega, rho, D)$G
    print(paste0("Reduction:", sum(abs(G_new)) < sum(abs(G_prev))))
    print(paste0("G_new: ", sum(abs(G_new)), ", G_old: ", sum(abs(G_prev))))
  }

  beta <- beta_new
  print("Summary of beta")
  print(summary(beta))
  return(list(
    beta = beta,
    G = sum(abs(G_prev)),
    G_new = sum(abs(G_new))
  ))
}
