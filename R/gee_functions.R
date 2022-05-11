
#' Title
#'
#' @param beta
#' @param n
#' @param p
#' @param q
#' @param Y
#' @param X
#' @param hess
#' @param omega
#' @param rho
#' @param D
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
get_gee_equations <- function(beta, n, p, q, Y, X,
                              hess = T, omega, rho, D, lambda, fixed = F, A, R_inv) {
  H <- NULL

  eta <- get_eta(X, beta, n, p)
  alpha <- exp(eta)
  alpha0 <- rowSums(matrix(alpha, nrow = n, ncol = p, byrow = T))
  mu <- alpha / rep(alpha0, each = p)
  # # A is A^{-1/2} A is dirichlet variance
  # FIXME when not fixed
  if (!fixed) {
    A <- Diagonal(n * p)
    diag(A) <- sqrt(1 / get_dirichlet_var(alpha, n, p))
    R_inv <- get_R_inv(alpha, omega, rho, D, n, p)
  }


  # If independence  --------------------------------------------------------
  # R_inv <- diag(n*p)

  partials <- calculate_partials(alpha, alpha0, n, p, X)


  resid <- diag(A %*% Diagonal(x = Y - mu))
  # Overdispersion
  # phi = 1/(sum sum r^2)/(N-p)
  # sum of squared residuals
  # FIXME 
  phi <- 1
  #phi <- 1 / (as.numeric(sum(resid^2) * (1 / (n * p - (p * q - 1)))))

  # Save V inv
  # reminder A is A^{-1/2}, and A^t = A (A is diagonal)
  # V <- 1/phi * A %*% R %*% A
  V_inv <- phi * A %*% R_inv %*% A


  G <- partials %*% V_inv %*% as.matrix(Y - mu)
  if (hess) {
    # Since we have more parameters than samples
    # likely hessian matrix will be singular
    # add a small diagonal lambda to hessian matrix.
    # start line search
    # examine the hessian matrix.
    H <- -partials %*% V_inv %*% t(partials)
    # lambda based on eigenvalues?
    # res <- eigen(H)
    # lambda <- max(abs(res$values))/1000

    H <- H - diag(rep(lambda, q * p))
  }

  return(list(
    G = G,
    H = H
  ))
}






#' get_R_inv
#'
#' @param alpha
#' @param omega
#' @param rho
#' @param D
#' @param n
#' @param p
#'
#' @return
#' @export
#'
#' @examples
get_R_inv <- function(alpha, omega, rho, D, n, p) {
  cor_dirichlet_list <- get_dirichlet_cor(alpha, n, p)
  Rs <- cor_dirichlet_list
  # FIXME when both kinds of correlation
  # Rs <- purrr::map(cor_dirichlet_list, ~ .x*omega + (1-omega)*exp(-2*rho*D))
  # Use Moore-Penrose generalized inverse
  R_invs <- map(Rs, ginv)
  R_inv <- bdiag(R_invs)
}


#' calculate_partials
#'
#' @param alpha
#' @param alpha0
#' @param n
#' @param p
#' @param X
#'
#' @return
#' @export
#'
#' @examples
calculate_partials <- function(alpha, alpha0, n, p, X) {
  # Make the partials for each block.
  partiali <- list()
  for (i in 1:n) {
    alphai0 <- alpha0[i]
    xi <- X[i, ]
    alphai <- alpha[((i - 1) * p + 1):(i * p)]
    # this xi is the vector xi.
    partiali[[i]] <- (1 / alphai0)^2 * kronecker((alphai0 * diag(alphai) -
      alphai %*% t(alphai)), xi)
  }
  # dimension of partials is pq * np
  partials <- purrr::reduce(partiali, cbind)
  return(partials)
}
