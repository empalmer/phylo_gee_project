library(tidyverse)
library(Matrix)
# work out with mb example
formula <- Y ~ X 
id <- mb_data$sampleID
data <- mb_data
distance_matrix <- D_filtered

# for filtered data (50 samples)
test_results <- gee_mb(formula = formula, id = id, data = data,
                       family = gaussian(), distance_matrix = distance_matrix)
test_results

# for full data (100 samples)
test_results_full <- gee_mb(Y~X, id = mb_data_full$sampleID, data = mb_data_full,
                       family = binomial(), distance_matrix = D_filtered)
test_results_full



# GEE function
# 
# id identifies the clusters/subjects
gee_mb <- function(formula, id, data = parent.frame(), distance_matrix){
  # dont need family argument since we have decided to use dirichlet distribution 
  # although maybe change to allow dirichlet multinomial 
  data <- model.frame(formula, data, na.action=na.pass)
  data$id <- id
  
  dirichlet_family <- brms::dirichlet()
  # total number of observations (clusters x number in cluster)/ #subjects * #OTU
  N <- dim(data)[1]
  
  # n is the number of clusters - ie the number of subjects
  # 50 for mb example
  n <- length(unique(id))
  #p is the number of OTUs # 43 for mini example
  p <- nrow(distance_matrix)

  X <- model.matrix(formula,data)
  Y <- model.response(data)
  
  
  
  # The number of parameters in the model
  q <- dim(X)[2]
  q2 <- q*p
  # Change so it also includes dimension for beta for the p OTUs. 
  
  # Check if there is an intercept column. 
  interceptcol <- apply(X==1, 2, all)
  
  
  # Save link functions from specified family
  # eg. gaussian(), binomial() etc. 
  # link_mean_y <- family$linkfun(mean(Y))
  # InvLink <- function(x){family$linkinv(x)}
  # VarFun <- function(x){family$variance(x)}
  # InvLinkDeriv <- function(x){family$mu.eta(x)}
  
  # Save functions just for dirichlet distribution
  link_mean_y <- dirichlet_family$linkfun(mean(Y))
  
  # initialize alpha to just be the Ys
  alpha_init <- Y
  # get total sum of the alpha in each group ie read depth of the sample
  alpha_o <- data.frame(alpha_init, id) %>% 
    group_by(id) %>% 
    summarize(mean = mean(alpha_init)) %>% 
    pull(mean)
  # Initialize the beta column 
  # The intercept column will be initialized to the (link of) the mean of the Y
  # remaining betas will be 0 initially. 
  beta <- rep(0, dim(X)[2])
  beta_j <- rep(0, q2)
  beta[which(interceptcol)] <- link_mean_y
  
  # Number of observations in each cluster 
  # 43 for mb example
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  # Set initial rho value 
  # Is this reasonable??
  rho <- 1
  # Set overdispersion as 1 
  phi <- 1
  
  # Initialize the R working correlation matrix with the distance matrix and the initial value of rho 
  # R_{ij} = e^{rho d_{ij}}
  R_etm.init <- exp(-2*rho*distance_matrix)
  # How to create working dirichlet correlation matrix? 
  R_d <- - alpha_init*alpha_init
  
  # Set up matrix storage
  # A is the diagonal matrix with kth diagonal element equal to a(mu_ik) (stderr)
  # where V(Yik) = phi a(mu_{ik})
  A <- Diagonal(N)
  dInvLinkedEta <- Diagonal(N)
  derivmu <- dInvLinkedEta 
  Resid <- Diagonal(N)

  stop <- F
  converged <- F
  count <- 0
  beta.old <- beta
  unstable <- F
  phi.old <- phi
  # Main fisher scoring loop
  # dummy loop for now 
  #while(!stop){
  while(count < 30){
    count <- count + 1
    
    eta <- as.vector(X %*% beta)
    mu <- InvLink(eta)
    # A is actually A^{-1/2}. 
    diag(A) <- sqrt(1/VarFun(mu))
    
    
    # Step for updating phi   
    resid <- diag(A %*% sqrtW %*% Diagonal(x = Y - mu))
    phi <- (1/(sum(included)- q))*crossprod(resid, resid)
    phi.new <- phi
    
    
    resid_df <- data.frame(resid = resid, 
                           sample = rep(1:n,each = p))
    ggplot(resid_df, 
           aes(x = resid, group = sample)) + 
      geom_density()
    
    hist(resid)
    sd(resid)
    # Update alpha/rho step 
    # Set up regression of form eij*eik ~ e^-2rhodjk
    ut <- upper.tri(distance_matrix)
    d_jk <- distance_matrix[ut]
    
    # set up pairs to do the multiplication
    combo <- cbind(j = row(distance_matrix)[ut],
                   k = col(distance_matrix)[ut])
    add <- ((1:n)-1)*p
    index_combos <- map_dfr(add, ~ as.data.frame(combo + .x))
    
    # e_(ij)e_(ik):
    mult_e <- resid[index_combos[,1]] * resid[index_combos[,2]]
    
    
    
    hist(mult_e)
    
    alpha_df <- data.frame(y = mult_e, x = d_jk, # straight values
                           y2 = log(abs(mult_e)), x2 = -2*d_jk, # if we were doing linear
                           y3 = abs(mult_e), x3 = -2*d_jk) # take absolute value)
    #NaNs produced since we are taking the log... 
    # in this small case about half are negative. 
    #sum(is.nan(log(mult_e)))
    # model no intercept
    mod_lm <- lm(y2 ~ x2 - 1, data = alpha_df)
    rho <- summary(mod_lm)$coefficients[1]
    #mod_nls <-  nls(y3 ~ exp(x*a), data = alpha_df, start = c(a = 5))
    # update rho 
    #rho <- summary(mod)$coefficients[1]
    
    #updated correlation matrix based on new rho 
    R <- exp(-2*rho*distance_matrix)
    
    # matrix is not invertable? 
    # Error in solve.default(R) : 
    #system is computationally singular: reciprocal condition number = 5.97453e-27
    # Works in one case when rho is positive.
    # only for the case when cluster sizes are equal
    Rinv <- solve(R)
    R_inv_vec <- as.numeric(Rinv)
    
    # make into block matrix for everything else. 
    big_R_inv <- getBlockDiag(len, xvec = R_inv_vec)
    big_R_inv_bdiag <- big_R_inv$BDiag
    
    # Step for updating beta   
    beta.new <- beta
    # limit convergence steps for speedier algorithm? 
    for(i in 1:10){
      eta <- as.vector(X%*%beta.new) 
      diag(dInvLinkdEta) <- InvLinkDeriv(eta)
      mu <- InvLink(eta)
      diag(A) <- sqrt(1/VarFun(mu))
      
      # need to make sure that big_R_inv is a matrix? 
      # or sparse matrix. 
      # made it a bdiag format. seems to run. maybe not most efficient?
      
      # I dont think this would need to change?
      hess <- crossprod(A %*% dInvLinkdEta %*% X,
                        big_R_inv_bdiag %*% W %*% A %*% dInvLinkdEta %*% X)
      esteq <- crossprod(A %*% dInvLinkdEta %*% X,
                         big_R_inv_bdiag %*% W %*% A %*% as.matrix(Y - mu))
      update <- solve(hess, esteq)
      
      # Should this be minus??? 
      # no, this is when hessian optimization is used 
      beta.new <- beta.new + as.vector(update)
      
    }
    
    beta <- beta.new
    phi.old <- phi
    eta <- as.vector(X %*% beta)

  }
  return(list(beta = beta, rho = rho, phi = phi))
}



### Get a block diagonal matrix. Each block has dimension corresponding to
### each cluster size.  By default, each block is just a matrix filled with ones.
getBlockDiag <- function(len, xvec=NULL){
  K <- length(len)
  
  if(is.null(xvec)){
    xvec <- rep.int(1, sum(len^2))
  }
  
  row.vec <- col.vec <- vector("numeric", sum(len^2))
  add.vec <- cumsum(len) - len
  if(K == 1){
    index <- c(0, sum(len^2))
  }else{
    index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  }
  
  for(i in 1:K){
    row.vec[(index[i] + 1):(index[i+1])] <- rep.int( (1:len[i]) + add.vec[i], len[i])
    col.vec[(index[i] + 1):(index[i+1])] <- rep( (1:len[i]) + add.vec[i], each=len[i])
  }	
  BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
  
  if(!is.null(xvec)){
    testsymm <- abs(sum(skewpart(BlockDiag)))
    if(testsymm != 0) {
      warning("Correlation matrix is not computed to be exactly symmetric. Taking only the symmetric part.")
    }
  }
  return(list(BDiag = symmpart(BlockDiag), row.vec =row.vec, col.vec=col.vec))
}
