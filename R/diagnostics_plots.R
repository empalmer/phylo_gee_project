
# Diagnostics: ------------------------------------------------------------

fun_diagnostics <- function(model_output, descr = "no description provided"){
  diag_df <- data.frame(count = 1:model_output$num_iter, 
                        phi  = model_output$phi, 
                        rho  = model_output$rho,
                        omega = model_output$omega, 
                        diff  = model_output$diff, 
                        gee  = model_output$G)
  # phi 
  phi_plot <- diag_df %>% 
    ggplot(aes(x = count, y = phi)) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  
  
  # omega 
  omega_plot <- diag_df %>% 
    ggplot(aes(x = count, y = omega)) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  
  # rho 
  rho_plot <- diag_df %>% 
    ggplot(aes(x = count, y = rho)) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  
  # speed of convergence 
  # diff 
  diff_plot <- diag_df %>% 
    ggplot(aes(x = count, y = diff)) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  
  # GEE values 
  gee_plot <- diag_df %>% 
    ggplot(aes(x = count, y = gee)) + 
    geom_point() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  
  # Changing of betas
  beta_plot <- beta_iteration_compare(model_output, descr)
  
  # last residuals (standardized)
  #
  resid_plot <- ggplot(data.frame(x = model_output$resids)) + 
    geom_density(aes(x= x))
  
  # combine plots together. 
  patch <- (diff_plot + omega_plot + resid_plot) / 
    (rho_plot  + gee_plot + phi_plot) / beta_plot
  plot_annotation(descr)
  return(patch)
  
}

beta_iteration_compare <- function(model_output, descr){
  num_beta <- length(model_output$betas[[1]])/2
  betas <- data.frame(model_output$betas) 
  colnames(betas) <- 1:model_output$num_iter
  
  betas %>% 
    mutate(type = rep(c('int','x'), num_beta), 
           id = factor(1:(num_beta*2))) %>% 
    pivot_longer(1:model_output$num_iter, names_to = "iteration") %>% 
    ggplot(aes(x = as.numeric(iteration), y = value, group = id)) +
    geom_line() + 
    labs(x = "iteration", y = "beta")
}

beta_model_compare <- function(model1, model2){
  
  final_beta1 <- model1$betas[[model1$num_iter]]
  final_beta2 <- model2$betas[[model2$num_iter]]
  num_beta <- length(final_beta1)
  
  data_frame(final_beta1, final_beta2) %>% 
    mutate(id = factor(1:(num_beta))) %>% 
    pivot_longer(1:2) %>% 
    ggplot(aes(x = name, y = value, group = id)) + 
    geom_line() + 
    labs(x = "", y = "beta_value")
}
