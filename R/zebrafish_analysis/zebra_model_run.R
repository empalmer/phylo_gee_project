library(tidyverse)
library(patchwork)
# Load data  --------------------------------------------------------------
# 10 percent
#test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_10p.rds")
# 20 percent
#test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_20p.rds")
# 30 percent
test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_30p.rds")
# flavo test: 
#test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_flavo_10p.rds")


test_data$descr
dat <- test_data$dat
group <- test_data$group
D <- test_data$D
asv_names <- test_data$asv_names
ps <- test_data$phy



# Run the model  ----------------------------------------------------------
# Load most recent model code
source(here::here("R","dm_cor_gee_clean.R"))
model_output <- dm_cor_gee(Y = dat$Y, X = group,
                           sample_id = dat$sampleID, 
                           ASV_id = asv_names,
                           distance_matrix = D, 
                           max_iter = 100, tol = .00001, 
                           gamma = .1, lambda = .01, 
                           save_beta = T)
# plot and save plots and output 
descr <- "F.3, gamma .1"
plots <- fun_diagnostics(model_output, descr)
plots
write_rds(list(descr = descr, results = model_output, plot = plots),
          file = here::here("Output",paste0("model_run",Sys.time(),".rds")))


beta_iteration_compare(model_output, descr)


prev_model_dat[[1]]
prev_model <- prev_model_dat$results

beta_model_compare(model_output, prev_model)


# Diagnostics: ------------------------------------------------------------

fun_diagnostics <- function(model_output, descr){
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


# Extra Diagnositics exploration ------------------------------------------

cor_phy <- exp(-2*model_output$rhos[[200]]*D)
heatmap(cor_phy,Colv = NA, Rowv = NA)




# residuals 
last_resids <- model_output$st_resid[[1]][[model_output$num_iter]]
hist(last_resids)
summary(last_resids)

# betas 
last_beta <- model_output$beta[[1]]
last_beta


beta1 <- matrix(last_beta, ncol = 2,byrow=T )

# GEE sum absolute value. 
data.frame(old = model_output$G, new = model_output$G_new, x = 1:length(model_output$G)) %>% 
  pivot_longer(cols = c(old, new)) %>% 
  ggplot(aes(x = x, y = value, color = name)) + geom_point() + geom_line() 

cbind(tax_table(ps),round(beta1, digits = 3))[,-c(1,7)]

test <- do.call(paste, c(data[my_cols], sep = ""))

tidyr::unite_(tax_table(ps), paste(colnames(tax_table(ps))[-1], collapse="_"), colnames(data))

hist(last_beta)
summary(last_beta)


