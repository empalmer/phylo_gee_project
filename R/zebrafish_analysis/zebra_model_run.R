
# Load data  --------------------------------------------------------------
# 10 percent
# test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_10p.rds")
# 20 percent
# test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_20p.rds")
# 30 percent
# test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_30p.rds")

test_data$descr
dat <- test_data$dat
group <- test_data$group
D <- test_data$D
asv_names <- test_data$asv_names



# Run the model  ----------------------------------------------------------
# Load model code
source(here::here("R","dm_cor_gee_clean.R"))
model_output <- dm_cor_gee(Y = dat$Y, X = group, sample_id = dat$sampleID, 
                    ASV_id = asv_names, distance_matrix = D, 
                    max_iter = 50)
# Save output
descr <- "Model run with filter of .1 and start rho of 10"
write_rds(list(descr, results = model_output),
          file = here::here("Output",paste0("model_run",Sys.time(),".rds")))







# Diagnostics: ------------------------------------------------------------

library(tidyverse)
diag_df <- data.frame(count = 1:model_output$num_iter, 
           phi  = model_output$phis, 
           rho  = model_output$rhos,
           omega = unname(model_output$omegas), 
           diff  = model_output$differences)

# phi 
diag_df %>% 
  ggplot(aes(x = count, y = phi)) + 
  geom_point()
summary(diag_df$phi)

# omega 
diag_df %>% 
  ggplot(aes(x = count, y = omega)) + 
  geom_point()

# rho 
diag_df %>% 
  ggplot(aes(x = count, y = rho)) + 
  geom_point()

# speed of convergence 
# diff 
diag_df %>% 
  ggplot(aes(x = count, y = diff)) + 
  geom_point()

# residuals 
last_resids <- model_output$st_resid[[model_output$num_iter]]
hist(last_resids)
summary(last_resids)

# betas 
last_beta <- model_output$betas[[1]][[model_output$num_iter]]
hist(last_beta)
summary(last_beta)


