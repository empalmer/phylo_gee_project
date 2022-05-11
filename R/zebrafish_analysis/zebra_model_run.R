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
                           gamma = .05, lambda = .01, 
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


