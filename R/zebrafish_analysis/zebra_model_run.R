library(tidyverse)
# Load data  --------------------------------------------------------------
# 10 percent
# test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_10p.rds")
# 20 percent
test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_20p.rds")
# 30 percent
#test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_30p.rds")

test_data$descr
dat <- test_data$dat
group <- test_data$group
D <- test_data$D
asv_names <- test_data$asv_names
ps <- test_data$ps



# Run the model  ----------------------------------------------------------
# Load most recent model code
source(here::here("R","dm_cor_gee_clean.R"))
model_output <- dm_cor_gee(Y = dat$Y, X = group, sample_id = dat$sampleID, 
                    ASV_id = asv_names, distance_matrix = D, 
                    max_iter = 30, tol = .0001)
# Save output
descr <- "F.3, lambda max(eigen)/1000, 30 rep"
write_rds(list(descr, results = model_output),
          file = here::here("Output",paste0("model_run",Sys.time(),".rds")))




b1 <- model_output1$beta
b2 <- model_output2$beta
b1
b2

plot(b1[[1]],b2[[1]])

# Diagnostics: ------------------------------------------------------------
diag_df <- data.frame(count = 1:model_output$num_iter, 
           phi  = model_output$phis, 
           rho  = model_output$rhos,
           omega = unname(model_output$omegas), 
           diff  = model_output$differences)

# GEE sum absolute value. 
data.frame(old = model_output$G, new = model_output$G_new, x = 1:length(model_output$G)) %>% 
  pivot_longer(cols = c(old, new)) %>% 
  ggplot(aes(x = x, y = value, color = name)) + geom_point() + geom_line()


# phi 
diag_df %>% 
  ggplot(aes(x = count, y = phi)) + 
  geom_point()
summary(diag_df$phi)

# omega 
diag_df %>% 
  ggplot(aes(x = count, y = omega)) + 
  geom_point()


model_output

# rho 
diag_df %>% 
  ggplot(aes(x = count, y = rho)) + 
  geom_point()




cor_phy <- exp(-2*model_output$rhos[[200]]*D)
heatmap(cor_phy,Colv = NA, Rowv = NA)


# speed of convergence 
# diff 
diag_df %>% 
  ggplot(aes(x = count, y = diff)) + 
  geom_point()

# residuals 
last_resids <- model_output$st_resid[[1]][[model_output$num_iter]]
hist(last_resids)
summary(last_resids)

# betas 
last_beta <- model_output$beta[[1]]
last_beta


beta1 <- matrix(last_beta, ncol = 2,byrow=T )


cbind(tax_table(ps),round(beta1, digits = 3))[,-c(1,7)]

test <- do.call(paste, c(data[my_cols], sep = ""))

tidyr::unite_(tax_table(ps), paste(colnames(tax_table(ps))[-1], collapse="_"), colnames(data))

hist(last_beta)
summary(last_beta)


