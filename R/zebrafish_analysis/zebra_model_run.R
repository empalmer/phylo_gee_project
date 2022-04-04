 # load data 
source(here::here("R","zebrafish_analysis","zebrafish_data.R"))
# function to make a normal OTU table into 'long' format 
source(here::here("R","make_long_data.R"))

# Use antibiotic use as the sole covariate. 
group <- sample_data(zebrafish_ps)$group
# convert Y to long format 
dat <- make_long_data(group, zebrafish_ps)
# convert distance into matrix 
taxa_names(zebrafish_ps)

D <- as.matrix(zebra_dist)

# look at histogram of proportions - still heavily skewed
hist(dat$Y)


# model! 
source(here::here("R","dm_cor_gee_clean.R"))
zebra_res <- dm_cor_gee(Y = dat$Y, X = group, sample_id = dat$sampleID, 
                    ASV_id = taxa_names(zebrafish_ps), distance_matrix = D)


# Diagnostics: 
library(tidyverse)
diag_df <- data.frame(count = 1:zebra_res$num_iter, 
           phi  = zebra_res$phis, 
           rho  = zebra_res$rhos,
           omega = unname(zebra_res$omegas), 
           diff  = zebra_res$differences)

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


