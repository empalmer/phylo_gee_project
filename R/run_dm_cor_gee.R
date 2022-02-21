
# Can load practice data from https://github.com/lichen-lab/glmmTree in file glmmTree_data.R
source(here::here("R","glmmTree_sim_data.R"))
# Filtered to have 50 samples and 48 OTUs so computation goes quickly 



# default <- options()
# options(error = recover)
# options(default)

# load functions 
source(here::here("R","dm_cor_gee.R"))
source(here::here("R","dm_cor_gee_clean.R"))
# try on smaller dataset with 50 samples 
res1 <- dm_cor_gee2(Y = Y_data_small,X =  X_data_small, id = id_small, 
           distance_matrix = D_filtered)

# Try w/ geem
library(geeM)
# need to first make X the ip version so there are enough betas.
n = 50
X_intercept <- data.frame(int = rep(1, length(X_data_small)), 
                          X = X_data_small)
X_mat <- as.matrix(X_intercept, nrow = n)
# Start setting up the design matrix for GEE
Ip <- diag(43)
# For use in the beta loop
Ip_n <- bdiag(rep(list(Ip),n))
X <- kronecker(X_mat, Ip)
X2 <- kronecker(Ip, X_mat)
res_geeM <- geem(Y_data_small ~ X -1, id = id_small, init.beta = res_geepack$coefficients)
summary(res_geeM)

# now try the geepack package 
library(geepack)
res_geepack <- geeglm(Y_data_small ~ X -1, id = id_small)
res_geepack$coefficients



# Try on dataset with all samples but strongly filtered OTUs: 
# 100 samples n = 100
# one covariate - q =1
# Number of OTUs 
source(here::here("R","dm_cor_gee.R"))
all_samples_test <- dm_cor_gee(Y = Y_data_all_samples,X =  X_data_all_samples,
           id = id_all_samples, distance_matrix = D_filtered)



# Try on dataset with all samples and less filtered OTUs: 
# use 45% prevalence which gives 193 ASVs and 100 samples
source(here::here("R","dm_cor_gee.R"))
all_samples_test <- dm_cor_gee(Y = Y_data_45_all,X =  X_data_45_all,
                               id = id_45_all, distance_matrix = D_filtered_45)

