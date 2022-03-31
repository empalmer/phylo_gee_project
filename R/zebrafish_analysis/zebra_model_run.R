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

start.time <- Sys.time()
zebra_res <- dm_cor_gee(Y = dat$Y, X = group, sample_id = dat$sampleID, 
                    ASV_id = taxa_names(zebrafish_ps), distance_matrix = D)
end.time <- Sys.time()
time.taken <- end.time - start.time




