 
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
                    ASV_id = taxa_names(zebra2), distance_matrix = D)






test_x <- matrix(c(1,1,1,2,3,3), nrow = 3)

test_x <- c(1,2,3,4)
model.matrix(~x, data.frame(x = test_x))

beta_mat <- matrix(c(1,1,1,3,4,5), nrow = 3)
beta_mat

as.numeric(tcrossprod(test_x, beta_mat))


unlist(map(split(test_x,rep(1:nrow(test_x))), ~as.numeric(tcrossprod(.x, beta_mat))))


