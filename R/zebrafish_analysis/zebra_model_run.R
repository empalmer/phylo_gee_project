 
# load data 
source(here::here("R","zebrafish_analysis","zebrafish_data.R"))

# function to make a normal OTU table into 'long' format 
source(here::here("R","make_long_data.R"))

zebra2 <- prune_samples((rowSums(otu_table(zebrafish_ps) == 0)/19) < .4, zebrafish_ps)
zebra2

sum(otu_table(zebra2) == 0)/(19*47)

# Use antibiotic use as the sole covariate. 
group <- sample_data(zebra2)$group
# convert Y to long format 
dat <- make_long_data(group, zebra2)
# convert distance into matrix 

zebra_dist
D <- as.matrix(zebra_dist)

# look at histogram of proportions - still heavily skewed
hist(dat$Y)




# model! 
source(here::here("R","dm_cor_gee_clean.R"))
zebra_res <- dm_cor_gee2(Y = dat$Y,X =  group, id = dat$sampleID, 
                    distance_matrix = D)
