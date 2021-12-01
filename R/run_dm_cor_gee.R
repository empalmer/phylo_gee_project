
# Can load practice data from https://github.com/lichen-lab/glmmTree in file glmmTree_data.R
source(here::here("R","glmmTree_data.R"))
# Filtered to have 50 samples and 48 OTUs so computation goes quickly 



# default <- options()
# options(error = recover)
# options(default)


# try on smaller dataset with 50 samples 
dm_cor_gee(Y = Y_data_small,X =  X_data_small, id = id_small, 
           distance_matrix = D_filtered)


# Try on dataset with all samples but strongly filtered OTUs: 
dm_cor_gee(Y = Y_data_all_samples,X =  X_data_all_samples,
           id = id_all_samples, distance_matrix = D)

# Try on dataset with all samples and less filtered OTUs: 

