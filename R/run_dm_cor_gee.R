
# Can load practice data from https://github.com/lichen-lab/glmmTree in file glmmTree_data.R
source('~/Desktop/ResearchFall2021/R/glmmTree_data.R')
# Filtered to have 50 samples and 48 OTUs so computation goes quickly 



# default <- options()
options(error = recover)
options(default)

dm_cor_gee(Y = Y_data,X =  X_data, id = id, 
           distance_matrix = D_filtered)
# dm_cor_gee(formula = Y ~ X, data = mb_data, id = sampleID)
