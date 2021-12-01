# default <- options()
options(error = recover)
options(default)

dm_cor_gee(Y = Y_data,X =  X_data, id = id, 
           distance_matrix = D_filtered)
# dm_cor_gee(formula = Y ~ X, data = mb_data, id = sampleID)
