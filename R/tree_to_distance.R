
# Try using the distTips function 
# from adephylo 
#install.packages("adephylo")
library(adephylo)

# This tree is of class phylo so it should work? 
t1 <- Sys.time()
test_dist <- distTips(data.obj$tree, method = "patristic")
diff <- Sys.time() - t1