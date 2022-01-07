
# Try using the distTips function 
# from adephylo 
#install.packages("adephylo")
library(adephylo)

# This tree is of class phylo so it should work? 
# t1 <- Sys.time()
# test_dist <- distTips(data.obj$tree,tips = data.obj$tree$tip.label[1:50] method = "patristic")
# diff <- Sys.time() - t1
# This naturally takes forever since this tree is for the unfiltered dataset 
# 

test_dist <- distTips(data.obj$tree, tips=data.obj$tree$tip.label[1:20], method = "patristic")



View(head(abund.list$Order))

data.obj$tree$tip.label[1:50]


library(ape)
# tests: 
tips <- data.obj$tree$tip.label[1:600]
pruned.tree<-drop.tip(data.obj$tree,data.obj$tree$tip.label[-match(tips, data.obj$tree$tip.label)])
write.tree(pruned.tree)

test_dist2 <- distTips(pruned.tree, method = "patristic")
test_dist2
