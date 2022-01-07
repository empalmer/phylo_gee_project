load("~/Desktop/phylo_gee_project/Data/glmmTree.Rdata")

# I think this is full all OTU data
otu.tab <- data.obj$otu.tab

# I think this is the otu table aggregated to various taxa levels 
abund.list <- data.obj$abund.list

# Meta data - 39 variables on 528 samples 
meta.dat <- data.obj$meta.dat

# Phylogenetic tree 
tree <- data.obj$tree

# I think this is the taxa classification? 
otu.name <- data.obj$otu.name
# also otu.name.full 

# Size factor? Looks like data were normalized using TSS 
# So this would be the total count in each sample. 
size.fator <- data.obj$size.factor

