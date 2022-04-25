

# Packages ----------------------------------------------------------------

library(tidyverse)
library(phyloseq)

# Load data ---------------------------------------------------------------


# Use 10p filtered data: 
test_data <- readRDS("~/Desktop/phylo_gee_project/Data/test_data/test_1c_10p.rds")

phy <- test_data$phy


# Look at tree  -----------------------------------------------------------

plot_tree(phy,"treeonly", color="Class")


