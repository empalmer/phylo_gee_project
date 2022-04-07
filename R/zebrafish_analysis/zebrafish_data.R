
# package for tree stuff
require(ape)
require(tidyverse)
require(phyloseq)
# Use distTips function from adephylo 
#install.packages("adephylo")
require(adephylo)
# function to make a normal OTU table into 'long' format 
source(here::here("R","make_long_data.R"))

# Import files ------------------------------------------------------------

# Old tree 
#zebra_tree <- ape::read.tree(here::here("Data","Zebrafish",
#                                        "microbiome_data","alignedsequences.tre"))

# New tree 
zebra_tree <- ape::read.tree(here::here("Data","Zebrafish",
                                       "updated_tree","alignedsequences.tre"))

# This is from paired sequences Box folder 
zebra_phy <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/pre_made_phyloseq/full_ps.rds")

# convert tree tip labels to be lowercase
zebra_tree$tip.label <- tolower(zebra_tree$tip.label)

# add phylogenetic tree to phyloseq object
zebra_tree <- phy_tree(zebra_tree)
phy_tree(zebra_phy) <- zebra_tree

# Covariate Definitions -----------------------------------------------------
# setup covariates - use group - ie antibiotic or not 
# check with tom 
sample_data(zebra_phy)$group <- ifelse(
  sample_data(zebra_phy)$group == "NC" |
    sample_data(zebra_phy)$group == "ND",1,0)

# Start of filtering ------------------------------------------------------
## Filtering: 
# filter to only be last day: day = 32
day32 <- prune_samples(sample_data(zebra_phy)$day == 32, zebra_phy)
# 68 samples
day32


# Abundance filtering -----------------------------------------------------
# filter to for abundance 
# handy filtering functions
# remotes::install_github("vmikk/metagMisc")
require(metagMisc)
#original percent of data zeros 
sum(otu_table(day32) == 0)/(3895*68)


# _ 10 percent  -----------------------------------------------------------
# Include only taxa present in at least 10% samples
ps_filter_10 <- phyloseq_filter_prevalence(day32, prev.trh = 0.1, abund.trh = NULL) 
ps_filter_10
# percent of data zeros: 
sum(otu_table(ps_filter_10) == 0)/(nsamples(ps_filter_10)*ntaxa(ps_filter_10))

# Calculate distance matrix from filtered dataset
D_10 <- as.matrix(distTips(phy_tree(ps_filter_10), method = "patristic"))

# Calculate relative abundances from filtered samples
ps_filter_10 <- phyloseq::transform_sample_counts(ps_filter_10,
                                               function(x) x / sum(x))
# save ASV names 
asv_names_10 <- taxa_names(ps_filter_10)
# Objects needed are D_10 and filter_10 and names_10

# ___ ND vs NC  --------------------------------
# Use ND/NC use as the sole covariate, saved as "group"
group_10 <- sample_data(ps_filter_10)$group
# convert Y to long format 
dat_10 <- make_long_data(group_10, ps_filter_10)

descr = "filtered at 10p, nc/nd covariate"
write_rds(list(descr = descr, 
               D = D_10, 
               group = group_10, 
               dat = dat_10, 
               asv_names = asv_names_10), 
          here::here("Data", "test_data", "test_1c_10p.rds"))



# _ 20 percent ------------------------------------------------------------
# Include only taxa present in at least 20% samples
ps_filter_20 <- phyloseq_filter_prevalence(day32,
                                           prev.trh = 0.2,
                                           abund.trh = NULL) 
ps_filter_20
# percent of data zeros: 
sum(otu_table(ps_filter_20) == 0)/(nsamples(ps_filter_20)*ntaxa(ps_filter_20))

# Calculate distance matrix from filtered dataset
D_20 <- as.matrix(distTips(phy_tree(ps_filter_20), method = "patristic"))

# Calculate relative abundances from filtered samples
ps_filter_20 <- phyloseq::transform_sample_counts(ps_filter_20,
                                                  function(x) x / sum(x))
# save ASV names 
asv_names_20 <- taxa_names(ps_filter_20)

# ___ ND vs NC  --------------------------------
# Use ND/NC use as the sole covariate, saved as "group"
group_20 <- sample_data(ps_filter_20)$group
# convert Y to long format 
dat_20 <- make_long_data(group_20, ps_filter_20)

descr = "filtered at 20p, nc/nd covariate"
write_rds(list(descr = descr, 
               D = D_20, 
               group = group_20, 
               dat = dat_20, 
               asv_names = asv_names_20), 
          here::here("Data", "test_data", "test_1c_20p.rds"))


# _ 30 percent ------------------------------------------------------------

# Include only taxa present in at least 30% samples
ps_filter_30 <- phyloseq_filter_prevalence(day32,
                                           prev.trh = 0.3,
                                           abund.trh = NULL) 
ps_filter_30
# percent of data zeros: 
sum(otu_table(ps_filter_30) == 0)/(nsamples(ps_filter_30)*ntaxa(ps_filter_30))

# Calculate distance matrix from filtered dataset
D_30 <- as.matrix(distTips(phy_tree(ps_filter_30), method = "patristic"))

# Calculate relative abundances from filtered samples
ps_filter_30 <- phyloseq::transform_sample_counts(ps_filter_30,
                                                  function(x) x / sum(x))
# save ASV names 
asv_names_30 <- taxa_names(ps_filter_30)

# ___ ND vs NC  --------------------------------
# Use ND/NC use as the sole covariate, saved as "group"
group_30 <- sample_data(ps_filter_30)$group
# convert Y to long format 
dat_30 <- make_long_data(group_30, ps_filter_30)

descr = "filtered at 30p, nc/nd covariate"
write_rds(list(descr = descr, 
               D = D_30, 
               group = group_30, 
               dat = dat_30, 
               asv_names = asv_names_30), 
          here::here("Data", "test_data", "test_1c_30p.rds"))
