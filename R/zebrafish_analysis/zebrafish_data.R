
# package for tree stuff
require(ape)
require(tidyverse)
require(phyloseq)

# Old tree 
#zebra_tree <- ape::read.tree(here::here("Data","Zebrafish",
#                                        "microbiome_data","alignedsequences.tre"))

# New tree 
zebra_tree <- ape::read.tree(here::here("Data","Zebrafish",
                                       "updated_tree","alignedsequences.tre"))


# phyloseq object already created 
# tax_final <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/microbiome_data/tax_final.rds")
# sequence_table <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/microbiome_data/sequence_table.rds")
# tax_table <- tax_table(tax_final)
# otu_table <- otu_table(sequence_table, taxa_are_rows = F)

zebra_phy <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/pre_made_phyloseq/full_ps.rds")

# convert tree tip labels to be lowercase
zebra_tree$tip.label <- tolower(zebra_tree$tip.label)

# add phylogenetic tree to phyloseq object
zebra_tree <- phy_tree(zebra_tree)
phy_tree(zebra_phy) <- zebra_tree



# Different metadata files 
# Not needed - already loaded in the full_data.rds file. 
# What to use as covariates? 
# pcap_np_metadata <- read_excel("Data/Zebrafish/metadata/pcap_np_metadata.xlsx")
# meta1 <- sample_data(pcap_np_metadata)
# zebra_phyloseq <- phyloseq(tax_table, otu_table, meta1)
# sample_names(meta1)
# 
# gaulke_metabolomics_metadata <- read_excel("Data/Zebrafish/metadata/gaulke_metabolomics_metadata.xlsx")
# meta2 <- sample_data(gaulke_metabolomics_metadata)
# zebra_phyloseq <- phyloseq(tax_table, otu_table, meta2)
# sample_names(meta2)
# 
# pcap_metabolomics_metadata <- read_excel("Data/Zebrafish/metadata/pcap_metabolomics_metadata.xlsx")



## Filtering: 
# filter to only be last day: day = 32

day32 <- prune_samples(sample_data(zebra_phy)$day == 32,zebra_phy)
# 68 samples
day32


# filter to for abundance 
# handy filtering functions
# remotes::install_github("vmikk/metagMisc")
require(metagMisc)

#original percent of data zeros 
sum(otu_table(day32) == 0)/(3895*68)
# Include only taxa with more than 1 reads (on average) in at least 50% samples
filtered_zebras <- phyloseq_filter_prevalence(day32, prev.trh = 0.3, abund.trh = NULL) 
filtered_zebras

# percent of data zeros: 
sum(otu_table(filtered_zebras) == 0)/(39*68)



# tests: 
# tips <- zebra_tree$tip.label[1:100]
# pruned.tree<-drop.tip(zebra_tree,zebra_tree$tip.label[-match(tips, zebra_tree$tip.label)])
# write.tree(pruned.tree)

# Try using the distTips function 
# from adephylo 
#install.packages("adephylo")
require(adephylo)

filtered_zebra_tree <- phy_tree(filtered_zebras)
zebra_dist <- distTips(filtered_zebra_tree, method = "patristic")



rm(day32)
rm(filtered_zebra_tree)
rm(zebra_tree)
rm(zebra_phy)


# setup covariates - use group - ie antibiotic or not 
# check with tom 
sample_data(filtered_zebras)$group <- 
  ifelse(sample_data(filtered_zebras)$group == "NC" | 
         sample_data(filtered_zebras)$group == "ND",
         1,0)

zebrafish_ps <- phyloseq::transform_sample_counts(filtered_zebras,
                                                  function(x) x / sum(x))

rm(filtered_zebras)