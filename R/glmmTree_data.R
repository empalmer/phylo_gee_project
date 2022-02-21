load("~/Desktop/phylo_gee_project/Data/glmmTree.Rdata")
library(tidyverse)
library(phyloseq)
# I think this is full all OTU data
otu.tab <- otu_table(data.obj$otu.tab, taxa_are_rows = T)
# I think this is the otu table aggregated to various taxa levels 
abund.list <- data.obj$abund.list
# Meta data - 39 variables on 528 samples 
meta.dat <- sample_data(data.obj$meta.dat)
# Phylogenetic tree 
tree <- phy_tree(data.obj$tree)
# I think this is the taxa classification? 
otu.name <- tax_table(data.obj$otu.name)
# also otu.name.full 
# Size factor? Looks like data were normalized using TSS 
# So this would be the total count in each sample. 
size.fator <- data.obj$size.factor



age_ps <- phyloseq(otu.tab, meta.dat, otu.name, tree)

# Arbitrary sample sums and family subsetting 
# Choose a good one 
#table(tax_table(age_ps)[,"Family"])
fam_ps <- subset_taxa(age_ps, Family=="Veillonellaceae")

# fam_ps <- prune_samples(sample_sums(fa_ps)>=1000, age_ps)
# fewer samples arbitrarily for test run
fam_ps <- prune_samples(sample_names(fam_ps)[1:100], fam_ps)
fam_ps
fam_ps <- prune_samples(!is.na(sample_data(fam_ps)$age), fam_ps)

fam_ps <- prune_taxa(apply(data.frame(otu_table(fam_ps, taxa_are_rows = T)),
                           1,
                           function(taxa) 
                             return(sum(taxa > 0) > .25 * phyloseq::nsamples(fam_ps))), 
                     fam_ps)

fam_ps

sum(otu_table(fam_ps, taxa_are_rows = T) == 0)/sum(otu_table(fam_ps, taxa_are_rows = T) >= 0)
# 37% of table entries are zero. 

# add pseudocount. Try 1. 
fam_ps_nonzero <- phyloseq::transform_sample_counts(fam_ps, function(x) x + 1)


norm_TSS <- function(ps, keep_prop = T){
  # keep as proportions or convert to counts per million?
  scale <- ifelse(keep_prop, 1, 1e6)
  # TSS function
  ps_normed <- phyloseq::transform_sample_counts(ps, function(x) x * scale / sum(x))
  return(ps_normed)
}

fam_ps_norm <- norm_TSS(fam_ps_nonzero)
head(otu_table(fam_ps_norm))

# Now build distance matrix: 
library(adephylo)
#library(ape)
# Saves as a dist object, need to convert to matrix. 
age_dist <- distTips(phy_tree(fam_ps_norm), method = "patristic")
D_age <- as.matrix(age_dist)

otu_table(fam_ps) %>% View()
sample_data(fam_ps) %>% View()


Y_df <- data.frame(t(otu_table(fam_ps_norm, taxa_are_rows = F)))%>%
  rownames_to_column("sampleID") %>%
  pivot_longer(-c(sampleID), names_to = "OTU_ID", values_to = "Y") 
  
Y <- Y_df %>% pull(Y)
id <- Y_df %>% pull(sampleID)
# use age as X 
X <- data.frame(sample_data(fam_ps_norm)) %>% 
  rownames_to_column("sampleID") %>% 
  select(sampleID, age) %>% 
  pull(age)

sam_data(fam_ps) %>% View()
taxa_names(fam_ps)
sample_names(fam_ps)

# load functions 
source(here::here("R","dm_cor_gee_clean.R"))
age_res1 <- dm_cor_gee2(Y = Y, X =  X, id = id, 
                   distance_matrix = D_age)
#n = 528
#p = 114
#q = 2