### Ohio respiratory data from geepack

data("ohio", package="geepack")
geem_ex1 <- geem(resp ~ age + smoke + age:smoke, id=id, data = ohio, family = binomial,
                    corstr = "m-dependent", Mv = 1 )

geem_ex1


mat <- matrix(c(1,2,0,0,
                2,1,2,0,
                0,2,1,2,
                0,0,2,1), nrow = 4)
geem_ex2 <- geem(resp ~ age + smoke + age:smoke, id=id, data = ohio, family = binomial,
                  corstr = "userdefined", corr.mat = mat )
geem_ex2
geem_ex1



## Microbiome dataset (from glmmtree)
library(tidyverse)
load("data1.rda")
# This has several included objects 
# z.tr: 778 OTU proportions for 100 samples in training set
# y.tr: Continous outcome for 100 samples in training set
# z.te: 778 OTU proportions for 200 samples in testing set
# y.te: Continous outcome for 200 samples in training set
# D: Patristic distance matrix among 778 OTUs

# Use the train data 
mb_data_otu <- z.tr
mb_data_sample <- data.frame(id = 1:200,y = y.tr)


# Convert into phyloseq for easier filtering
library(phyloseq)
otu <- otu_table(mb_data_otu, taxa_are_rows = F)
sam <- sample_data(mb_data_sample)
ps <- phyloseq(otu,sam)

colnames(mb_data_otu) <- paste0("OTU",1:778)
mb_data <- cbind(X = y.tr,
                sampleID = 1:200,
                data.frame(mb_data_otu)) %>% 
  pivot_longer(-c(X,sampleID), names_to = "OTU_ID", values_to = "Y") %>% 
  separate(OTU_ID,3, into = c("chr","OTU_ID")) %>% 
  mutate(OTU_ID = as.integer(OTU_ID))


geem_ind_mb <- geem(formula = Y ~ X, family = "binomial", id = sampleID,data = mb_data, 
                    corstr = "independence") 
# not done in 4 min... 

