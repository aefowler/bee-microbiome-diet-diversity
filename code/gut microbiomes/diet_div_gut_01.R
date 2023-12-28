# Pollen Diet Diversity does not Affect Gut Bacterial Communities or Melanization in a Social and Solitary Bee Species
# Alison E. Fowler · Quinn S. McFrederick · Lynn S. Adler

##### The effect of diet diversity on bee gut microbes
##### Script 1 
##### Loading and formatting community data 

# load libraries 
library(tidyverse)
library(plyr)

# import OTU table
lab_feat_table <- read.table("data/gut-microbe-diversity-data/diet_diversity_feature_table.txt", header = T, sep = "") 
# had to first remove hashtag in front of first col and change OTU ID to OTU_ID 
# this table has the OTUs as the first column and samples as the rest of the columns

# import metadata 
# # metadata from sequenced samples 
# metadata_partial <- read_csv("data/gut-microbe-diversity-data/lab_metadata_partial.csv") 
# # and metadata from the whole experiment (including samples that did not get sequenced )
# experiment_data <- read_csv("data/gut-microbe-diversity-data/lab_experiment_data.csv")
# # merge 
#lab_metadata <- merge(metadata_partial, experiment_data, by = "Bee_ID")

lab_metadata <- read_csv("data/gut-microbe-diversity-data/lab_metadata_full.csv")

# fix some variables  
#lab_metadata$Age <- revalue(lab_metadata$Age, c("adult" = "non-callow"))
lab_metadata$Age <- as.factor(lab_metadata$Age)
lab_metadata$Colony <- as.factor(lab_metadata$Colony)
lab_metadata$diet <- as.factor(lab_metadata$diet)
lab_metadata$plate <- as.factor(lab_metadata$plate)
lab_metadata$host <- as.factor(lab_metadata$host)
lab_metadata$No_pollen_spp <- as.factor(lab_metadata$No_pollen_spp)

# relevel diets  
lab_metadata$diet <- relevel(lab_metadata$diet, "Dandelion + Sumac + Hawthorn")
lab_metadata$diet <- relevel(lab_metadata$diet, "Dandelion + Sumac")
lab_metadata$diet <- relevel(lab_metadata$diet, "Sumac + Hawthorn")
lab_metadata$diet <- relevel(lab_metadata$diet, "Dandelion + Hawthorn")
lab_metadata$diet <- relevel(lab_metadata$diet, "Sumac")
lab_metadata$diet <- relevel(lab_metadata$diet, "Hawthorn")
lab_metadata$diet <- relevel(lab_metadata$diet, "Dandelion")

lab_metadata$diet_initials <- as.factor(lab_metadata$diet_initials)
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "D + S + H")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "D + S")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "S + H")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "D + H")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "S")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "H")
lab_metadata$diet_initials <- relevel(lab_metadata$diet_initials, "D")

# add index as rownames for metadata file 
lab_metadata_index <- column_to_rownames(lab_metadata, "index")
  
# format feature table 

# transpose OTU table 
lab_feat_table_t <- data.frame(t(lab_feat_table[-1]))
# make the column names the OTUs
colnames(lab_feat_table_t) <- lab_feat_table[,1]

########### for analyses comparing host species, use lab_feat_table_t #########

# separate the two host species: 
# separate the metadata 
bimp_metadata <- lab_metadata %>% dplyr::filter(host == "B. impatiens")
megarot_metadata <- lab_metadata %>% dplyr::filter(host == "M. rotundata")

# now for the OTU table
# make row names first column
lab_feat_table_t_index <- rownames_to_column(lab_feat_table_t, 'index')

# create df with just host and index columns 
host_index <- lab_metadata %>% dplyr::select(host, index)

# add this to the feature table so we can filter features by host
lab_feat_table_t_index_hosts <- merge(host_index, lab_feat_table_t_index, by = "index")

# separate by host species 
bimp_feat_table <- lab_feat_table_t_index_hosts %>% filter(host == "B. impatiens")
megarot_feat_table <- lab_feat_table_t_index_hosts %>% filter(host == "M. rotundata")

# remove host col & make index col rownames 
bimp_feat_table <- bimp_feat_table %>% dplyr::select(-host)
bimp_feat_table <- column_to_rownames(bimp_feat_table, 'index')

# remove host col & make index col rownames 
megarot_feat_table <- megarot_feat_table %>% dplyr::select(-host)
megarot_feat_table <- column_to_rownames(megarot_feat_table, 'index')

########## now they are ready to be used in analyses! ############ 

# do the same for diet
diet_index <- lab_metadata %>% dplyr::select(diet, index)

host_diet_index <- lab_metadata %>% dplyr::select(diet, host, index)

# add this to the feature table so we can filter features by host
lab_feat_table_t_index_diet_host <- merge(host_diet_index, lab_feat_table_t_index, by = "index")

bimp_feat_table_index_diet <- lab_feat_table_t_index_diet_host %>% filter(host == "B. impatiens")

megarot_feat_table_index_diet <- lab_feat_table_t_index_diet_host %>% filter(host == "M. rotundata")
