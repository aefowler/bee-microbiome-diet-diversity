##### The effect of diet diversity on bee gut microbes
##### Script 3
##### Comparing beta diversity of microbiomes between hosts

# load libraries and other scripts 
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(ape)
library(ggpubr)
library(vegan)
library(lattice)
library(permute)
library(ggpubr)

# create a vector with the host group information 
host_groups <- factor(lab_feat_table_t_index_hosts$host)

# calculate bray-curtis distance matrix 
dist <- vegdist(lab_feat_table_t, method="bray")

# simper 
host_simper <- with(lab_metadata_index, simper(lab_feat_table_t, host))
summary(host_simper)

# betadisper model with bray curtis distance matrix 
host_dist_mod <- betadisper(dist, host_groups)

anova(host_dist_mod)
# Response: Distances
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      1 0.15102 0.151019  11.345 0.0009602 ***
#   Residuals 151 2.01008 0.013312  

permutest(host_dist_mod)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups      1 0.15102 0.151019 11.345    999  0.001 ***
#   Residuals 151 2.01008 0.013312        

# plot 
betadisper_host_plot1 <- plot(host_dist_mod, col = c("orchid4", "deepskyblue"))

# adonis 
adonis2(dist ~ host_groups, perm = 10000) 
# adonis2(formula = dist ~ host_groups, permutations = 10000)
# Df SumOfSqs      R2     F    Pr(>F)    
# host_groups   1    9.992 0.21733 41.93 9.999e-05 ***
#   Residual    151   35.985 0.78267                    
# Total       152   45.977 1.00000          

# PHYLOSEQ

# import data 
physeq <- qza_to_phyloseq(
  features="data/gut-microbe-diversity-data/lab_filttable3.qza",
  taxonomy="data/gut-microbe-diversity-data/lab_16STaxonomy.qza",
  metadata = "data/gut-microbe-diversity-data/lab_metadata_full_phyloseq.tsv"
)

physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334 taxa and 153 samples ]
# sample_data() Sample Data:       [ 153 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 334 taxa by 7 taxonomic ranks ]

# bar plot
plot_bar(physeq, fill = "Family")

# make a phylogeny 
random_tree <- rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))

# merge with phylseq object 
physeq1 <- merge_phyloseq(physeq, random_tree)
physeq1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 334 taxa and 153 samples ]
# sample_data() Sample Data:       [ 153 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 334 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 334 tips and 333 internal nodes ]

#### Change NAs to be named as unclassified to the lowest possible rank 

for(i in 1:nrow(physeq1@tax_table)) { # go through each row  
  #  print(colnames(temp)[i]) # for debugging 
  for (j in 1:ncol(physeq1@tax_table)) { # then go through each column 
    #  print(temp[i,j]) # for debugging
    if (is.na(physeq1@tax_table[i,j])) { # if it is NA 
      if(startsWith(physeq1@tax_table[i,j-1], # AND if cell in previous column starts with "unclassified" 
                    "Unclassified")){
        v <- strsplit(physeq1@tax_table[i,j-1], split = " ") # then split the previous cell contents into two strings 
        physeq1@tax_table[i,j] <- paste("Unclassified", unlist(v[[1]][2])) # assign as "unclassified" + last taxanomic level 
        #  print(paste("print this", unlist(v[[1]][2]))) # for debugging
      } else {physeq1@tax_table[i,j] <- paste("Unclassified", physeq1@tax_table[i,j-1]) # if NA and not unclassified in previous column, assign as unclassified + taxonomic level in previous cell 
      }
    }
  }
}

tax_table <- as.data.frame(physeq1@tax_table)
tax_table <- rownames_to_column(tax_table, "OTU")
genus_table <- tax_table %>% select("OTU", "Genus")
species_table <- tax_table %>% select("OTU", "Genus", "Species")

# plot heatmap 
plot_heatmap(physeq1, sample.order = "host", low = "darkblue", high = "gold")

### Create data frome from phyloseq object that we will then use for plotting 
host_genera <- physeq1 %>%
 # tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Abundance) 

# further classify with BLAST and then rename 

host_genera <- within(host_genera, Genus[OTU == "4a8e08ef613b9fc1d1523112af5aa01e"] <- "Bombus (host)")

host_genera <- within(host_genera, Genus[OTU == "f72afeb78337970564ef6737a753a359"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "f72afeb78337970564ef6737a753a359"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "d9746fb0b36e4afb76dfec7e3e8925d3"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "e09bc90831bb219468010611697a8a90"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "c9fd02cbcdb891506aa96e92e6bdbbb3"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "ca310449cd52d0005b5281d9657d3963"] <- "Pantoea")
host_genera <- within(host_genera, Genus[OTU == "c89b42ea6dc37965a4a48a7361a313cb"] <- "Pantoea")

host_genera <- within(host_genera, Genus[OTU == "35b0bc6e20cc31eb6108a12530d1800b"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "a220dbde7245dfb7899ec7fe484525ed"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "ff2293ea3b41482e6d11f94b080b7ae7"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "f271e6e6d977397e39f1d55be04bf3ec"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "099f235697942b6e9afd171a3200e7b7"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "4a98edcaa9b90656232fd99aa8bd0c51"] <- "Unclassified Enterobacteriaceae")

host_genera <- within(host_genera, Genus[OTU == "bc7a67162e397d0b74121e9fe9602d08"] <- "Klebsiella") ## 
host_genera <- within(host_genera, Genus[OTU == "df208d4225db4d84097d4635026c0ada"] <- "Klebsiella") ## 
host_genera <- within(host_genera, Genus[OTU == "646be12e02a43f09ab0c5065ade76df3"] <- "Klebsiella") ##
host_genera <- within(host_genera, Genus[OTU == "5c80ac19f5637f31515c5ff3f1c59307"] <- "Klebsiella")
host_genera <- within(host_genera, Genus[OTU == "e47d852d8b0df2936195622ef95be11b"] <- "Klebsiella")
host_genera <- within(host_genera, Genus[OTU == "fe041460dc4756ea5ff0edf699e6c7e3"] <- "Klebsiella")
host_genera <- within(host_genera, Genus[OTU == "8f28a8b45e0c7e53c1ea4138de6d731b"] <- "Klebsiella")

host_genera <- within(host_genera, Genus[OTU == "29977f6c24ddd107d3db22f59a133da9"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "4353407195b1c7f4ee257f8b88fd532d"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "384112dd093eba5a3b56a44acd76ec02"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "b67c138c5ec29bb1a13dcefbb476b2d9"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "16ac73b8c3dacc6564afb87b5497bd35"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "be9cf2697f7bbf30993c3ba0129c9109"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "88a0786cb772eb069cbb66edd8419b0d"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "0787b72217e5d87504b8489206b45fe0"] <- "Rosenbergiella")
host_genera <- within(host_genera, Genus[OTU == "0f55fdd781c613fc15be4998b47856b3"] <- "Rosenbergiella")

host_genera <- within(host_genera, Genus[OTU == "28bd1c3ab72d270aa0b6d52795c5bae5"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "ed4e41ef76195901b7c71c1aa10b899a"] <- "Unclassified Enterobacteriaceae")
host_genera <- within(host_genera, Genus[OTU == "bdb4997736a918b100d59f0cdc5bcf6c"] <- "Unclassified Enterobacteriaceae")
# host_genera <- within(host_genera, Genus[OTU == "22c4f5e0b68d3359a46ccbf60c379b43"] <- "Unclassified Enterobacteriaceae") # leave as is 

host_genera$Genus <- as.factor(host_genera$Genus)
levels(host_genera$Genus)[levels(host_genera$Genus) == "Candidatus_Schmidhempelia" ] <- "Candidatus Schmidhempelia"
levels(host_genera$Genus)[levels(host_genera$Genus) == "Unclassified d__Bacteria" ] <- "Unclassified"

# reorder treatments 
host_genera$diet_initials <- as.factor(host_genera$diet_initials)
host_genera$diet_initials <- relevel(host_genera$diet_initials, "D + S + H")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "D + S")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "S + H")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "D + H")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "S")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "H")
host_genera$diet_initials <- relevel(host_genera$diet_initials, "D")

bimp_genera <- host_genera %>% filter(host == "B. impatiens")
mrot_genera <- host_genera %>% filter(host == "M. rotundata")

# relative abundance plots 

b.taxa.plot<-
  ggplot(bimp_genera, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_grid(.~diet_initials, scales = "free", space = "free") +
  geom_bar(stat = "identity") +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        plot.title = element_text(size = 25),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'black')) + 
  scale_fill_manual(values = c("firebrick", "lightpink", "lightskyblue", "mediumseagreen",
                               "red", "hotpink3", "olivedrab3", "grey", "gold", "maroon4",
                               "dodgerblue", "tan", "turquoise3", "royalblue3", "orchid")) +
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, nrow = 3, byrow = T)) +
  ylab("Relative Abundance (>1%)\n") +
  ggtitle(expression(paste(italic("B. impatiens ")))) 

b.taxa.plot

m.taxa.plot<-
  ggplot(mrot_genera, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_grid(.~diet_initials, 
             scales = "free", 
             space = "free", 
             #ncol = 3
  ) +
  geom_bar(stat = "identity") +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        plot.title = element_text(size = 25),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'black')) + 
  scale_fill_manual(values = c("firebrick", "lightpink", "olivedrab3", "forestgreen", "lightslateblue", 
                               "gold", "maroon4", "sienna2", "royalblue3", "orchid", "grey")) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (>1%)\n") +
  ggtitle(expression(paste(italic("M. rotundata ")))) + 
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

m.taxa.plot

taxa_plots <- ggarrange(b.taxa.plot, m.taxa.plot,
                        ncol = 1, nrow = 2, 
                        labels = c("A", "B"))
taxa_plots

ggsave(
  "fig3_7-2.jpg",
  plot = taxa_plots,
  width = 3000,
  height =2000,
  units = "px",
  dpi = 300
)

### DISTANCE METRICS RESULTS #### 
host.bray <- distance(physeq1, method = "bray")
host.wunifrac <- distance(physeq1, method = "wunifrac")
host.unifrac <- distance(physeq1, method = "unifrac")

host.df <- data.frame(sample_data(physeq1))

# betadisper bray
host.bray.betadisper <- betadisper(host.bray, host.df$host)
anova(host.bray.betadisper)

# betadisper weighted unifrac
host.wuni.betadisper <- betadisper(host.wunifrac, host.df$host)
anova(host.wuni.betadisper)

# betadisper unweighted unifrac
host.uni.betadisper <- betadisper(host.unifrac, host.df$host)
anova(host.uni.betadisper)

# adonis 
adonis2(host.bray ~host, data = host.df) 
adonis2(host.wunifrac ~host, data = host.df) 
adonis2(host.unifrac ~host, data = host.df) 

#### MAKE ORDINATION PLOTS ##### 

# ordinate using PCoA and bray curtis 
host_ord <- ordinate(physeq1, method = "PCoA", distance = "bray")
host_ord_plot <- plot_ordination(physeq1, host_ord, color = "host")
host_bc <- host_ord_plot + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95,
               lwd = 1.75) + 
  geom_point(size = 4.5) +
  theme_classic() + 
  theme(legend.position = c(.18,.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 22,
                                   face = "italic"),
        axis.title = element_text(size = 22, color = "black"),
        axis.text = element_text(size = 20,
                                 color = "black")) + 
  scale_color_manual(values = c("orchid4", "deepskyblue")) 

# ordinate using PCoA and weighted unifrac 
host_ord1 <- ordinate(physeq1, method = "PCoA", distance = "wunifrac")
host_ord_plot1 <- plot_ordination(physeq1, host_ord1, color = "host")
host_wuni <- host_ord_plot1 + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  ggtitle("Weighted UniFrac") + 
  theme_classic() + 
  scale_color_manual(values = c("orchid4", "deepskyblue")) + 
  theme(legend.position = "none")

# ordinate using PCoA and unweighted unifrac 
host_ord2 <- ordinate(physeq1, method = "PCoA", distance = "unifrac")
host_ord_plot2 <- plot_ordination(physeq1, host_ord2, color = "host")
host_uni <- host_ord_plot2 + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  ggtitle("Unweighted UniFrac") + 
  theme_classic() + 
  scale_color_manual(values = c("orchid4", "deepskyblue")) + 
  theme(legend.position = "none")


host_pcoas <- ggarrange(host_bc, host_wuni, host_uni,
                        nrow = 1, ncol = 3)
host_pcoas
