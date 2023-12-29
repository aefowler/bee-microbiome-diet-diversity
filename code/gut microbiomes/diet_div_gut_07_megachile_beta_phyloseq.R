##### The effect of diet diversity on bee gut microbes
##### Script 9 
##### Microbiomes in M. rotundata using PHYLOSEQ

# libraries 
library(tidyverse)
library(phyloseq)
library(ape)
library(qiime2R)
library(vegan)
library(plyr)

m.physeq <- qza_to_phyloseq(
  features="data/gut-microbe-diversity-data/lab_megarot_filttable.qza",
  taxonomy="data/gut-microbe-diversity-data/lab_16STaxonomy.qza",
  metadata = "data/gut-microbe-diversity-data/lab_megarot_metadata_phyloseq.tsv"
)

# create tree 
m_tree <- rtree(ntaxa(m.physeq), rooted=TRUE, tip.label=taxa_names(m.physeq))

# combine 
m.physeq1 <- merge_phyloseq(m.physeq, m_tree)
m.physeq1

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 163 taxa and 68 samples ]
# sample_data() Sample Data:       [ 68 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 163 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 163 tips and 162 internal nodes ]

### DISTANCE MATRICES RESULTS 

# Compare outputs with different distance metrics 
m.df <- data.frame(sample_data(m.physeq1))
m.df$No_pollen_spp <- as.factor(m.df$No_pollen_spp)
m.df$diet <- as.factor(m.df$diet)

# weighted unifrac 
m.wunifrac <- distance(m.physeq1, method = "wunifrac")

# Betadisper with diet comp 
m.wunifrac.betadisper <- betadisper(m.wunifrac, m.df$diet)
anova(m.wunifrac.betadisper)

# Betadisper with diet div 
m_div.wunifrac.betadisper <- betadisper(m.wunifrac, m.df$No_pollen_spp)
anova(m_div.wunifrac.betadisper)

# adonis 
adonis2(m.wunifrac ~ diet, data = m.df) 
adonis2(m.wunifrac ~ No_pollen_spp, data = m.df) 

# bray curtis
m.bray <- distance(m.physeq1, method = "bray")

# Betadisper with diet comp 
m.bray.betadisper <- betadisper(m.bray, m.df$diet)
anova(m.bray.betadisper)

# Betadisper with diet div 
m.div_bray.betadisper <- betadisper(m.bray, m.df$No_pollen_spp)
anova(m.div_bray.betadisper)

# adonis 
m.perms.n <- with(m.df, how(nperm = 1000, blocks = Emerge_start_date))
adonis2(m.bray ~ diet, permutations = m.perms.n, data = m.df, pairwise = TRUE)
adonis2(m.bray ~ No_pollen_spp, permutations = m.perms.n, data = m.df, pairwise = TRUE)

# unweighted unifrac
m.uni <- distance(m.physeq1, method = "unifrac")

# Betadisper with diet comp 
m.uni.betadisper <- betadisper(m.uni, m.df$diet)
anova(m.uni.betadisper)

# Betadisper with diet div 
m.div_uni.betadisper <- betadisper(m.uni, m.df$No_pollen_spp)
anova(m.div_uni.betadisper)

# adonis 
adonis2(m.uni ~ diet, data = m.df) 
adonis2(m.uni ~ No_pollen_spp, data = m.df) 

### ORDINATION PLOTS ####
m.ord <- ordinate(m.physeq1, method = "PCoA", distance = "bray")
m.ord2 <- ordinate(m.physeq1, method = "PCoA", distance = "wunifrac")
m.ord3 <- ordinate(m.physeq1, method = "PCoA", distance = "unifrac")

m.physeq1@sam_data$No_pollen_spp <- as.factor(m.physeq1@sam_data$No_pollen_spp)

### diet comp ###
# bray 
mega_bc <- plot_ordination(m.physeq1, m.ord, color = "diet")
mega_bc <- mega_bc + stat_ellipse(geom = "polygon",
                          aes(color = diet),
                          alpha = 0,
                          level = 0.95,
                          lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
  ggtitle("Bray-Curtis") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# weighted unifrac
wuni_mega <- plot_ordination(m.physeq1, m.ord2, color = "diet")
wuni_mega <- wuni_mega + stat_ellipse(geom = "polygon",
                          aes(color = diet),
                          alpha = 0,
                          level = 0.95,
                          lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
  ggtitle("Weighted UniFrac") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# unweighted unifrac
uni_mega <- plot_ordination(m.physeq1, m.ord3, color = "diet")
uni_mega <- uni_mega + stat_ellipse(geom = "polygon",
                           aes(color = diet),
                           alpha = 0,
                           level = 0.95,
                           lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
  ggtitle("Unweighted UniFrac") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

### diet div ###
# bray 
mega_div_bc <- plot_ordination(m.physeq1, m.ord, color = "No_pollen_spp")
mega_div_bc <- mega_div_bc + stat_ellipse(geom = "polygon",
                                  aes(color = No_pollen_spp),
                                  alpha = 0,
                                  level = 0.95,
                                  lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
#  ggtitle("Bray-Curtis") +
  theme(legend.position = "none")


# weighted unifrac
wuni_mega_div <- plot_ordination(m.physeq1, m.ord2, color = "No_pollen_spp")
wuni_mega_div <- wuni_mega_div + stat_ellipse(geom = "polygon",
                                      aes(color = No_pollen_spp),
                                      alpha = 0,
                                      level = 0.95,
                                      lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
#  ggtitle("Weighted UniFrac") +
  theme(legend.position = "none")

# unweighted unifrac
uni_mega_div <- plot_ordination(m.physeq1, m.ord3, color = "No_pollen_spp")
uni_mega_div <- uni_mega_div + stat_ellipse(geom = "polygon",
                                    aes(color = No_pollen_spp),
                                    alpha = 0,
                                    level = 0.95,
                                    lwd = 0.7) + 
  geom_point(size = 2) +
  theme_classic() + 
#  ggtitle("Unweighted UniFrac") +
  theme(legend.position = "none")

mega_pcoas <- ggarrange(mega_bc, wuni_mega, uni_mega,
                        mega_div_bc, wuni_mega_div, uni_mega_div,
                        nrow = 2, ncol = 3)
mega_pcoas

# BAR PLOTS 
##### Modify the physeq1 object to rename unclassified taxa 
for(i in 1:nrow(m.physeq1@tax_table)) { # go through each row  
  for (j in 1:ncol(m.physeq1@tax_table)) { # then go through each column 
    if (is.na(m.physeq1@tax_table[i,j])) { # if it is NA 
      if(startsWith(m.physeq1@tax_table[i,j-1], # AND if cell in previous column starts with "unclassified" 
                    "Unclassified")){
        v <- strsplit(m.physeq1@tax_table[i,j-1], split = " ") # then split the previous cell contents into two strings 
        m.physeq1@tax_table[i,j] <- paste("Unclassified", unlist(v[[1]][2])) # assign as "unclassified" + last taxanomic level 
      } else {m.physeq1@tax_table[i,j] <- paste("Unclassified", m.physeq1@tax_table[i,j-1]) # if NA and not unclassified in previous column, assign as unclassified + taxonomic level in previous cell 
      }
    }
  }
}

m.tax_table <- as.data.frame(m.physeq1@tax_table)
m.tax_table <- rownames_to_column(m.tax_table, "OTU")
m.genus_table <- m.tax_table %>% select("OTU", "Genus")
m.species_table <- m.tax_table %>% select("OTU", "Genus", "Species")

m_genus_plot <- m.physeq1 %>%
 # tax_glom(taxrank = "Genus") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                         
  arrange(Genus) 

m_genus_plot$diet <- as.factor(m_genus_plot$diet)
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Dandelion + Sumac + Hawthorn")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Dandelion + Sumac")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Sumac + Hawthorn")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Dandelion + Hawthorn")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Sumac")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Hawthorn")
m_genus_plot$diet <- relevel(m_genus_plot$diet, "Dandelion")

m_genus_plot$diet_initials <- as.factor(m_genus_plot$diet_initials)
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "D + S + H")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "D + S")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "S + H")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "D + H")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "S")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "H")
m_genus_plot$diet_initials <- relevel(m_genus_plot$diet_initials, "D")

# blast and rename taxa that were unclassified 
m_genus_plot2 <- within(m_genus_plot, Genus[OTU == "35b0bc6e20cc31eb6108a12530d1800b"] <- "Enterobacter*")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "a220dbde7245dfb7899ec7fe484525ed"] <- "Enterobacter*")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "ff2293ea3b41482e6d11f94b080b7ae7"] <- "Enterobacter*")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "29977f6c24ddd107d3db22f59a133da9"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "be9cf2697f7bbf30993c3ba0129c9109"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "f271e6e6d977397e39f1d55be04bf3ec"] <- "Enterobacter*")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "099f235697942b6e9afd171a3200e7b7"] <- "Enterobacter*")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "88a0786cb772eb069cbb66edd8419b0d"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "4a98edcaa9b90656232fd99aa8bd0c51"] <- "Enterobacter*")

m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "f72afeb78337970564ef6737a753a359"] <- "Pantoea")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "d9746fb0b36e4afb76dfec7e3e8925d3"] <- "Pantoea")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "e09bc90831bb219468010611697a8a90"] <- "Pantoea")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "c9fd02cbcdb891506aa96e92e6bdbbb3"] <- "Pantoea")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "ca310449cd52d0005b5281d9657d3963"] <- "Pantoea")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "c89b42ea6dc37965a4a48a7361a313cb"] <- "Pantoea")

m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "4353407195b1c7f4ee257f8b88fd532d"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "384112dd093eba5a3b56a44acd76ec02"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "b67c138c5ec29bb1a13dcefbb476b2d9"] <- "Rosenbergiella")
m_genus_plot2 <- within(m_genus_plot2, Genus[OTU == "16ac73b8c3dacc6564afb87b5497bd35"] <- "Rosenbergiella")

# Plot 
m.taxa.plot<-
  ggplot(m_genus_plot2, aes(x = Sample, y = Abundance, fill = Genus)) + 
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
   scale_fill_manual(values = c("firebrick", "hotpink3", "dodgerblue", "royalblue4", "darkseagreen", 
                                "gold", "orchid", "sienna2", "lightpink", "forestgreen", "grey")) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (>1%)\n") +
  ggtitle(expression(paste(italic("M. rotundata ")))) + 
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

m.taxa.plot

ggsave(
  "fig3b.jpg",
  plot = m.taxa.plot,
  width = 3000,
  height =1100,
  units = "px",
  dpi = 300
)
