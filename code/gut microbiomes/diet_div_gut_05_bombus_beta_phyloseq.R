##### The effect of diet diversity on bee gut microbes
##### Script 5 
##### Microbiomes in B. impatiens using PHYLOSEQ

# libraries 
library(tidyverse)
library(phyloseq)
library(ape)
library(qiime2R)
library(vegan)
library(plyr)
library(ggpubr)
library(MASS)
library(DHARMa)
library(glmmTMB)
library(car)

# import data
b.physeq <- qza_to_phyloseq(
  features="data/gut-microbe-diversity-data/lab_bimp_filttable.qza",
  taxonomy="data/gut-microbe-diversity-data/lab_16STaxonomy.qza",
  metadata = "data/gut-microbe-diversity-data/lab_bimp_metadata3.tsv"
)

b.physeq

# create tree 
bimp_tree = rtree(ntaxa(b.physeq), rooted=TRUE, tip.label=taxa_names(b.physeq))
plot(bimp_tree)

# combine 
b.physeq1 <- merge_phyloseq(b.physeq, bimp_tree)
b.physeq1

#### Change NAs to be named as unclassified to the lowest possible rank 

for(i in 1:nrow(b.physeq1@tax_table)) { # go through each row  
  #  print(colnames(temp)[i]) # for debugging 
  for (j in 1:ncol(b.physeq1@tax_table)) { # then go through each column 
    #  print(temp[i,j]) # for debugging
    if (is.na(b.physeq1@tax_table[i,j])) { # if it is NA 
      if(startsWith(b.physeq1@tax_table[i,j-1], # AND if cell in previous column starts with "unclassified" 
                    "Unclassified")){
        v <- strsplit(b.physeq1@tax_table[i,j-1], split = " ") # then split the previous cell contents into two strings 
        b.physeq1@tax_table[i,j] <- paste("Unclassified", unlist(v[[1]][2])) # assign as "unclassified" + last taxanomic level 
        #  print(paste("print this", unlist(v[[1]][2]))) # for debugging
      } else {b.physeq1@tax_table[i,j] <- paste("Unclassified", b.physeq1@tax_table[i,j-1]) # if NA and not unclassified in previous column, assign as unclassified + taxonomic level in previous cell 
      }
    }
  }
}

# make things factors 
b.physeq1@sam_data$No_pollen_spp <- as.factor(b.physeq1@sam_data$No_pollen_spp)
b.physeq1@sam_data$diet <- as.factor(b.physeq1@sam_data$diet)
b.physeq1@sam_data$Age <- as.factor(b.physeq1@sam_data$Age)
b.physeq1@sam_data$Colony <- as.factor(b.physeq1@sam_data$Colony)
b.physeq1@sam_data$Emerge_start_date <- as.factor(b.physeq1@sam_data$Emerge_start_date)

# ORDINATIONS BY DIET 

# ordinations using PCoA 
bimp_ord <- ordinate(b.physeq1, method = "PCoA", distance = "bray")
bimp_ord1 <- ordinate(b.physeq1, method = "PCoA", distance = "wunifrac")
bimp_ord2 <- ordinate(b.physeq1, method = "PCoA", distance = "unifrac")

### plot by diet 
bc_diet <- plot_ordination(b.physeq1, bimp_ord, color = "diet")
bc_diet <- bc_diet + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  ggtitle("Bray-Curtis") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

wuni_diet <- plot_ordination(b.physeq1, bimp_ord1, color = "diet")
wuni_diet <- wuni_diet + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  ggtitle("Weighted UniFrac") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

uni_diet <- plot_ordination(b.physeq1, bimp_ord2, color = "diet")
uni_diet <- uni_diet + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  ggtitle("Unweighted UniFrac") + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) 

#### PLOT BY DIET DIVERSITY 

bc_diet_div <- plot_ordination(b.physeq1, bimp_ord, color = "No_pollen_spp")
bc_diet_div <- bc_diet_div + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
#  ggtitle("Bray-Curtis") + 
  theme_classic() +
  theme(legend.position = "none")

wuni_diet_div <- plot_ordination(b.physeq1, bimp_ord1, color = "No_pollen_spp")
wuni_diet_div <- wuni_diet_div + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
#  ggtitle("Weighted UniFrac") + 
  theme_classic() +
  theme(legend.position = "none")

uni_diet_div <- plot_ordination(b.physeq1, bimp_ord2, color = "No_pollen_spp")
uni_diet_div <- uni_diet_div + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
#  ggtitle("Unweighted UniFrac") + 
  theme_classic() +
  theme(legend.position = "none")

# PLOT BY AGE ####

bc_age <- plot_ordination(b.physeq1, bimp_ord, color = "Age")
bc_age<-bc_age + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
 # ggtitle("Bray-Curtis") + 
  theme_classic() +
  scale_color_manual(values = c("mediumorchid",
                                "darkorchid4")) + 
  theme(legend.position = "none")

wuni_age <- plot_ordination(b.physeq1, bimp_ord1, color = "Age")
wuni_age<-wuni_age + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
#  ggtitle("Weighted UniFrac") + 
  theme_classic() +
  scale_color_manual(values = c("mediumorchid",
                                "darkorchid4")) + 
  theme(legend.position = "none")

uni_age <- plot_ordination(b.physeq1, bimp_ord2, color = "Age")
uni_age<-uni_age + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
 # ggtitle("Unweighted UniFrac") + 
  theme_classic() +
  scale_color_manual(values = c("mediumorchid",
                                "darkorchid4")) + 
  theme(legend.position = "none")

# PLOT BY COLONY #####

bc_colony <- plot_ordination(b.physeq1, bimp_ord, color = "Colony")
bc_colony<-bc_colony + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
 # ggtitle("Bray-Curtis") + 
  theme_classic() + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick","chocolate1", "forestgreen",
                                "yellowgreen","royalblue4", "darkorchid")) 

wuni_colony <- plot_ordination(b.physeq1, bimp_ord1, color = "Colony")
wuni_colony<-wuni_colony + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
 # ggtitle("Weighted UniFrac") + 
  theme_classic() + 
  theme(legend.position = "none")  +
  scale_color_manual(values = c("firebrick","chocolate1", "forestgreen",
                                "yellowgreen","royalblue4", "darkorchid"))

uni_colony <- plot_ordination(b.physeq1, bimp_ord2, color = "Colony")
uni_colony<-uni_colony + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
 # ggtitle("Unweighted UniFrac") + 
  theme_classic() + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick","chocolate1", "forestgreen",
                                "yellowgreen","royalblue4", "darkorchid"))

bc_date <- plot_ordination(b.physeq1, bimp_ord, color = "Emerge_start_date")
bc_date<-bc_date + 
  stat_ellipse(geom = "polygon",
               aes(fill = host),
               alpha = 0,
               level = 0.95) + 
  # ggtitle("Bray-Curtis") + 
  theme_classic() +
  # scale_color_manual(values = c("mediumorchid",
  #                               "darkorchid4")) + 
  theme(legend.position = "right")


#### COMBINE 12 PLOTS FOR FIGURE 
bombus_pcoa_plots <- ggarrange(bc_diet, wuni_diet, uni_diet,
                               bc_diet_div, wuni_diet_div, uni_diet_div,
                               bc_colony, wuni_colony, uni_colony,
                               bc_age, wuni_age, uni_age,
                               nrow = 4, ncol = 3)
bombus_pcoa_plots

# Compare outputs with different distance metrics 
b.df <- data.frame(sample_data(b.physeq1))

# WEIGHTED UNIFRAC
bimp.wunifrac <- distance(b.physeq1, method = "wunifrac")

# Betadisper with diet comp 
b.wunifrac.betadisper <- betadisper(bimp.wunifrac, b.df$diet)
anova(b.wunifrac.betadisper)
summary(b.wunifrac.betadisper)

# Betadisper with diet div 
b_div.wunifrac.betadisper <- betadisper(bimp.wunifrac, b.df$No_pollen_spp)
anova(b_div.wunifrac.betadisper)
summary(b_div.wunifrac.betadisper)

# betadisper with age
b_age.wunifrac.betadisper <- betadisper(bimp.wunifrac, b.df$Age)
anova(b_age.wunifrac.betadisper)

# betadisper with colony
b_col.wunifrac.betadisper <- betadisper(bimp.wunifrac, b.df$Colony)
anova(b_col.wunifrac.betadisper)

# adonis 
adonis2(bimp.wunifrac ~ diet + Age + Colony, data = b.df) 
adonis2(bimp.wunifrac ~ No_pollen_spp + Age + Colony, data = b.df) 

# UNWEIGHTED UNIFRAC 
bimp.unifrac <- distance(b.physeq1, method="unifrac")

# betadisper with diet comp 
b.unifrac.betadisper <- betadisper(bimp.unifrac, b.df$diet)
anova(b.unifrac.betadisper)

# betadisper with diet div
b_div.unifrac.betadisper <- betadisper(bimp.unifrac, b.df$No_pollen_spp)
anova(b_div.unifrac.betadisper)

# betadisper with age 
b_age.unifrac.betadisper <- betadisper(bimp.unifrac, b.df$Age)
anova(b_age.unifrac.betadisper)

# betadisper with colony
b_col.unifrac.betadisper <- betadisper(bimp.unifrac, b.df$Colony)
anova(b_col.unifrac.betadisper)

# betadisper with adonis 
adonis2(bimp.unifrac ~ diet + Age + Colony, data = b.df) 
adonis2(bimp.unifrac ~ No_pollen_spp + Age + Colony, data = b.df) 

# BRAY CURTIS
bimp.bray <- distance(b.physeq1, method="bray")

# betadisper with diet comp
b.bray.betadisper <- betadisper(bimp.bray, b.df$diet)
anova(b.bray.betadisper)

# betadisper with diet div 
b.bray.betadisper <- betadisper(bimp.bray, b.df$No_pollen_spp)
anova(b.bray.betadisper)

# betadisper with age 
b_age.bray.betadisper <- betadisper(bimp.bray, b.df$Age)
anova(b_age.bray.betadisper)

# betadisper with colony 
b_col.bray.betadisper <- betadisper(bimp.bray, b.df$Colony)
anova(b_col.bray.betadisper)

# adonis 
adonis2(bimp.bray ~ diet + Age + Colony, data = b.df)
# Df SumOfSqs      R2      F Pr(>F)  
# diet      6   1.6169 0.07366 1.0729  0.363  
# Age       1   0.4483 0.02042 1.7849  0.068 .
# Colony    5   1.7999 0.08200 1.4332  0.050 *
# Residual 72  18.0850 0.82391                
# Total    84  21.9502 1.00000                

adonis2(bimp.bray ~ No_pollen_spp + Age + Colony, data = b.df)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bimp.bray ~ No_pollen_spp + Age + Colony, data = b.df)
# Df SumOfSqs      R2      F Pr(>F)  
# No_pollen_spp  2   0.7274 0.03314 1.4576  0.101  
# Age            1   0.3911 0.01782 1.5674  0.123  
# Colony         5   1.8673 0.08507 1.4966  0.035 *
#   Residual      76  18.9644 0.86397                
# Total         84  21.9502 1.00000

# BAR PLOTS ####
b_genus_plot <- b.physeq1 %>%
  #tax_glom(taxrank = "Genus") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                         
  arrange(Genus) 

b_genus_plot$diet <- as.factor(b_genus_plot$diet)
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Dandelion + Sumac + Hawthorn")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Dandelion + Sumac")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Sumac + Hawthorn")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Dandelion + Hawthorn")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Sumac")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Hawthorn")
b_genus_plot$diet <- relevel(b_genus_plot$diet, "Dandelion")

b_genus_plot$diet_initials <- as.factor(b_genus_plot$diet_initials)
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "D + S + H")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "D + S")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "S + H")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "D + H")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "S")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "H")
b_genus_plot$diet_initials <- relevel(b_genus_plot$diet_initials, "D")

# Blast and rename the taxa that did not get classified 

b.tax_table <- as.data.frame(b.physeq1@tax_table)
b.tax_table <- rownames_to_column(b.tax_table, "OTU")
b.genus_table <- b.tax_table %>% select("OTU", "Genus")
b.species_table <- b.tax_table %>% select("OTU", "Genus", "Species")

b_genus_plot <- within(b_genus_plot, Genus[OTU == "4a8e08ef613b9fc1d1523112af5aa01e"] <- "Bombus")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "f72afeb78337970564ef6737a753a359"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "a220dbde7245dfb7899ec7fe484525ed"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "29977f6c24ddd107d3db22f59a133da9"] <- "Rosenbergiella")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "f72afeb78337970564ef6737a753a359"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "d9746fb0b36e4afb76dfec7e3e8925d3"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "e09bc90831bb219468010611697a8a90"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "c9fd02cbcdb891506aa96e92e6bdbbb3"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "ca310449cd52d0005b5281d9657d3963"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "c89b42ea6dc37965a4a48a7361a313cb"] <- "Pantoea")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "35b0bc6e20cc31eb6108a12530d1800b"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "a220dbde7245dfb7899ec7fe484525ed"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "ff2293ea3b41482e6d11f94b080b7ae7"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "f271e6e6d977397e39f1d55be04bf3ec"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "099f235697942b6e9afd171a3200e7b7"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "4a98edcaa9b90656232fd99aa8bd0c51"] <- "Enterobacter*")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "bc7a67162e397d0b74121e9fe9602d08"] <- "Klebsiella") 
b_genus_plot <- within(b_genus_plot, Genus[OTU == "df208d4225db4d84097d4635026c0ada"] <- "Klebsiella") 
b_genus_plot <- within(b_genus_plot, Genus[OTU == "646be12e02a43f09ab0c5065ade76df3"] <- "Klebsiella") 
b_genus_plot <- within(b_genus_plot, Genus[OTU == "5c80ac19f5637f31515c5ff3f1c59307"] <- "Klebsiella")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "e47d852d8b0df2936195622ef95be11b"] <- "Klebsiella")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "fe041460dc4756ea5ff0edf699e6c7e3"] <- "Klebsiella")
b_genus_plot <- within(b_genus_plot, Genus[OTU == "8f28a8b45e0c7e53c1ea4138de6d731b"] <- "Klebsiella")

# remove the underscore from Candidatus Schmidhempelia 
levels(b_genus_plot$Genus)[levels(b_genus_plot$Genus) == "Candidatus_Schmidhempelia" ] <- "Candidatus Schmidhempelia"

# Plot 
b.taxa.plot<-
  ggplot(b_genus_plot, aes(x = Sample, y = Abundance, fill = Genus)) + 
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
  scale_fill_manual(values = c("grey", "lightpink", "lightskyblue", "hotpink3",
                               "royalblue3", "gold", "forestgreen", "olivedrab3", "dodgerblue",
                               "tan", "red", "green", "sienna2", "brown", "orchid", "lightblue", "orange")) +
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, nrow = 2, byrow = T)) +
  ylab("Relative Abundance (>1%)\n") +
  ggtitle(expression(paste(italic("B. impatiens ")))) 

b.taxa.plot

### ggarrange with Megarot taxa plot (script 7)

taxa_plots <- ggarrange(b.taxa.plot, m.taxa.plot,
                            ncol = 1, nrow = 2, 
                            labels = c("A", "B"))
taxa_plots



