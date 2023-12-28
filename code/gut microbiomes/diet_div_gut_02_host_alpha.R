##### The effect of diet diversity on bee gut microbes
##### Script 2 
##### Comparing microbe alpha diversity between two host species 

# load libraries and previous script 
library(tidyverse)
library(vegan)
library(emmeans)
library(bbmle)
library(DHARMa)
library(car)
source("scripts/diet_div_01.R")

#### SHANNON #### 

lab_metadata$No_pollen_spp <- as.factor(lab_metadata$No_pollen_spp)

# calculate shannon index 
lab_shannon <- as.data.frame(diversity(lab_feat_table_t, index = "shannon"))
# make row names a new column
lab_shannon <- rownames_to_column(lab_shannon, "index")
# merge shannon index with metadatafile
lab_shannon <- merge(lab_shannon, lab_metadata, by = "index")

# explore shannon diversity 
# by host
boxplot(lab_shannon$`diversity(lab_feat_table_t, index = "shannon")` ~ lab_shannon$host)
# by diet 
boxplot(lab_shannon$`diversity(lab_feat_table_t, index = "shannon")` ~ lab_shannon$diet)
boxplot(lab_shannon$`diversity(lab_feat_table_t, index = "shannon")` ~ lab_shannon$No_pollen_spp)
# distribution
hist(lab_shannon$`diversity(lab_feat_table_t, index = "shannon")`)

# model shannon diversity 
shannon_mod <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host * diet, data = lab_shannon)
shannon_mod2 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host + diet, data = lab_shannon)
shannon_mod3 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host, data = lab_shannon)
shannon_mod4 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ diet, data = lab_shannon)

shannon_mod5 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host * No_pollen_spp, data = lab_shannon)
shannon_mod6 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host + No_pollen_spp, data = lab_shannon)
shannon_mod7 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp, data = lab_shannon)

AICtab(shannon_mod, shannon_mod2, shannon_mod3, shannon_mod4,
       shannon_mod5, shannon_mod6, shannon_mod7) # 3 best although 6 and 5 are not far behind 
# dAIC df
# shannon_mod3  0.0 3 
# shannon_mod6  1.6 5 
# shannon_mod5  4.9 7 
# shannon_mod2  5.6 9 
# shannon_mod7  9.5 4 
# shannon_mod4 13.7 8 
# shannon_mod  15.3 15
simulateResiduals(shannon_mod3, plot = T) # levene test for homogeneity of variance significant 
simulateResiduals(shannon_mod6, plot = T) # all good 
simulateResiduals(shannon_mod5, plot = T) # also all good 

Anova(shannon_mod3)
# Response: diversity(lab_feat_table_t, index = "shannon")
# Sum Sq  Df F value   Pr(>F)   
# host       1.6007   1  9.8505 0.002042 **
#   Residuals 24.5379 151 

Anova(shannon_mod6)
# Response: diversity(lab_feat_table_t, index = "shannon")
# Sum Sq  Df F value   Pr(>F)   
# host           1.5806   1  9.6995 0.002208 **
#   No_pollen_spp  0.0949   1  0.5822 0.446644   
# Residuals     24.4430 150  

Anova(shannon_mod5, type = "III")
# Response: diversity(lab_feat_table_t, index = "shannon")
# Sum Sq  Df  F value  Pr(>F)    
# (Intercept)        67.781   1 414.8639 < 2e-16 ***
#   host                0.600   1   3.6732 0.05721 .  
# No_pollen_spp       0.000   1   0.0013 0.97182    
# host:No_pollen_spp  0.099   1   0.6060 0.43752    
# Residuals          24.344 149     

# plot shannon diversity 
host_shannon_mod <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ host, data = lab_shannon)
host_shannon_mod_means <- emmeans(host_shannon_mod, ~host) 
host.shannon.means.to.plot <- as.data.frame(summary(host_shannon_mod_means))
host.shannon.means.to.plot$tfupper <- host.shannon.means.to.plot$emmean + host.shannon.means.to.plot$SE
host.shannon.means.to.plot$tlower <- host.shannon.means.to.plot$emmean - host.shannon.means.to.plot$SE

host_shannon_mod_plot<-
  ggplot(host.shannon.means.to.plot, aes(x = host, color = host, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_point(stat ="identity", size = 6.5, aes(color = host)) +
  geom_errorbar(width=0.09,  aes(color = host), lwd = 1.2) +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 25),
        axis.text = element_text(size = 25, color = "black"),
        axis.title.y = element_text(size = 25)) +
  scale_color_manual(values = c("orchid4", "deepskyblue")) +
  theme(legend.position = "none")  +
  ylab("Shannon Index \n") +
  xlab("")

host_shannon_mod_plot

#### SPECIES RICHNESS ####
lab_rich <- as.data.frame(specnumber(lab_feat_table_t))
lab_rich <- rownames_to_column(lab_rich, "index")
lab_rich <- merge(lab_rich, lab_metadata, by = "index")

### explore data
boxplot(lab_rich$`specnumber(lab_feat_table_t)` ~ lab_rich$host)
boxplot(lab_rich$`specnumber(lab_feat_table_t)` ~ lab_rich$diet)
boxplot(lab_rich$`specnumber(lab_feat_table_t)` ~ lab_rich$No_pollen_spp)
hist(lab_rich$`specnumber(lab_feat_table_t)`)

### model richness 
rich.mod <- lm(`specnumber(lab_feat_table_t)` ~ host * diet, data = lab_rich)
rich.mod2 <- lm(`specnumber(lab_feat_table_t)` ~ host + diet, data = lab_rich)
rich.mod3 <- lm(`specnumber(lab_feat_table_t)` ~ host, data = lab_rich)
rich.mod4 <- lm(`specnumber(lab_feat_table_t)` ~ diet, data = lab_rich)

rich.mod5 <- lm(`specnumber(lab_feat_table_t)` ~ host * No_pollen_spp, data = lab_rich)
rich.mod6 <- lm(`specnumber(lab_feat_table_t)` ~ host + No_pollen_spp, data = lab_rich)
rich.mod7 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp, data = lab_rich)

AICtab(rich.mod, rich.mod2, rich.mod3, rich.mod4,
       rich.mod5, rich.mod6, rich.mod7) # mod 3 is best, although six is really close  
# dAIC df
# rich.mod3  0.0 3 
# rich.mod6  0.3 5 
# rich.mod5  2.2 7 
# rich.mod2  6.1 9 
# rich.mod  13.2 15
# rich.mod7 28.2 4 
# rich.mod4 34.3 8  
simulateResiduals(rich.mod3, plot = T) # looks good 

Anova(rich.mod3)
# Response: specnumber(lab_feat_table_t)
# Sum Sq  Df F value    Pr(>F)    
# host      1449.4   1  31.887 7.877e-08 ***
#   Residuals 6863.6 151       

Anova(rich.mod6)
# Response: specnumber(lab_feat_table_t)
# Sum Sq  Df F value    Pr(>F)    
# host          1448.2   1 32.2119 6.994e-08 ***
#   No_pollen_spp  164.6   2  1.8305    0.1639    
# Residuals     6699.0 149                  

### plot richness 
host_rich_mod <- lm(`specnumber(lab_feat_table_t)` ~ host, data = lab_rich)
host_rich_mod_means <- emmeans(host_rich_mod, ~host) 
host.means.to.plot <- as.data.frame(summary(host_rich_mod_means))
host.means.to.plot$tfupper <- host.means.to.plot$emmean + host.means.to.plot$SE
host.means.to.plot$tlower <- host.means.to.plot$emmean - host.means.to.plot$SE

host_rich_mod_plot<-
  ggplot(host.means.to.plot, aes(x = host, color = host,
                                 y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_point(stat ="identity", size = 6.5, aes(color = host)) +
  geom_errorbar(width=0.09, aes(color = host), lwd = 1.2) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 25, color = "black"),
      #  axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_color_manual(values = c("orchid4", "deepskyblue")) +
  theme(legend.position = "none")  +
  ylab(" ") + 
  xlab("")

host_rich_mod_plot

### Combine the shannon and richness plots 
host_div_plots <- ggarrange(host_rich_mod_plot, host_shannon_mod_plot,
                            ncol = 2, nrow = 1, 
                            labels = c("A", "B"))

host_div_plots_full <- annotate_figure(host_div_plots, bottom = text_grob("\n Bee Species", size = 18))

host_div_plots_full
