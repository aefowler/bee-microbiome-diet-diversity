##### The effect of diet diversity on bee gut microbes
##### Script 6 
##### Microbial alpha diversity in M. rotundata 

library(tidyverse)
library(vegan)
library(lme4)
library(emmeans)
library(ggpubr)

## Shannon diversity

# separate b. impatiens from full data set
megarot_shannon <- lab_shannon %>% filter(host == "M. rotundata")

# explore data set 
hist(megarot_shannon$`diversity(lab_feat_table_t, index = "shannon")`)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ diet, data = megarot_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp, data = megarot_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ Emerge_start_date, data = megarot_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ plate, data = megarot_shannon)
plot(`diversity(lab_feat_table_t, index = "shannon")` ~ Wing, data = megarot_shannon)

megarot_shannon$plate <- as.factor(megarot_shannon$plate)
megarot_shannon$Emerge_start_date <- as.factor(megarot_shannon$Emerge_start_date)

# model
megarot_shannon_mod1 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ diet + Wing + (1|plate) + (1|Emerge_start_date), 
                            data = megarot_shannon)

megarot_shannon_mod2 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ diet + Wing + (1|Emerge_start_date), 
                             data = megarot_shannon)

megarot_shannon_mod3 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ diet + Wing + (1|plate), 
                             data = megarot_shannon)

megarot_shannon_mod4 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ diet + Wing, data = megarot_shannon)

megarot_shannon_mod5 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ diet, data = megarot_shannon)

megarot_shannon_mod6 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp + Wing + (1|plate) + (1|Emerge_start_date), 
                            data = megarot_shannon)

megarot_shannon_mod7 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp + Wing + (1|Emerge_start_date), 
                             data = megarot_shannon)

megarot_shannon_mod8 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp + Wing + (1|plate), 
                             data = megarot_shannon)

megarot_shannon_mod9 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp + Wing, data = megarot_shannon)

megarot_shannon_mod10 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp, data = megarot_shannon)

AICtab(megarot_shannon_mod1, megarot_shannon_mod2, 
       megarot_shannon_mod3, megarot_shannon_mod4,
       megarot_shannon_mod5, megarot_shannon_mod6, 
       megarot_shannon_mod7, megarot_shannon_mod8, 
       megarot_shannon_mod9, megarot_shannon_mod10) 

# dAIC df
# megarot_shannon_mod10  0.0 4 <- use for diet div
# megarot_shannon_mod9   1.4 5 
# megarot_shannon_mod5   3.0 8 <- use for diet comp 
# megarot_shannon_mod4   4.8 9 
# megarot_shannon_mod7   6.8 6 
# megarot_shannon_mod8   6.8 6 
# megarot_shannon_mod6   8.8 7 
# megarot_shannon_mod2  16.0 10
# megarot_shannon_mod3  16.0 10
# megarot_shannon_mod1  18.0 11

simulateResiduals(megarot_shannon_mod10, plot = T) # looks good 
simulateResiduals(megarot_shannon_mod5, plot = T) # looks good 

Anova(megarot_shannon_mod10)
# Response: diversity(lab_feat_table_t, index = "shannon")
# Sum Sq Df F value Pr(>F)
# No_pollen_spp  0.3279  2  0.9529 0.3909
# Residuals     11.1845 65

Anova(megarot_shannon_mod5)
# Response: diversity(lab_feat_table_t, index = "shannon")
#            Sum Sq Df F value Pr(>F)
# diet       1.1274  6  1.1037 0.3707
# Residuals 10.3850 61 

summary(glht(megarot_shannon_mod5, linfct = mcp(diet = "Tukey")))

### plot

megarot_shannon_mod5 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ diet_initials, data = megarot_shannon)
megarot_shannon_mod_means <- emmeans(megarot_shannon_mod5, ~diet_initials) 
megarot.shannon.means.to.plot <- as.data.frame(summary(megarot_shannon_mod_means))
megarot.shannon.means.to.plot$tfupper <- megarot.shannon.means.to.plot$emmean + megarot.shannon.means.to.plot$SE
megarot.shannon.means.to.plot$tlower <- megarot.shannon.means.to.plot$emmean - megarot.shannon.means.to.plot$SE

megarot_shannon_mod_plot<-
  ggplot(megarot.shannon.means.to.plot, aes(x = diet_initials, color = diet_initials, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_point(stat ="identity", size = 6, aes(color = diet_initials)) +
  geom_errorbar(width=0.065,  aes(color = diet_initials)) +
  theme(text = element_text(size=18)) +
  scale_color_manual(values = c("deepskyblue",
                                "deepskyblue",
                                "deepskyblue",
                                "dodgerblue3",
                                "dodgerblue3",
                                "dodgerblue3",
                                "dodgerblue4")) +
  theme(legend.position = "none")  +
  ylab("Shannon Index \n") +
  xlab(" \n Diet Treatment") + 
  ylim(1.7,2.75) + 
  geom_signif(y_position = 2.6, 
              xmin = 0.75, 
              xmax = 3.25,
              annotation = "One Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") +
  geom_signif(y_position = 2.6, 
              xmin = 3.75, 
              xmax = 6.25,
              annotation = "Two Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") +
  geom_signif(y_position = 2.6, 
              xmin = 6.55, 
              xmax = 7.45,
              annotation = "Three Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") 

megarot_shannon_mod_plot

## Species richness 

# separate megarot from full data set
megarot_rich <- lab_rich %>% filter(host == "M. rotundata")

# explore data set 
hist(megarot_rich$`specnumber(lab_feat_table_t)`)
boxplot(`specnumber(lab_feat_table_t)` ~ diet, data = megarot_rich)
boxplot(`specnumber(lab_feat_table_t)` ~ plate, data = megarot_rich)
boxplot(`specnumber(lab_feat_table_t)` ~ Emerge_start_date, data = megarot_rich)
plot(`specnumber(lab_feat_table_t)` ~ Wing, data = megarot_rich)

# model
megarot_rich_mod1 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + (1|plate) + (1|Emerge_start_date), 
                         data = megarot_rich)

megarot_rich_mod2 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + (1|Emerge_start_date), 
                          data = megarot_rich)

megarot_rich_mod3 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + (1|plate), 
                          data = megarot_rich)

megarot_rich_mod4 <- lm(`specnumber(lab_feat_table_t)` ~ diet + Wing, data = megarot_rich)

megarot_rich_mod5 <- lm(`specnumber(lab_feat_table_t)` ~ diet, data = megarot_rich)

megarot_rich_mod3.5 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + (1|plate), 
                            data = megarot_rich)

megarot_rich_mod6 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + (1|plate) + (1|Emerge_start_date), 
                          data = megarot_rich)

megarot_rich_mod7 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + (1|Emerge_start_date), 
                          data = megarot_rich)

megarot_rich_mod8 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + (1|plate), 
                          data = megarot_rich)

megarot_rich_mod9 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing, data = megarot_rich)

megarot_rich_mod10 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp, data = megarot_rich)

megarot_rich_mod8.5 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + (1|plate), 
                          data = megarot_rich)

AICtab(megarot_rich_mod1, megarot_rich_mod2, 
       megarot_rich_mod3, megarot_rich_mod4, 
       megarot_rich_mod5, megarot_rich_mod6, 
       megarot_rich_mod7, megarot_rich_mod8,
       megarot_rich_mod9, megarot_rich_mod10,
       megarot_rich_mod3.5, megarot_rich_mod8.5) 

# dAIC df
# megarot_rich_mod3    0.0 10 <- use for diet comp 
# megarot_rich_mod2    1.2 10
# megarot_rich_mod1    2.0 11
# megarot_rich_mod8   12.7 6  <- use for diet div
# megarot_rich_mod7   13.5 6 
# megarot_rich_mod6   14.4 7 
# megarot_rich_mod3.5 17.0 9 
# megarot_rich_mod9   29.6 5 
# megarot_rich_mod8.5 29.7 5 
# megarot_rich_mod4   33.1 9 
# megarot_rich_mod10  35.2 4 
# megarot_rich_mod5   38.5 8 

simulateResiduals(megarot_rich_mod3, plot = T) # looks good
Anova(megarot_rich_mod3)
# Response: specnumber(lab_feat_table_t)
#       Chisq Df Pr(>Chisq)
# diet 7.2942  6     0.2945
# Wing 0.0032  1     0.9545
summary(glht(megarot_rich_mod3, linfct = mcp(diet = "Tukey")))

simulateResiduals(megarot_rich_mod8, plot = T) # looks good
Anova(megarot_rich_mod8)
# Response: specnumber(lab_feat_table_t)
# Chisq Df Pr(>Chisq)
# No_pollen_spp 2.4696  2     0.2909
# Wing          0.1142  1     0.7354

# plot model
megarot_rich_mod3 <- lmer(`specnumber(lab_feat_table_t)` ~ diet_initials + Wing + (1|plate), 
                          data = megarot_rich)
megarot_rich_mod_means <- emmeans(megarot_rich_mod3, ~diet_initials) 
megarot.rich.means.to.plot <- as.data.frame(summary(megarot_rich_mod_means))
megarot.rich.means.to.plot$tfupper <- megarot.rich.means.to.plot$emmean + megarot.rich.means.to.plot$SE
megarot.rich.means.to.plot$tlower <- megarot.rich.means.to.plot$emmean - megarot.rich.means.to.plot$SE

megarot.rich.means.to.plot$diet_initials2 <- factor(c("D", "H", "S", "DH", "SH", "DS", "DSH"))
megarot.rich.means.to.plot$diet_initials2 <- factor(megarot.rich.means.to.plot$diet_initials2, 
                                                    levels = c("D", "H", "S", "DH", "SH", "DS", "DSH"))

# bar plot
megarot_rich_mod_barplot<-
  ggplot(megarot.rich.means.to.plot, aes(x = diet_initials2, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_bar(stat ="identity", 
           color = "black", 
           fill = "deepskyblue",
           width = 0.6) +
  geom_errorbar(width=0.065, lwd = 0.75) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size=25, 
                                 color = "black")) +
     #   axis.title.x = element_blank(),
      #  axis.title.y = element_blank(),
       # plot.title = element_text(size = 22, face = "bold")
  theme(legend.position = "none")  +
#  ggtitle(expression(paste(italic("M. rotundata")))) +
  ylab(" ") +
  xlab(" \n Pollen Diet Type") +
  ylim(0,30)  +
  geom_text(x = 1.5, y = 28, label = expression(paste("Diet: ", italic("P"), "= 0.294")), size = 6)

megarot_rich_mod_barplot

## Plot for diet diversity 
megarot_rich_div_mod_means <- emmeans(megarot_rich_mod8, ~No_pollen_spp) 
megarot.rich.means.to.plot2 <- as.data.frame(summary(megarot_rich_div_mod_means))
megarot.rich.means.to.plot2$tfupper <- megarot.rich.means.to.plot2$emmean + megarot.rich.means.to.plot2$SE
megarot.rich.means.to.plot2$tlower <- megarot.rich.means.to.plot2$emmean - megarot.rich.means.to.plot2$SE

megarot_rich_div_mod_barplot<-
  ggplot(megarot.rich.means.to.plot2, aes(x = No_pollen_spp, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_bar(stat ="identity", 
           color = "black", 
           fill = "deepskyblue",
           width = 0.6) +
  geom_errorbar(width=0.065, lwd = 0.75) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size=25,
                                 color = "black")) +
     #   axis.title.x = element_blank(),
       # axis.title.y = element_blank()) +
  theme(legend.position = "none")  +
 # ggtitle(expression(paste(italic("M. rotundata")))) +
  ylab(" ") +
  xlab(" \n Pollen Diet Diversity (# species)") +
  ylim(0,30) + 
  geom_text(x = 1, y = 28, label = expression(paste("Diet: ", italic("P"), "= 0.735")), size = 6)

megarot_rich_div_mod_barplot

## arrange with bombus 
rich_div_plots <- ggarrange(bimp_rich_div_mod_barplot,
                            megarot_rich_div_mod_barplot,
                            nrow = 2, ncol = 1)

rich_div_plots <- annotate_figure(rich_div_plots, left = text_grob("Species Richness \n", rot = 90, size = 18))
rich_div_plots

### ARRANGE ALL FOUR PLOTS 
rich_plots <- ggarrange(bimp_rich_mod_barplot, bimp_rich_div_mod_barplot,
                        megarot_rich_mod_barplot, megarot_rich_div_mod_barplot,
                        nrow = 2, ncol = 2)

rich_plots <- annotate_figure(rich_plots, left = text_grob("Species Richness \n", rot = 90, size = 18))
rich_plots
