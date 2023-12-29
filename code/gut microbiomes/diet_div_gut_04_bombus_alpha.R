##### The effect of diet diversity on bee gut microbes
##### Script 4 
##### Microbial alpha diversity in B. impatiens 

library(tidyverse)
library(vegan)
library(lme4)

source("diet_div_gut_01.R")
source("diet_div_gut_02_host_alpha.R")

#### SHANNON #### 

# filter out bombus from the second script where we calculated shannon for all samples 
bimp_shannon <- lab_shannon %>% filter(host == "B. impatiens")

# explore shannon 
hist(bimp_shannon$`diversity(lab_feat_table_t, index = "shannon")`)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ diet, data = bimp_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ No_pollen_spp, data = bimp_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ Age, data = bimp_shannon)
boxplot(`diversity(lab_feat_table_t, index = "shannon")` ~ Emerge_start_date, data = bimp_shannon)
plot(`diversity(lab_feat_table_t, index = "shannon")` ~ Wing, data = bimp_shannon)

# model shannon with diet composition 
bimp_shannon_mod <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           diet + Wing + Age + (1|Colony) + (1|Emerge_start_date), 
                         data = bimp_shannon)

bimp_shannon_mod1 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                            diet * Age + Wing + (1|Colony) + (1|Emerge_start_date), 
                          data = bimp_shannon)

bimp_shannon_mod2 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                            diet + Wing + Age + (1|Emerge_start_date), 
                          data = bimp_shannon)

bimp_shannon_mod3 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                            diet + Wing + Age + (1|Colony), 
                          data = bimp_shannon)

bimp_shannon_mod4 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                          diet + Wing + Age, data = bimp_shannon)

bimp_shannon_mod5 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                          diet + Wing, data = bimp_shannon)

bimp_shannon_mod6 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                          diet + Age, data = bimp_shannon)

bimp_shannon_mod6int <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             diet * Age, data = bimp_shannon)

bimp_shannon_mod7 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                          diet, data = bimp_shannon)

# model shannon with diet diversity 
bimp_shannon_mod8 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           No_pollen_spp + Wing + Age + (1|Colony) + (1|Emerge_start_date), 
                         data = bimp_shannon)

bimp_shannon_mod9 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                            No_pollen_spp * Age + Wing + (1|Colony) + (1|Emerge_start_date), 
                          data = bimp_shannon)

bimp_shannon_mod10 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             No_pollen_spp + Wing + Age + (1|Emerge_start_date), 
                          data = bimp_shannon)

bimp_shannon_mod11 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             No_pollen_spp + Wing + Age + (1|Colony), 
                          data = bimp_shannon)

bimp_shannon_mod12 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           No_pollen_spp + Wing + Age, data = bimp_shannon)

bimp_shannon_mod13 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           No_pollen_spp + Wing, data = bimp_shannon)

bimp_shannon_mod14 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           No_pollen_spp + Age, data = bimp_shannon)

bimp_shannon_mod15 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             diet * Age, data = bimp_shannon)

bimp_shannon_mod16 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                           No_pollen_spp, data = bimp_shannon)

bimp_shannon_mod17 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             No_pollen_spp + (1|Colony) + (1|Emerge_start_date), data = bimp_shannon)

bimp_shannon_mod18 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             No_pollen_spp + (1|Emerge_start_date), data = bimp_shannon)

bimp_shannon_mod19 <- lmer(`diversity(lab_feat_table_t, index = "shannon")` ~ 
                             No_pollen_spp + (1|Colony), data = bimp_shannon)

AICtab(bimp_shannon_mod, 
       bimp_shannon_mod1,
       bimp_shannon_mod2,
       bimp_shannon_mod3,
       bimp_shannon_mod4,
       bimp_shannon_mod5, 
       bimp_shannon_mod6,
       bimp_shannon_mod6int,
       bimp_shannon_mod7,
       bimp_shannon_mod8,
       bimp_shannon_mod9,
       bimp_shannon_mod10,
       bimp_shannon_mod11,
       bimp_shannon_mod12,
       bimp_shannon_mod13,
       bimp_shannon_mod14,
       bimp_shannon_mod15,
       bimp_shannon_mod16,
       bimp_shannon_mod17, 
       bimp_shannon_mod18,
       bimp_shannon_mod19)

# dAIC df
# bimp_shannon_mod16    0.0 4 -> use for diet div
# bimp_shannon_mod13    1.6 5 
# bimp_shannon_mod14    1.9 5 
# bimp_shannon_mod12    3.5 6 
# bimp_shannon_mod7     7.5 8 -> use for diet comp
# bimp_shannon_mod5     9.1 9 
# bimp_shannon_mod6     9.5 9 
# bimp_shannon_mod4    11.0 10
# bimp_shannon_mod10   14.1 7 
# bimp_shannon_mod6int 14.6 15
# bimp_shannon_mod15   14.6 15
# bimp_shannon_mod8    14.9 8 
# bimp_shannon_mod11   15.2 7 
# bimp_shannon_mod9    19.0 10
# bimp_shannon_mod2    28.2 11
# bimp_shannon_mod     29.1 12
# bimp_shannon_mod3    29.4 11
# bimp_shannon_mod1    36.4 18

simulateResiduals(bimp_shannon_mod16, plot = T)# looks good
simulateResiduals(bimp_shannon_mod7, plot = T)# looks good

Anova(bimp_shannon_mod16)
# Response: diversity(lab_feat_table_t, index = "shannon")
# Sum Sq Df F value Pr(>F)
# No_pollen_spp  0.1661  2  0.5295 0.5909
# Residuals     12.8593 82

Anova(bimp_shannon_mod7)
# Response: diversity(lab_feat_table_t, index = "shannon")
#            Sum Sq Df F value Pr(>F)
# diet       0.2354  6  0.2393 0.9622
# Residuals 12.7900 78  

summary(glht(bimp_shannon_mod7, linfct = mcp(diet = "Tukey")))

# plot shannon 
bimp_shannon_mod7 <- lm(`diversity(lab_feat_table_t, index = "shannon")` ~ diet_initials, data = bimp_shannon)
bimp_shannon_mod_means <- emmeans(bimp_shannon_mod7, ~diet_initials) 
bimp.shannon.means.to.plot <- as.data.frame(summary(bimp_shannon_mod_means))
bimp.shannon.means.to.plot$tfupper <- bimp.shannon.means.to.plot$emmean + bimp.shannon.means.to.plot$SE
bimp.shannon.means.to.plot$tlower <- bimp.shannon.means.to.plot$emmean - bimp.shannon.means.to.plot$SE

bimp_shannon_mod_plot<-
  ggplot(bimp.shannon.means.to.plot, aes(x = diet_initials, color = diet_initials, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_point(stat ="identity", size = 6, aes(color = diet_initials)) +
  geom_errorbar(width=0.065,  aes(color = diet_initials)) +
  theme(text = element_text(size=18)) +
  scale_color_manual(values = c("mediumorchid",
                                "mediumorchid",
                                "mediumorchid",
                                "darkorchid3",
                                "darkorchid3",
                                "darkorchid3",
                                "darkorchid4")) +
  theme(legend.position = "none")  +
  # ggtitle(expression(paste(italic("B. impatiens")))) +
  ylab("Shannon Index \n") +
  xlab(" \n Diet Treatment") + 
  ylim(1.85,3) + 
  geom_signif(y_position = 2.85, 
              xmin = 0.75, 
              xmax = 3.25,
              annotation = "One Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") +
  geom_signif(y_position = 2.85, 
              xmin = 3.75, 
              xmax = 6.25,
              annotation = "Two Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") +
  geom_signif(y_position = 2.85, 
              xmin = 6.55, 
              xmax = 7.45,
              annotation = "Three Species", 
              tip_length = 0.05,
              textsize = 4,
              color = "black") 

bimp_shannon_mod_plot

#### SPECIES RICHNESS ####

# separate b. impatiens from full richness 
bimp_rich <- lab_rich %>% filter(host == "B. impatiens")

# explore data set 
hist(bimp_rich$`specnumber(lab_feat_table_t)`)
boxplot(`specnumber(lab_feat_table_t)` ~ diet, data = bimp_rich)
boxplot(`specnumber(lab_feat_table_t)` ~ No_pollen_spp, data = bimp_rich)
boxplot(`specnumber(lab_feat_table_t)` ~ Age, data = bimp_rich)
boxplot(`specnumber(lab_feat_table_t)` ~ Emerge_start_date, data = bimp_rich)
plot(`specnumber(lab_feat_table_t)` ~ Wing, data = bimp_rich)

# model

bimp_rich_mod <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + Age + (1|Colony) + (1|Emerge_start_date), 
                      data = bimp_rich)

bimp_rich_mod1 <- lmer(`specnumber(lab_feat_table_t)` ~ diet * Age + Wing + (1|Colony) + (1|Emerge_start_date), 
                       data = bimp_rich)

bimp_rich_mod2 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + Age + (1|Emerge_start_date), 
                       data = bimp_rich)

bimp_rich_mod3 <- lmer(`specnumber(lab_feat_table_t)` ~ diet + Wing + Age + (1|Colony), 
                       data = bimp_rich)

bimp_rich_mod4 <- lm(`specnumber(lab_feat_table_t)` ~ diet + Wing + Age, data = bimp_rich)

bimp_rich_mod5 <- lm(`specnumber(lab_feat_table_t)` ~ diet + Wing, data = bimp_rich)

bimp_rich_mod6 <- lm(`specnumber(lab_feat_table_t)` ~ diet + Age, data = bimp_rich)

bimp_rich_mod6int <- lm(`specnumber(lab_feat_table_t)` ~ diet * Age, data = bimp_rich)

bimp_rich_mod7 <- lm(`specnumber(lab_feat_table_t)` ~ diet, data = bimp_rich)

bimp_rich_mod8 <- lmer(`specnumber(lab_feat_table_t)` ~ diet * Age + (1|Colony) + (1|Emerge_start_date), 
                       data = bimp_rich)

# again with diet diversity 
bimp_rich_mod9 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + Age + (1|Colony) + (1|Emerge_start_date), 
                      data = bimp_rich)

bimp_rich_mod10 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp * Age + Wing + (1|Colony) + (1|Emerge_start_date), 
                       data = bimp_rich)

bimp_rich_mod11 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + Age + (1|Emerge_start_date), 
                       data = bimp_rich)

bimp_rich_mod12 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + Age + (1|Colony), 
                       data = bimp_rich)

bimp_rich_mod13 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing + Age, data = bimp_rich)

bimp_rich_mod14 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Wing, data = bimp_rich)

bimp_rich_mod15 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp + Age, data = bimp_rich)

bimp_rich_mod16 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp * Age, data = bimp_rich)

bimp_rich_mod17 <- lm(`specnumber(lab_feat_table_t)` ~ No_pollen_spp, data = bimp_rich)

bimp_rich_mod18 <- lmer(`specnumber(lab_feat_table_t)` ~ No_pollen_spp * Age + (1|Colony) + (1|Emerge_start_date), 
                       data = bimp_rich)

AICtab(bimp_rich_mod, 
       bimp_rich_mod1,
       bimp_rich_mod2,
       bimp_rich_mod3,
       bimp_rich_mod4,
       bimp_rich_mod5, 
       bimp_rich_mod6,
       bimp_rich_mod6int,
       bimp_rich_mod7,
       bimp_rich_mod8,
       bimp_rich_mod9,
       bimp_rich_mod10,
       bimp_rich_mod11,
       bimp_rich_mod12,
       bimp_rich_mod13,
       bimp_rich_mod14,
       bimp_rich_mod15,
       bimp_rich_mod16,
       bimp_rich_mod17,
       bimp_rich_mod18) # 1 is best

# dAIC df
# bimp_rich_mod1     0.0 18 -> use for diet comp
# bimp_rich_mod8     7.7 17
# bimp_rich_mod3    33.1 11
# bimp_rich_mod10   33.7 10 -> use for diet div 
# bimp_rich_mod2    34.3 11
# bimp_rich_mod     35.1 12
# bimp_rich_mod12   42.2 7 
# bimp_rich_mod11   43.5 7 
# bimp_rich_mod18   43.6 9 
# bimp_rich_mod9    44.2 8 
# bimp_rich_mod16   59.8 7 
# bimp_rich_mod13   60.1 6 
# bimp_rich_mod14   60.2 5 
# bimp_rich_mod17   61.1 4 
# bimp_rich_mod15   61.6 5 
# bimp_rich_mod6int 63.4 15
# bimp_rich_mod5    67.0 9 
# bimp_rich_mod4    67.1 10
# bimp_rich_mod7    68.1 8 
# bimp_rich_mod6    68.8 9 

simulateResiduals(bimp_rich_mod1, plot = T) # looks good
Anova(bimp_rich_mod1, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: specnumber(lab_feat_table_t)
# Chisq Df Pr(>Chisq)   
# (Intercept)  3.9920  1   0.045717 * 
#   diet         8.4391  6   0.207668   
# Age          2.0154  1   0.155708   
# Wing         0.7857  1   0.375390   
# diet:Age    19.5411  6   0.003341 **

summary(glht(bimp_rich_mod1, linfct = mcp(Age = "Tukey")))

simulateResiduals(bimp_rich_mod10, plot = T) # looks good 
Anova(bimp_rich_mod10, type = "III")
# Chisq Df Pr(>Chisq)  
# (Intercept)       2.5699  1    0.10891  
# No_pollen_spp     4.8936  2    0.08657 .
# Age               1.4830  1    0.22331  
# Wing              3.1170  1    0.07748 .
# No_pollen_spp:Age 5.9256  2    0.05168 .

# test effect of diet within age variables
joint_tests(bimp_rich_mod1, by= "Age")

# plot model for diet composition

bimp_rich_mod1 <- lmer(`specnumber(lab_feat_table_t)` ~ 
                         diet_initials * Age + Wing + (1|Colony) + (1|Emerge_start_date), 
                       data = bimp_rich)
bimp_rich_mod_means <- emmeans(bimp_rich_mod1, ~diet_initials*Age) 
bimp.rich.means.to.plot <- as.data.frame(summary(bimp_rich_mod_means))
bimp.rich.means.to.plot$tfupper <- bimp.rich.means.to.plot$emmean + bimp.rich.means.to.plot$SE
bimp.rich.means.to.plot$tlower <- bimp.rich.means.to.plot$emmean - bimp.rich.means.to.plot$SE
bimp.rich.means.to.plot$diet_initials2 <- factor(c("D", "H", "S", "DH", "SH", "DS", "DSH"))
bimp.rich.means.to.plot$diet_initials2 <- factor(bimp.rich.means.to.plot$diet_initials2, 
                                                    levels = c("D", "H", "S", "DH", "SH", "DS", "DSH"))
# bar plot - currently what I'm using in the ms 

bimp_rich_mod_barplot<-
  ggplot(bimp.rich.means.to.plot, aes(x = diet_initials2, fill = Age, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_bar(stat ="identity", 
           aes(fill = Age), 
           width = 0.6, 
           position = position_dodge(width = .6), 
           color = "black") +
  geom_errorbar(width=0.1, 
                position = position_dodge(width = .6), 
                lwd = .75) +
  theme(text = element_text(size=25),
        axis.title.x = element_blank(),
        axis.text = element_text(size=25,
                                 color = "black")) +
  scale_fill_manual(values = c("mediumorchid",
                               "darkorchid4")) +
  theme(legend.position = "none")  +
 # ggtitle(expression(paste(italic("B. impatiens")))) +
  ylab(" ") +
  xlab(" \n Diet Treatment") +
  ylim(0,60) + 
  geom_text(x = 1, hjust = 0, y = 54, label = expression(paste("Diet: ", italic("P"), "= 0.208")), size = 6) +
  geom_text(x = 1, hjust = 0, y = 49, label = expression(paste("Age: ", italic("P"), "= 0.156")), size = 6) +
  geom_text(x = 1, hjust = 0, y = 44, label = expression(paste("Diet x Age: ", italic("P"), "= 0.003")), size = 6) +
  geom_signif(y_position = 60,
              xmin = 0.75,
              xmax = 3.25,
              annotation = "1 sp.",
              tip_length = 0.05,
              textsize = 6,
              color = "black") +
  geom_signif(y_position = 60,
              xmin = 3.75,
              xmax = 6.25,
              annotation = "2 spp.",
              tip_length = 0.05,
              textsize = 6,
              color = "black") +
  geom_signif(y_position = 60,
              xmin = 6.55,
              xmax = 7.45,
              annotation = "3 spp.",
              tip_length = 0.05,
              textsize = 6,
              color = "black")
bimp_rich_mod_barplot

### PLOT MODEL WITH DIET DIVERSITY 

bimp_rich_div_mod_means <- emmeans(bimp_rich_mod10, ~No_pollen_spp*Age) 
bimp.rich.means.to.plot2 <- as.data.frame(summary(bimp_rich_div_mod_means))
bimp.rich.means.to.plot2$tfupper <- bimp.rich.means.to.plot2$emmean + bimp.rich.means.to.plot2$SE
bimp.rich.means.to.plot2$tlower <- bimp.rich.means.to.plot2$emmean - bimp.rich.means.to.plot2$SE

bimp_rich_div_mod_barplot<-
  ggplot(bimp.rich.means.to.plot2, aes(x = No_pollen_spp, fill = Age, y = emmean, ymin = tlower, ymax = tfupper)) + 
  theme_classic()+
  geom_bar(stat ="identity", 
           aes(fill = Age), 
           width = 0.6, 
           position = position_dodge(width = .6), 
           color = "black") +
  geom_errorbar(width=0.1, 
                position = position_dodge(width = .6), 
                lwd = .75) +
  theme(text = element_text(size=25),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=25, color = "black")) +
  scale_fill_manual(values = c("mediumorchid",
                               "darkorchid4")) +
  theme(legend.position = c(.75,.85))  +
 # ggtitle(expression(paste(italic("B. impatiens")))) +
  ylab(" ") +
  xlab(" \n Diet Diversity (# species)") +
  ylim(0,60) +
  geom_text(x = 0.75, hjust = 0, y = 54, label = expression(paste("Diet: ", italic("P"), "= 0.087")), size = 6) +
  geom_text(x = 0.75, hjust = 0, y = 49, label = expression(paste("Age: ", italic("P"), "= 0.223")), size = 6) +
  geom_text(x = 0.75, hjust = 0, y = 44, label = expression(paste("Diet x Age: ", italic("P"), "= 0.052")), size = 6)

bimp_rich_div_mod_barplot
