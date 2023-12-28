
# Pollen Diet Diversity does not Affect Gut Bacterial Communities or Melanization in a Social and Solitary Bee Species

# Alison E. Fowler · Quinn S. McFrederick · Lynn S. Adler

# Includes analysis of percent melanization, pollen consumption, survival by pollen diet treatment in Bombus impatiens 

# Libraries 
library(lattice)
library(psych)
library(car)
library(AER)
library(MASS)
library(corrplot)
library(AICcmodavg)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(mgcv)
library(scatterplot3d)
library(bbmle)
library(emmeans)
library(ggplot2)
library(tidyverse)
library(ggalt)
library(ggsignif)
library(multcomp)
library(coxme) 
library(survival) 
library(DHARMa)
library(aod)
library(pscl) 
library(writexl)
library(glmmTMB)
library(survminer)
library(fitdistrplus)

# Percent melanization 
## load & set up data

setwd("/Users/alisonfowler/Docs/ch4-diet-diversity")
data <- read_csv("final-data-and-scripts/final-data/melanization/melanization_data_full.csv")

data$Treat_composition <- as.factor(data$Treat_composition)
data$Treat_no_spp <- as.factor(data$Treat_no_spp)
data$Colony <- as.factor(data$Colony)
data$Date_start <- as.factor(data$Date_start)

# relevel diets 
levels(data$Treat_composition)
data$Treat_composition <- relevel(data$Treat_composition, "Dandelion + Sumac + Hawthorn")
data$Treat_composition <- relevel(data$Treat_composition, "Dandelion + Sumac")
data$Treat_composition <- relevel(data$Treat_composition, "Sumac + Hawthorn")
data$Treat_composition <- relevel(data$Treat_composition, "Dandelion + Hawthorn")
data$Treat_composition <- relevel(data$Treat_composition, "Sumac")
data$Treat_composition <- relevel(data$Treat_composition, "Hawthorn")
data$Treat_composition <- relevel(data$Treat_composition, "Dandelion")
levels(data$Treat_composition)

# make a new col with initials of diets for easier plotting 
data$diet_initials <- NA

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Dandelion") {
    data$diet_initials[i] <- "D"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Sumac") {
    data$diet_initials[i] <- "S"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Hawthorn") {
    data$diet_initials[i] <- "H"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Sumac + Hawthorn") {
    data$diet_initials[i] <- "S + H"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Dandelion + Hawthorn") {
    data$diet_initials[i] <- "D + H"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Dandelion + Sumac") {
    data$diet_initials[i] <- "D + S"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

for(i in 1:nrow(data)){
  if(data$Treat_composition[i] == "Dandelion + Sumac + Hawthorn") {
    data$diet_initials[i] <- "D + S + H"
  }  else {
    data$diet_initials[i] <- data$diet_initials[i]
  }
} 

# relevel 
data$diet_initials <- as.factor(data$diet_initials)
data$diet_initials <- relevel(data$diet_initials, "D + S + H")
data$diet_initials <- relevel(data$diet_initials, "D + S")
data$diet_initials <- relevel(data$diet_initials, "S + H")
data$diet_initials <- relevel(data$diet_initials, "D + H")
data$diet_initials <- relevel(data$diet_initials, "S")
data$diet_initials <- relevel(data$diet_initials, "H")
data$diet_initials <- relevel(data$diet_initials, "D")

## explore data

boxplot(perc_area_halved ~ Treat_composition, data = data)
boxplot(perc_area_halved ~ Treat_no_spp, data = data)
boxplot(perc_area_halved ~ Colony, data = data)
boxplot(perc_area_halved ~ Date_start, data = data)
plot(perc_area_halved ~ Wing, data = data)

boxplot(perc_area ~ Treat_composition, data = data)
boxplot(perc_area ~ Treat_no_spp, data = data)
boxplot(perc_area ~ Colony, data = data)
boxplot(perc_area ~ Date_start, data = data)
plot(perc_area ~ Wing, data = data)

hist(data$perc_area_halved)
# let's figure out what distribution to use for percent area halved 
# data bound by zero, non-integer. Gamma or beta? 
data_noNAs <- data %>% filter(perc_area_halved != is.na(perc_area_halved))
plotdist(data_noNAs$perc_area_halved, histo = TRUE, demp = TRUE)
descdist(data_noNAs$perc_area_halved, discrete=F, boot=500)
# this suggests beta 

## models 
# global model 
mod <- glmmTMB(perc_area_halved_decimal ~ Treat_composition + Wing + (1|Date_start) + (1|Colony), data = data, family = beta_family())
simresid <- simulateResiduals(mod, plot = T)

# exclude colony RE 
mod2 <- glmmTMB(perc_area_halved_decimal ~ Treat_composition + Wing + (1|Date_start), data = data, family = beta_family())

# exclude date RE 
mod3 <- glmmTMB(perc_area_halved_decimal ~ Treat_composition + Wing + (1|Colony), data = data, family = beta_family())

# exclude both 
mod4 <- glmmTMB(perc_area_halved_decimal ~ Treat_composition + Wing, data = data, family = beta_family())

AICtab(mod, mod2, mod3, mod4) 
# dAIC df
# mod   0.0 11
# mod3  7.1 10
# mod2  7.8 10
# mod4 13.0 9 

# global mod does best 

# remove wing from global model 
mod5 <- glmmTMB(perc_area_halved_decimal ~ Treat_composition + (1|Date_start) + (1|Colony), data = data, family = beta_family())

AICtab(mod, mod5) 
# dAIC df
# mod5  0.0 10
# mod   8.1 11

# model without wing does better 
simresid <- simulateResiduals(mod5, plot = T) # looks good

Anova(mod5)
# Response: perc_area_halved_decimal
# Chisq Df Pr(>Chisq)
# Treat_composition 2.6931  6     0.8463

summary(mod5)
# Family: beta  ( logit )
# Formula:          perc_area_halved_decimal ~ Treat_composition + (1 | Date_start) +      (1 | Colony)
# Data: data
# 
# AIC      BIC   logLik deviance df.resid 
# -283.0   -256.8    151.5   -303.0       92 
# 
# Random effects:
#   
#   Conditional model:
#   Groups     Name        Variance Std.Dev.
# Date_start (Intercept) 0.1946   0.4411  
# Colony     (Intercept) 0.1698   0.4120  
# Number of obs: 102, groups:  Date_start, 5; Colony, 3
# 
# Dispersion parameter for beta family (): 8.14 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                   -2.29895    0.41737  -5.508 3.62e-08 ***
#   Treat_compositionHawthorn                      0.34586    0.38295   0.903    0.366    
# Treat_compositionSumac                         0.10356    0.35966   0.288    0.773    
# Treat_compositionDandelion + Hawthorn          0.01700    0.38299   0.044    0.965    
# Treat_compositionSumac + Hawthorn             -0.25799    0.37921  -0.680    0.496    
# Treat_compositionDandelion + Sumac            -0.03744    0.39148  -0.096    0.924    
# Treat_compositionDandelion + Sumac + Hawthorn -0.05091    0.30788  -0.165    0.869    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(glht(mod5,linfct=mcp(Treat_composition="Tukey")))

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: glmmTMB(formula = perc_area_halved_decimal ~ Treat_composition + 
#                (1 | Date_start) + (1 | Colony), data = data, family = beta_family(), 
#              ziformula = ~0, dispformula = ~1)
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# Hawthorn - Dandelion == 0                                 0.34586    0.38295   0.903    0.971
# Sumac - Dandelion == 0                                    0.10356    0.35966   0.288    1.000
# Dandelion + Hawthorn - Dandelion == 0                     0.01700    0.38299   0.044    1.000
# Sumac + Hawthorn - Dandelion == 0                        -0.25799    0.37921  -0.680    0.993
# Dandelion + Sumac - Dandelion == 0                       -0.03744    0.39148  -0.096    1.000
# Dandelion + Sumac + Hawthorn - Dandelion == 0            -0.05091    0.30788  -0.165    1.000
# Sumac - Hawthorn == 0                                    -0.24231    0.37412  -0.648    0.995
# Dandelion + Hawthorn - Hawthorn == 0                     -0.32887    0.39335  -0.836    0.980
# Sumac + Hawthorn - Hawthorn == 0                         -0.60386    0.39630  -1.524    0.724
# Dandelion + Sumac - Hawthorn == 0                        -0.38330    0.40443  -0.948    0.963
# Dandelion + Sumac + Hawthorn - Hawthorn == 0             -0.39678    0.32119  -1.235    0.876
# Dandelion + Hawthorn - Sumac == 0                        -0.08656    0.35645  -0.243    1.000
# Sumac + Hawthorn - Sumac == 0                            -0.36155    0.36448  -0.992    0.954
# Dandelion + Sumac - Sumac == 0                           -0.14100    0.39015  -0.361    1.000
# Dandelion + Sumac + Hawthorn - Sumac == 0                -0.15447    0.28506  -0.542    0.998
# Sumac + Hawthorn - Dandelion + Hawthorn == 0             -0.27499    0.37462  -0.734    0.990
# Dandelion + Sumac - Dandelion + Hawthorn == 0            -0.05444    0.40231  -0.135    1.000
# Dandelion + Sumac + Hawthorn - Dandelion + Hawthorn == 0 -0.06791    0.30362  -0.224    1.000
# Dandelion + Sumac - Sumac + Hawthorn == 0                 0.22055    0.40243   0.548    0.998
# Dandelion + Sumac + Hawthorn - Sumac + Hawthorn == 0      0.20708    0.31189   0.664    0.994
# Dandelion + Sumac + Hawthorn - Dandelion + Sumac == 0    -0.01347    0.33686  -0.040    1.000
# (Adjusted p values reported -- single-step method)

# repeat with number of pollen spp as predictor 

# global model 
mod1.1 <- glmmTMB(perc_area_halved_decimal ~ Treat_no_spp + Wing + (1|Date_start) + (1|Colony), data = data, family = beta_family())

# exclude colony RE 
mod2.1 <- glmmTMB(perc_area_halved_decimal ~ Treat_no_spp + Wing + (1|Date_start), data = data, family = beta_family())

# exclude date RE 
mod3.1 <- glmmTMB(perc_area_halved_decimal ~ Treat_no_spp + Wing + (1|Colony), data = data, family = beta_family())

# exclude both 
mod4.1 <- glmmTMB(perc_area_halved_decimal ~ Treat_no_spp + Wing, data = data, family = beta_family())

# remove wing
mod5.1 <- glmmTMB(perc_area_halved_decimal ~ Treat_no_spp + (1|Date_start) + (1|Colony), data = data, family = beta_family())

AICtab(mod1.1, mod2.1, mod3.1, mod4.1, mod5.1)
# dAIC df
# mod5.1  0.0 6 
# mod1.1  7.9 7 
# mod3.1 16.5 6 
# mod2.1 16.7 6 
# mod4.1 24.6 5 

Anova(mod5.1)
# Response: perc_area_halved_decimal
# Chisq Df Pr(>Chisq)
# Treat_no_spp 1.2193  2     0.5436

## plots

### boxplot raw means  

plot<-
  ggplot(data, aes(x = diet_initials, y = perc_area_halved)) +
  geom_boxplot(aes(fill = diet_initials), position = position_dodge(1)) + 
  #  ggtitle(expression(paste(italic("M. rotundata")))) + 
  scale_fill_manual(values = c("lightskyblue", 
                               "lightskyblue",
                               "lightskyblue",
                               "dodgerblue",
                               "dodgerblue",
                               "dodgerblue",
                               "dodgerblue4")) + 
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=18)) + 
  ylab("Percent melanization") +
  xlab("Diet Treatment") + 
  ylim(0,50) +
  geom_signif(y_position = 48, 
              xmin = 0.75, 
              xmax = 3.25,
              annotation = "One Species", 
              tip_length = 0.05,
              textsize = 4,
              #color = "lightskyblue"
  ) +
  geom_signif(y_position = 48, 
              xmin = 3.75, 
              xmax = 6.25,
              annotation = "Two Species", 
              tip_length = 0.05,
              textsize = 4,
              #color = "dodgerblue"
  ) +
  geom_signif(y_position = 48, 
              xmin = 6.55, 
              xmax = 7.45,
              annotation = "Three Species", 
              tip_length = 0.05,
              textsize = 4,
              #color = "dodgerblue4"
  ) 

plot


### plot with model estimates 

# run with diet initials column for the plot 
mod5 <- glmmTMB(perc_area_halved_decimal ~ diet_initials + (1|Date_start) + (1|Colony), data = data, family = beta_family())
# extract means
means <- emmeans(mod5, ~ diet_initials, type = "response")
means.to.plot<-as.data.frame(summary(means))
means.to.plot$upperSE <- means.to.plot$response + means.to.plot$SE
means.to.plot$lowerSE <- means.to.plot$response - means.to.plot$SE
means.to.plot


# plot model values -- what we're using for the manuscript figure 
mod_plot<- 
  ggplot(means.to.plot,aes(x=diet_initials, y=response, ymin=lowerSE, ymax=upperSE)) + 
  theme_classic() +
  geom_bar(stat="identity",
           color = "black",
           fill = "gray",
           width = 0.6) +
  geom_errorbar(width = 0.065, 
                color = "black",
                lwd = 0.75) +
  theme(legend.position = "none",
        text = element_text(size = 15)) + 
  ylab("Proportion melanized \n") +
  xlab("\n Diet Treatment") + 
  # labs(fill = "Diet Treatment") + 
  ylim(0,0.25) +
  geom_signif(y_position = 0.215, 
              xmin = 0.75, 
              xmax = 3.25,
              annotation = "1 sp.", 
              tip_length = 0.05,
              textsize = 5,
              color = "black"
  ) +
  geom_signif(y_position = 0.215, 
              xmin = 3.75, 
              xmax = 6.25,
              annotation = "2 spp.", 
              tip_length = 0.05,
              textsize = 5,
              color = "black"
  ) +
  geom_signif(y_position = 0.215, 
              xmin = 6.55, 
              xmax = 7.45,
              annotation = "3 spp.", 
              tip_length = 0.05,
              textsize = 5,
              color = "black"
  ) 

mod_plot


# Categorical response variable - low/med/high melanization

## Make the 3 categories 

# make a new col
data_noNAs$mel_level3 <- NA

# less than 10 =  low 
for(i in 1:nrow(data_noNAs)){
  if(data_noNAs$perc_area[i] < 10) {
    data_noNAs$mel_level3[i] <- "Low"
  }  else {
    data_noNAs$mel_level3[i] <- data_noNAs$mel_level3[i] 
  }
} 

# between 10-26 = medium 
for(i in 1:nrow(data_noNAs)){
  if(data_noNAs$perc_area[i] >=10 & data_noNAs$perc_area[i] < 26) {
    data_noNAs$mel_level3[i] <- "Medium"
  }  else {
    data_noNAs$mel_level3[i] <- data_noNAs$mel_level3[i] 
  }
} 

# higher than 26 = high 
for(i in 1:nrow(data_noNAs)){
  if(data_noNAs$perc_area[i] >= 26) {
    data_noNAs$mel_level3[i] <- "High"
  }  else {
    data_noNAs$mel_level3[i] <- data_noNAs$mel_level3[i] 
  }
} 

summary(data_noNAs$mel_level3)


## I made a new spreadsheet in excel with a frequency table and then made a spreadsheet with the proportions of total

freq_table <- read_csv("final-data-and-scripts/final-data/melanization/mel_levels_props_initials_3cat_final.csv")

freq_table_long <- freq_table %>%  
  pivot_longer(!Diet, names_to = "mel_level", values_to = "proportion") 

freq_table_long$Diet <- as.factor(freq_table_long$Diet)
freq_table_long$mel_level <- as.factor(freq_table_long$mel_level)

freq_table_long$Diet <- relevel(freq_table_long$Diet, "D + S + H")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "D + S")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "S + H")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "D + H")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "S")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "H")
freq_table_long$Diet <- relevel(freq_table_long$Diet, "D")

freq_table_long$mel_level <- relevel(freq_table_long$mel_level, "High")
freq_table_long$mel_level <- relevel(freq_table_long$mel_level, "Medium")
freq_table_long$mel_level <- relevel(freq_table_long$mel_level, "Low")


## melanization categories plot 

categories_plot<- 
  ggplot(freq_table_long, aes(x = Diet, fill = mel_level, y = proportion)) + 
  geom_bar(aes(fill = mel_level),
           stat = "identity",
           #  position = position_dodge(width = 0.6),
           color = "black",
           width = 0.6) + 
  theme_classic() + 
  scale_fill_manual(values = c("green4",
                               "chartreuse3",
                               "darkolivegreen1")) + 
  ylab("Proportion of samples \n") + 
  xlab("\n Diet Treatment") + 
  labs(fill = "Melanization") +
  theme(text = element_text(size = 15),
        legend.position = "none") 

categories_plot


# Consumption ####

cons <- read_csv("final-data-and-scripts/final-data/melanization/melanization_consumption.csv")
cons_full <- merge(cons, data, by = "Bee_ID")

## Evaporation 

## Import evaporation data 

evap <- read_csv("final-data-and-scripts/final-data/melanization/melanization_evaporation.csv")
evap$Treat <- as.factor(evap$Treat)

levels(evap$Treat)
evap$Treat <- relevel(evap$Treat, "Dandelion + Sumac + Hawthorn")
evap$Treat <- relevel(evap$Treat, "Dandelion + Sumac")
evap$Treat <- relevel(evap$Treat, "Sumac + Hawthorn")
evap$Treat <- relevel(evap$Treat, "Dandelion + Hawthorn")
evap$Treat <- relevel(evap$Treat, "Sumac")
evap$Treat <- relevel(evap$Treat, "Hawthorn")
evap$Treat <- relevel(evap$Treat, "Dandelion")
levels(evap$Treat)

boxplot(evap$net_pollen ~evap$Treat)

hist(evap$net_pollen)

## Setting up the data 

# separate the  treatments. 
evap_D <- evap %>% filter(Treat == "Dandelion")
evap_H <- evap %>% filter(Treat == "Hawthorn")
evap_S <- evap %>% filter(Treat == "Sumac")
evap_DH <- evap %>% filter(Treat == "Dandelion + Hawthorn")
evap_SH <- evap %>% filter(Treat == "Sumac + Hawthorn")
evap_DS <- evap %>% filter(Treat == "Dandelion + Sumac")
evap_DSH <- evap %>% filter(Treat == "Dandelion + Sumac + Hawthorn")

## Calculations for each diet treatment

### Dandelion (includes instructions)

#Now let's make a linear regression with the old (final) & new (initial) weights with intercept constrained to zero for just the dandelion treatment. 
evap_D_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_D)

# now take the actual values from the dandelion-fed bees 
cons_D <- cons_full %>% filter(Treat_composition == "Dandelion")

# now predict what those values should be after accounting for evaporation: 
cons_D$predicted_new <- predict(evap_D_lg, newdata = cons_D, type = "response")

#So now we use these predicted values as the new initial pollen weights. Now we subtract old (final) weights from this to get amount consumed.  
cons_D <- cons_D %>% mutate(new_consumed = predicted_new - final_pollen_wt)


hist(cons_D$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_D)
plot(net_pollen ~ new_consumed, data = cons_D)


### Hawthorn 

evap_H_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_H)

cons_H <- cons_full %>% filter(Treat_composition == "Hawthorn")

cons_H$predicted_new <- predict(evap_H_lg, newdata = cons_H, type = "response")

cons_H <- cons_H %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_H$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_H)
plot(net_pollen ~ new_consumed, data = cons_H)


### Sumac

evap_S_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_S)

cons_S <- cons_full %>% filter(Treat_composition == "Sumac")

cons_S$predicted_new <- predict(evap_S_lg, newdata = cons_S, type = "response")

cons_S <- cons_S %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_S$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_S)
plot(net_pollen ~ new_consumed, data = cons_S)


### D + H

evap_DH_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_DH)

cons_DH <- cons_full %>% filter(Treat_composition == "Dandelion + Hawthorn")

cons_DH$predicted_new <- predict(evap_DH_lg, newdata = cons_DH, type = "response")

cons_DH <- cons_DH %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_DH$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_DH)
plot(net_pollen ~ new_consumed, data = cons_DH)


### S + H

evap_SH_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_SH)

cons_SH <- cons_full %>% filter(Treat_composition == "Sumac + Hawthorn")

cons_SH$predicted_new <- predict(evap_SH_lg, newdata = cons_SH, type = "response")

cons_SH <- cons_SH %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_SH$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_SH)
plot(net_pollen ~ new_consumed, data = cons_SH)


### D + S

evap_DS_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_DS)

cons_DS <- cons_full %>% filter(Treat_composition == "Dandelion + Sumac")

cons_DS$predicted_new <- predict(evap_DS_lg, newdata = cons_DS, type = "response")

cons_DS <- cons_DS %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_DS$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_DS)
plot(net_pollen ~ new_consumed, data = cons_DS)

### D + S + H

evap_DSH_lg <- lm(final_pollen_wt ~ initial_pollen_wt + 0, data = evap_DSH)

cons_DSH <- cons_full %>% filter(Treat_composition == "Dandelion + Sumac + Hawthorn")

cons_DSH$predicted_new <- predict(evap_DSH_lg, newdata = cons_DSH, type = "response")

cons_DSH <- cons_DSH %>% mutate(new_consumed = predicted_new - final_pollen_wt)

hist(cons_DSH$new_consumed)
plot(predicted_new ~ initial_pollen_wt, data = cons_DSH)
plot(net_pollen ~ new_consumed, data = cons_DSH)


## COMBINE ALL TREATMENTS 
cons_full_predicts <- rbind(cons_D, cons_H, cons_S, cons_DH, cons_SH, cons_DS, cons_DSH)

## Now explore new data 
boxplot(new_consumed ~ Treat_composition, data = cons_full_predicts)
boxplot(new_consumed ~ Colony, data = cons_full_predicts)
plot(new_consumed ~ Wing, data = cons_full_predicts)
hist(cons_full_predicts$new_consumed)

## Model consumption 
cons_mod <- lmer(new_consumed ~ Treat_composition + Wing + (1|Colony) + (1|Date_start), data = cons_full_predicts)
cons_mod2 <- lmer(new_consumed ~ Treat_composition + Wing + (1|Colony), data = cons_full_predicts)
cons_mod3 <- lmer(new_consumed ~ Treat_composition + Wing + (1|Date_start), data = cons_full_predicts)
cons_mod4 <- lm(new_consumed ~ Treat_composition + Wing, data = cons_full_predicts)
cons_mod5 <- lm(new_consumed ~ Treat_composition, data = cons_full_predicts)
cons_mod6 <- lmer(new_consumed ~ Treat_composition + (1|Date_start), data = cons_full_predicts)

AICtab(cons_mod, cons_mod2, cons_mod3, cons_mod4, cons_mod5, cons_mod6) # 5 best 

#           dAIC  df
# cons_mod5   0.0 8 
# cons_mod6  61.9 9 
# cons_mod4 129.1 9 
# cons_mod3 192.0 10
# cons_mod  194.0 11
# cons_mod2 203.5 10

simulateResiduals(cons_mod5, plot = T)
Anova(cons_mod5)

#Anova Table (Type II tests)

# Response: new_consumed
#                      Sum Sq  Df F value Pr(>F)
# Treat_composition 0.0013782   6  1.7114 0.1239
# Residuals         0.0162394 121   

summary(cons_mod5)
summary(glht(cons_mod5,linfct=mcp(Treat_composition="Tukey")))

# Survival ######

# make a survival object
event <- data$Dead
time<-data$Days_til_dead
survival<-Surv(time, event)

# model survival as the response 
mort<- coxme(survival ~ Treat_composition + (1|Colony), data = data) 
Anova(mort)
# Analysis of Deviance Table (Type II tests)
# 
# Response: survival
# Df  Chisq Pr(>Chisq)
# Treat_composition  6 0.3555     0.9992
