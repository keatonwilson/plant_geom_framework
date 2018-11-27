#Initial Figures and Analyses for Plant Geometric Framework Project
#Keaton Wilson
#2018-1107
#keatonwilson@me.com
#
#Libraries
library(tidyverse)
library(lubridate)

#Need to pull in and run data manipulation on raw data from google sheets first
#source("./scripts/standard_curves_merging.R")

#Making sure the dataframe is loaded
plants = plant_master_master
plants

#Let's address outliers
ggplot(plants, aes(x = protein_percent)) +
  geom_histogram()

plants %>%
  filter(protein_weight > 10000 | protein_weight < 0) %>%
  select(protein_weight, mean_protein, protein_percent, plant_id, age)

#75% protein by weight is pretty high, but, we'll leave it for now. Get rid of negative value though.
plants = plants %>%
  filter(protein_weight > 0)

ggplot(plants, aes(x = carb_percent)) +
  geom_histogram()

#More problems here. Everything above 100% doesn't make sense. 
View(plants %>%
  filter(carb_percent > 0.75) %>%
  select(plant_id, age, species, mean_carb, carb_weight, carb_percent))

plants_filtered = plants %>%
  filter(protein_weight > 0 & protein_percent < 0.75 & carb_percent < 0.80 )

#Also generating a time period between transplant and collection (proxy for age)
plants_filtered = plants_filtered %>%
  mutate(time_lag_interval = interval(tp_date, collect_date),
         time_lag = time_lag_interval %/% days(1)) %>%
  select(-tp_date, -collect_date, -time_lag_interval)


#Preliminary Visualization

#Protein
plants_filtered %>%
  filter(mean_protein > 0) %>%
  ggplot(aes(x = species, y = protein_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()

#carbs
plants_filtered %>%
  ggplot(aes(x = species, y = carb_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()

#messing around with models - probably something nested
library(lme4)
library(nlme)

#Just out of curiosity - effect of size and time-lag
lm_height = lm(protein_percent ~ time_lag , data = plants_filtered)
summary(lm_height)

plot(protein_percent ~ time_lag, data = plants_filtered)


#Some data exploration
library(car)
library(MASS)

qqp(logit(plants_filtered$protein_percent), "norm")
logit(plants_filtered$protein_percent)

lmm = lmer(logit(protein_percent) ~ species + age + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm)
Anova(lmm)

#Example lme code from Meck's stuff
##Random slopes and intercepts - this model is actually better based on AIC
# fm2 = lme(log(CatWeight)~Time*Treatment, data = subset(manducaalive, manducaalive$Round != "9"),
#           random = ~1+Time|TrueID, na.action = na.omit)

lmm_0 = lme(protein_percent ~ 1, random = ~1|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_0)
lmm_1 = lme(protein_percent ~ species, random = ~1+species|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_1)
lmm_2 = lme(protein_percent ~ species + age, random = ~1+species|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_2)
lmm_3 = lme(protein_percent ~ species * age, random = ~1+species|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_3)

AIC(lmm_0)
AIC(lmm_1)
AIC(lmm_2)
AIC(lmm_3)

lmm0 = lmer(protein_percent ~ 1 + (1|plant_id), data = plants_filtered)
summary(lmm0)
AIC(lmm0)

lmm1 = lmer(protein_percent ~ species + (1|plant_id), data = plants_filtered)
summary(lmm1)
AIC(lmm1)

lmm_0_c = lme(carb_percent ~ 1, random = ~1|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_0_c)
lmm_1_c = lme(carb_percent ~ species, random = ~1|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_1_c)
lmm_2_c = lme(carb_percent ~ species + age, random = ~1|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_2_c)
lmm_3_c = lme(carb_percent ~ species * age, random = ~1|plant_id, na.action = na.omit, data = plants_filtered)
summary(lmm_3_c)

lmm = lmer(logit(carb_percent) ~ species + age + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm)
Anova(lmm)

#Correlation between protein and carbohydrates?
lm_both = lm(protein_percent ~ carb_percent, data = plants_filtered)
summary(lm_both)

#Nutrient space figure
library(car)
#install.packages("hexbin")

#Making a summary data frame 
data_sub = plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  group_by(species, age) %>%
  summarize(protein_percent_median = median(protein_percent),
            protein_percent_1st = quantile(protein_percent, probs = 0.25),
            protein_percent_3rd = quantile(protein_percent, probs = 0.75),
            protein_percent_max = max(protein_percent),
            protein_percent_min = min(protein_percent),
            carb_percent_median = median(carb_percent),
            carb_percent_1st = quantile(carb_percent, probs = 0.25),
            carb_percent_3rd = quantile(carb_percent, probs = 0.75),
            carb_percent_max = max(carb_percent),
            carb_percent_min = min(carb_percent)
  )

ggplot(data_sub, aes(fill = species, color = species)) +
  geom_rect(aes(xmin = protein_percent_1st,
                xmax = protein_percent_3rd,
                ymin = carb_percent_1st, 
                ymax = carb_percent_3rd), alpha = 0.4) +
  #geom_segment(aes(y = carb_percent_median,
  #                 yend = carb_percent_median,
  #                 x = protein_percent_min,
  #                 xend = protein_percent_max)) +
  geom_segment(aes(x = protein_percent_median,
                   yend = carb_percent_max,
                   y = protein_percent_min,
                   xend = protein_percent_median)) +
  facet_wrap(~ age)

plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  ggplot(aes(x = protein_percent, y = carb_percent, fill = species)) +
  stat_density_2d(geom = "polygon", alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  facet_wrap(~ species) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_classic()

#Calculating p:c ratios

plants_filtered %>%
  mutate(p_c = protein_weight/carb_weight) %>%
  ggplot(aes(x = species, y = p_c, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()