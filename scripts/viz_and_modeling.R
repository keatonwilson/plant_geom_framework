#Initial Figures and Analyses for Plant Geometric Framework Project
#Keaton Wilson
#2018-1107
#keatonwilson@me.com
#
#Libraries
library(tidyverse)
library(lubridate)
library(skimr)

#Need to pull in and run data manipulation on raw data from google sheets first
#source("./scripts/standard_curves_merging.R")

#Making sure the dataframe is loaded
plants = plant_master_master

#Let's address outliers
ggplot(plants, aes(x = protein_percent)) +
  geom_histogram()

#Actually looks good after preprocessing in earlier script - but let's remove wonky stuff anyway
plants %>%
  filter(protein_weight > 10000 | protein_weight < 0) %>%
  select(protein_weight, mean_protein, protein_percent, plant_id, age)

#75% protein by weight is pretty high, but, we'll leave it for now. Get rid of negative value though.
plants = plants %>%
  filter(protein_weight > 0)

ggplot(plants, aes(x = carb_percent)) +
  geom_histogram()

#More problems here. Everything above 100% doesn't make sense - same cutoff as protein 
View(plants %>%
  filter(carb_percent > 0.75) %>%
  select(plant_id, age, species, mean_carb, carb_weight, carb_percent))

plants_filtered = plants %>%
  filter(protein_weight > 0 & protein_percent < 0.75 & carb_percent < 0.75 & carb_percent > 0)

ggplot(plants_filtered, aes(x = carb_percent)) +
  geom_histogram()

#Also generating a time period between transplant and collection (proxy for age)
plants_filtered = plants_filtered %>%
  mutate(time_lag_interval = interval(tp_date, collect_date),
         time_lag = time_lag_interval %/% days(1)) %>%
  dplyr::select(-tp_date, -collect_date, -time_lag_interval)

#Preliminary Visualization

#Protein
protein_percent = plants_filtered %>%
  filter(mean_protein > 0) %>%
  ggplot(aes(x = species, y = protein_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.2) +
  theme_classic()

ggsave(filename = "./output/protein_percent.png", protein_percent, device = "png")

#carbs
carb_percent = plants_filtered %>%
  ggplot(aes(x = species, y = carb_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.2) +
  theme_classic()

ggsave(filename = "./output/carb_percent.png", carb_percent, device = "png")

#messing around with models - probably something nested
library(lme4)
library(nlme)

#Just out of curiosity - effect of size and time-lag
lm_height_protein = lm(protein_percent ~ time_lag , data = plants_filtered)
lm_height_carb = lm(carb_percent ~ time_lag, data = plants_filtered)
summary(lm_height_protein)
summary(lm_height_carb)

#They both have an effect, but the effect size is relatively small. We can add it as a covariate.

#Serious modeling
library(car)
library(MASS)

#Look at transformation - because we want to model percent, a logit transformation is appropriate - makes our data closer to normal
qqp(plants_filtered$protein_percent, "norm")
qqp(logit(plants_filtered$protein_percent), "norm")
logit(plants_filtered$protein_percent)

#protein model
lmm = lmer(logit(protein_percent) ~ species + age + time_lag + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm)
Anova(lmm)

#Result is that there is a significant effect of age, but not species on protein. Not huge differences.


#Carbs
lmm = lmer(logit(carb_percent) ~ species + age + time_lag + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm)
Anova(lmm)

#Main result here is that there is a significant effect of species, but not age - also time-lag is more improtant here. The older plants had more carbs.

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

nutrient_space = plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  ggplot(aes(x = protein_percent, y = carb_percent, fill = species)) +
  stat_density_2d(geom = "polygon", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.4) +
  facet_wrap(~ species) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  coord_cartesian(xlim = c(0,0.75), ylim = c(0,0.75)) +
  theme_classic() + 
  geom_abline(slope = 37/31, intercept = 0, alpha = 0.4) + #Control
geom_abline(slope = 35/18, intercept = 0, alpha = 0.4) + #Low Protein
geom_abline(slope = 28/30, intercept = 0, alpha = 0.4) + #Low Carb
geom_abline(slope = 35/21, intercept = 0, alpha = 0.4) + #Medium Protein
geom_abline(slope = (25/30), intercept = 0, alpha = 0.4) #Ultra Low Carb

ggsave(filename = "./output/nutrient_space.png", nutrient_space, device = "png")

#Calculating p:c ratios

plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  mutate(p_c = protein_weight/carb_weight) %>%
  ggplot(aes(x = species, y = p_c, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()
