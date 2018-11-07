#Initial Figures and Analyses for Plant Geometric Framework Project
#Keaton Wilson
#2018-1107
#keatonwilson@me.com


#Some of these values are extreme - a value of 40 here corresponds to ~ 65% protein by dry weight. These are definitely outliers. Let's remove them. 

plant_master = plant_master %>%
  filter(mean_protein < 25 & mean_protein > 0) %>%
  mutate(protein_weight = (((mean_protein/60)*1000)/50)*1000,
         protein_percent = (protein_weight/1000)/20)

#messing around with models - probably something nested
library(lme4)
library(nlme)

plant_master = plant_master %>%
  mutate(plant_id = as.factor(plant_id))

lmm_1 = lme(protein_weight ~ species + age, random = ~1|plant_id, na.action = na.omit, data = plant_master)
summary(lmm_1)

lm_protein = lm(mean_protein ~ species + age, data = plant_master)

anova(lm_protein)

plant_master %>%
  filter(mean_protein > 0) %>%
  ggplot(aes(x = species, y = protein_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()

#What happens if we combine?

plant_master %>%
  ggplot(aes(x = species, y = mean_protein)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic()


plant_master %>%
  filter(carb_percent < 0.75) %>%
  ggplot(aes(x = species, y = carb_percent, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()

#Nutrient space figure
library(car)
install.packages("hexbin")
library(hexbin)

data_sub = plant_master %>%
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

plant_master %>%
  filter(carb_percent < 0.75) %>%
  ggplot(aes(x = protein_percent, y = carb_percent, fill = species)) +
  stat_density_2d(geom = "polygon", alpha = 0.3) +
  facet_wrap(~ species) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_classic()

