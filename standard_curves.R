#Initial Script for Analyzing Data from Plant Geometric Framework Project
#Keaton Wilson
#2018-10-30
#keatonwilson@me.com


#libraries
library(tidyverse)
library(googlesheets)
library(purrr)
library(broom)

#reading in the data from googlesheets
gs_ls()

plant = gs_title("Plant Geometric Framework") #registering the sheet I want

plant_master = gs_read(plant, ws = 1) #turning it into a dataframe/tibble
plant_carb_curves = gs_read(plant, ws = 2)
plant_carb_samples = gs_read(plant, ws = 3)
plant_protein_curves = gs_read(plant, ws = 4)
plant_protein_samples = gs_read(plant, ws = 5)

#let's visually inspect the standard curves
ggplot(plant_protein_curves, aes(x = abs, y = standard_protein)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  facet_wrap(~ run_id)

#They look curvilinear (2nd order polynomial), so let's use that. 

#generating linear models for plant protein
nested_protein_curves = plant_protein_curves %>%
  filter(!is.na(abs)) %>%
  group_by(run_id) %>%
  nest()

#custom lm function to feed into map
protein_lm = function(df) {
  lm(standard_protein ~ poly(abs, 2), data = df)
}

#iterating across the nested data
models = nested_protein_curves %>%
  mutate(models = map(data, protein_lm))

get_rsq = function(mod) glance(mod)$r.squared

models = models %>%
  mutate(r.squared = map_dbl(models, get_rsq))

#joining this data (including the models) onto the protein samples data frame
protein_master = left_join(plant_protein_samples, models, by = "run_id")

#Might be able to do this with predict(), particularly if we end up using some non-linear standard curves
#predicting protein values
#
#Example from the web

exampleTable <- data.frame(
  x = c(1:5, 1:5),
  y = c((1:5) + rnorm(5), 2*(5:1)),
  groups = rep(LETTERS[1:2], each = 5)
)

exampleTable %>% 
  group_by(groups) %>%
  nest() %>% 
  mutate(model = data %>% map(~lm(y ~ x, data = .))) %>% 
  mutate(Pred = map2(model, data, predict)) %>% 
  unnest(Pred, data)

#Loop to do predictions on a row by row basis
protein_pred = c()
for (i in 1:length(protein_master$sample_id)) {
  protein_pred[i] = predict(protein_master$models[[i]], newdata = protein_master[i,])
}

protein_master_sub = protein_master %>%
  bind_cols(protein = protein_pred) %>%
  select(sample_id, run_id, plant_id, age, protein) %>%
  group_by(plant_id, age) %>%
  summarize(mean_protein = mean(protein, na.rm = TRUE))

#plant_master %>%
 # select(-carb, -protein) %>%
#  left_join(protein_master_sub) %>%
#  mutate(age = as.factor(age),
#         species = as.factor(species)) %>%
#  ggplot(aes(x = species, y = mean_protein)) +
#  geom_boxplot() +
#  geom_jitter(width = 0.1) +
#  facet_wrap(~ age) +
#  theme_classic()

plant_master = plant_master %>%
  #select(-carb, -protein) %>%
  left_join(protein_master_sub) %>%
  mutate(age = as.factor(age),
         species = as.factor(species)) %>%
  filter(!is.na(mean_protein)) %>%
  group_by(plant_id, age, species) %>%
  summarize(mean_protein = mean(mean_protein))

plant_master = ungroup(plant_master)

#Some of these values are extreme - a value of 40 here corresponds to ~ 65% protein by dry weight. These are definitely outliers. Let's remove them. 

plant_master = plant_master %>%
  filter(mean_protein < 25 & mean_protein > 0)

#messing around with models - probably something nested
library(lme4)
library(nlme)

plant_master = plant_master %>%
  mutate(plant_id = as.factor(plant_id))

lmm_1 = lme(mean_protein ~ species + age, random = ~1|plant_id, na.action = na.omit, data = plant_master)
summary(lmm_1)

lm_protein = lm(mean_protein ~ species + age, data = plant_master)

anova(lm_protein)

plant_master %>%
  filter(mean_protein > 0) %>%
  ggplot(aes(x = species, y = mean_protein, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()

#What happens if we combine?

plant_master %>%
ggplot(aes(x = species, y = mean_protein)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_classic()
