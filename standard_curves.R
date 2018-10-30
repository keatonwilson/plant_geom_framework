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

#generating linear models for plant protein
nested_protein_curves = plant_protein_curves %>%
  group_by(run_id) %>%
  nest()

#custom lm function to feed into map
protein_lm = function(df) {
  lm(standard_protein ~ abs, data = df)
}

#iterating across the nested data
models = nested_protein_curves %>%
  mutate(models = map(data, protein_lm))

get_intercept = function(mod) tidy(mod)$estimate[1]
get_slope = function(mod) tidy(mod)$estimate[2]
get_rsq = function(mod) glance(mod)$r.squared

models = models %>%
  mutate(r.squared = map_dbl(models, get_rsq),
         intercept = map_dbl(models, get_intercept),
         slope = map_dbl(models, get_slope))

#joining this data (including the models) onto the protein samples data frame
protein_master = left_join(plant_protein_samples, models, by = "run_id")

#predicting protein values
protein_master_sub = protein_master %>%
  mutate(protein = abs*slope + intercept) %>%
  select(sample_id, run_id, plant_id, age, protein) %>%
  group_by(plant_id, age) %>%
  summarize(mean_protein = mean(protein, na.rm = TRUE))

plant_master %>%
  select(-carb, -protein) %>%
  left_join(protein_master_sub) %>%
  mutate(age = as.factor(age),
         species = as.factor(species)) %>%
  ggplot(aes(x = species, y = mean_protein)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_wrap(~ age) +
  theme_classic()

plant_master = plant_master %>%
  #select(-carb, -protein) %>%
  left_join(protein_master_sub) %>%
  mutate(age = as.factor(age),
         species = as.factor(species)) %>%
  filter(!is.na(mean_protein)) %>%
  group_by(plant_id, age, species) %>%
  summarize(mean_protein = mean(mean_protein))

plant_master = ungroup(plant_master)

#messing around with models - probably something nested
library(lme4)
library(nlme)

plant_master = plant_master %>%
  mutate(plant_id = as.factor(plant_id))

lm_protein = lm(mean_protein ~ species*age, data = plant_master)

anova(lm_protein)

plant_master %>%
  ggplot(aes(x = species, y = mean_protein, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = age), position = position_dodge(width = 0.75), alpha = 0.5) +
  theme_classic()
