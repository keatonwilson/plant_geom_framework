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

#Some of these are not good - let's not use any samples with an R^2 less than 0.

#joining this data (including the models) onto the protein samples data frame
protein_master = left_join(plant_protein_samples, models, by = "run_id")

#Getting rid of bad runs (#16 runID above - horrible R^2)
protein_master = protein_master %>% 
  filter(r.squared > 0.6)

#Might be able to do this with predict(), particularly if we end up using some non-linear standard curves
#predicting protein values
#
#Example from the web

# exampleTable <- data.frame(
# x = c(1:5, 1:5),
#  y = c((1:5) + rnorm(5), 2*(5:1)),
#  groups = rep(LETTERS[1:2], each = 5)
#)

#exampleTable %>% 
#  group_by(groups) %>%
#  nest() %>% 
#  mutate(model = data %>% map(~lm(y ~ x, data = .))) %>% 
#  mutate(Pred = map2(model, data, predict)) %>% 
#  unnest(Pred, data)

#Loop to do predictions on a row by row basis
protein_pred = c()
for (i in 1:length(protein_master$sample_id)) {
  protein_pred[i] = predict(protein_master$models[[i]], newdata = protein_master[i,])
}

protein_master_sub = protein_master %>%
  bind_cols(protein = protein_pred) %>%
  group_by(plant_id, age) %>%
  summarize(mean_protein = median(protein, na.rm = TRUE))

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
  #filter(!is.na(mean_protein)) %>%
  group_by(plant_id, age, species) %>%
  summarize(mean_protein = median(mean_protein))

plant_master = ungroup(plant_master)

plant_master = plant_master %>%
  mutate(protein_weight = (((mean_protein/60)*1000)/50)*1000,
         protein_percent = (protein_weight/1000)/20)


#Carb stuff
plant_carb_curves

#let's visually inspect the standard curves
ggplot(plant_carb_curves, aes(x = abs, y = standard_carbs)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ run_id)

#generating linear models for plant protein
nested_carb_curves = plant_carb_curves %>%
  filter(!is.na(abs)) %>%
  group_by(run_id) %>%
  nest()

#custom lm function to feed into map
carb_lm = function(df) {
  lm(standard_carbs ~ abs, data = df)
}

#iterating across the nested data
models_carb = nested_carb_curves %>%
  mutate(models = map(data, carb_lm))

get_rsq = function(mod) glance(mod)$r.squared

models_carb = models_carb %>%
  mutate(r.squared = map_dbl(models, get_rsq))

#joining this data (including the models) onto the protein samples data frame
carb_master = left_join(plant_carb_samples, models_carb, by = "run_id")

#Filtering out less than 0.7 R^2
carb_master = carb_master %>%
  filter(r.squared > 0.7)

carb_master = carb_master %>%
  filter(!is.na(abs))
#Loop to do predictions on a row by row basis
carb_pred = c()
for (i in 1:length(carb_master$sample_id)) {
  carb_pred[i] = predict(carb_master$models[[i]], newdata = carb_master[i,])
}

carb_master_sub = carb_master %>%
  bind_cols(carb = carb_pred) %>%
  group_by(plant_id, age) %>%
  summarize(mean_carb = median(carb, na.rm = TRUE)) %>%
  ungroup()

carb_master_sub$age = str_replace_all(carb_master_sub$age, c("old" = "O", "young" = "Y"))
carb_master_sub = carb_master_sub %>%
  filter(!is.na(age)) %>%
  mutate(age = as.factor(age)) %>%
  mutate(plant_id = as.factor(plant_id))

#warnings()
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
#  

carb_master_sub = carb_master_sub %>%
  mutate(carb_weight = (mean_carb/15)*1000,
         carb_percent = (carb_weight/1000)/20)


plant_master = plant_master %>%
  mutate(plant_id = as.factor(plant_id)) %>%
  #select(-carb, -protein) %>%
  left_join(carb_master_sub)

plant_master_gs = gs_read(plant, ws = 1) #turning it into a dataframe/tibble

plant_master_master = plant_master_gs %>%
  mutate(plant_id = as.character(plant_id)) %>%
  left_join(plant_master)

#Cleaning up
library(tidyverse)
library(lubridate)

plant_master_master = plant_master_master %>%
  mutate(plant_id = as.factor(plant_id),
         tp_date = ymd(tp_date),
         age = as.factor(age), 
         species = as.factor(species)) %>%
  select(-protein, -carb)

glimpse(plant_master_master)
summary(plant_master_master)
