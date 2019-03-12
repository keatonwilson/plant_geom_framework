#Initial Figures and Analyses for Plant Geometric Framework Project
#Keaton Wilson
#2018-1107
#keatonwilson@me.com
#
#Libraries
library(tidyverse)
library(lubridate)
library(skimr)
library(GGally)
library(gridExtra)
library(ggpubr)

#Need to pull in and run data manipulation on raw data from google sheets first
source("./scripts/standard_curves_merging.R", local = TRUE)

#Making sure the dataframe is loaded
plants = plant_master_master

#Let's address outliers
plants %>%
  dplyr::filter(plant_id == 9) %>%
  dplyr::select(carb_percent, protein_percent)

ggplot(plants, aes(x = protein_percent)) +
  geom_histogram()

ggplot(plants, aes(x = carb_percent)) +
  geom_histogram()
#
plants %>%
  dplyr::filter(protein_percent < 0 | carb_percent > 0.75 ) %>%
  dplyr::select(protein_percent, carb_percent, sample_id)

#Actually looks good after preprocessing in earlier script - but let's remove wonky stuff anyway
plants %>%
  dplyr::filter(protein_weight > 10000 | protein_weight < 0) %>%
  dplyr::select(protein_weight, mean_protein, protein_percent, plant_id, age)

#75% protein by weight is pretty high, but, we'll leave it for now. Get rid of negative value though.
plants = plants %>%
  dplyr::filter(protein_weight > 0)

ggplot(plants, aes(x = protein_percent)) +
  geom_histogram

#Much Better
#On to Carbs
ggplot(plants, aes(x = carb_percent)) +
  geom_histogram()

#More problems here. Everything above 100% doesn't make sense - same cutoff as protein 
plants %>%
  dplyr::filter(carb_percent > 0.75) %>%
  dplyr::select(plant_id, age, species, mean_carb, carb_weight, carb_percent)

plants_filtered = plants %>%
  dplyr::filter(protein_weight > 0 & protein_percent < 0.75 & carb_percent < 0.75 & carb_percent > 0)

ggplot(plants_filtered, aes(x = carb_percent)) +
  geom_histogram()

#Also generating a time period between transplant and collection (proxy for age)
plants_filtered = plants_filtered %>%
  dplyr::mutate(time_lag_interval = interval(tp_date, collect_date),
         time_lag = time_lag_interval %/% days(1)) %>%
  dplyr::select(-tp_date, -collect_date, -time_lag_interval)

plants_filtered = plants_filtered %>%
  dplyr::select(-sample_id) %>%
  distinct()

#Preliminary Visualization
#Protein
protein_percent = plants_filtered %>%
  filter(mean_protein > 0) %>%
  ggplot(aes(x = species, y = protein_percent*100, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = age), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2, size = 2) +
  theme_classic() +
  ylab("Percent Protein") +
  xlab("") +
  scale_x_discrete(labels = c("Nicotiana tabacum",
                              "Datura discolor", 
                              "Datura wrightii", 
                              "Capsicum annuum", 
                              "Probiscidea parviflora", 
                              "Lycopersicon esculentum", 
                              "Nicotiana attenuata")) +
  theme(axis.text.x = element_text(face = "italic", size = 11, angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_fill_discrete(name = "Leaf Age", 
                      labels = c("Old", "Young")) +
  annotate(geom = "text", x = 7, y = 65, label = "a", size = 8)

protein_percent

ggsave(filename = "./output/protein_percent.png", protein_percent, device = "png")
ggsave(filename = "./output/protein_percent.tiff", protein_percent, device = "tiff")

#carbs
carb_percent = plants_filtered %>%
  ggplot(aes(x = species, y = carb_percent*100, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = age), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2, size = 2) +
  theme_classic() +
  ylab("Percent Carbohydrates") +
  xlab("") +
  scale_x_discrete(labels = c("Nicotiana tabacum",
                              "Datura discolor", 
                              "Datura wrightii", 
                              "Capsicum annuum", 
                              "Probiscidea parviflora", 
                              "Lycopersicon esculentum", 
                              "Nicotiana attenuata")) +
  theme(axis.text.x = element_text(face = "italic", size = 11, angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_fill_discrete(name = "Leaf Age", 
                      labels = c("Old", "Young")) +
  annotate(geom = "text", x = 7, y = 65, label = "b", size = 8)


carb_percent

Fig_1 = ggarrange(protein_percent, carb_percent, common.legend = TRUE)
ggsave(filename = "./output/Fig_1.tiff", Fig_1, device = "tiff", width = 14, height = 8, dpi = 150)

#messing around with models - probably something nested
library(lme4)
library(nlme)
library(lmerTest)

#Just out of curiosity - effect of size and time-lag
lm_height_protein = lm(protein_percent ~ plant_height, data = plants_filtered)
lm_height_carb = lm(carb_percent ~ plant_height, data = plants_filtered)
summary(lm_height_protein)
summary(lm_height_carb)

lm_time_protein = lm(protein_percent ~ time_lag , data = plants_filtered)
lm_time_carb = lm(carb_percent ~ time_lag, data = plants_filtered)
summary(lm_time_protein)
summary(lm_time_carb)

#Time lag affects both (effect size is small, but we can include it as a covariate in the nested models)

#Serious modeling
library(car)
library(MASS)

#Look at transformation - because we want to model percent, a logit transformation is appropriate - makes our data closer to normal
qqp(plants_filtered$protein_percent, "norm")
qqp(logit(plants_filtered$protein_percent), "norm")
logit(plants_filtered$protein_percent)

#protein model
lmm_protein = lmer(logit(protein_percent) ~ species + age + time_lag + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm_protein)
anova(lmm_protein)

#Result is that there is a significant effect of age, but not species on protein. Not huge differences.

#Carbs
lmm_carbs = lmer(logit(carb_percent) ~ species + age + time_lag + (1|plant_id), data = plants_filtered, REML = FALSE)
summary(lmm_carbs)
anova(lmm_carbs)

#Main result here is that there is a significant effect of species, but not age - also time-lag is more improtant here. The older plants had more carbs.

#Correlation between protein and carbohydrates?
lm_both = lm(protein_percent ~ carb_percent, data = plants_filtered)
summary(lm_both)

#Nutrient space figure
library(car)
#install.packages("hexbin")

#Making a summary data frame 
# data_sub = plants_filtered %>%
#   filter(carb_percent < 0.75) %>%
#   group_by(species, age) %>%
#   summarize(protein_percent_median = median(protein_percent),
#             protein_percent_1st = quantile(protein_percent, probs = 0.25),
#             protein_percent_3rd = quantile(protein_percent, probs = 0.75),
#             protein_percent_max = max(protein_percent),
#             protein_percent_min = min(protein_percent),
#             carb_percent_median = median(carb_percent),
#             carb_percent_1st = quantile(carb_percent, probs = 0.25),
#             carb_percent_3rd = quantile(carb_percent, probs = 0.75),
#             carb_percent_max = max(carb_percent),
#             carb_percent_min = min(carb_percent)
#   )
# 
# ggplot(data_sub, aes(fill = species, color = species)) +
#   geom_rect(aes(xmin = protein_percent_1st,
#                 xmax = protein_percent_3rd,
#                 ymin = carb_percent_1st, 
#                 ymax = carb_percent_3rd), alpha = 0.4) +
#   #geom_segment(aes(y = carb_percent_median,
#   #                 yend = carb_percent_median,
#   #                 x = protein_percent_min,
#   #                 xend = protein_percent_max)) +
#   geom_segment(aes(x = protein_percent_median,
#                    yend = carb_percent_max,
#                    y = protein_percent_min,
#                    xend = protein_percent_median)) +
#   facet_wrap(~ age)

# nutrient_space = plants_filtered %>%
#   filter(carb_percent < 0.75) %>%
#   ggplot(aes(x = protein_percent, y = carb_percent, fill = species)) +
#   #stat_density_2d(geom = "polygon", alpha = 0.2) +
#   geom_bin2d(bins = 50) +
#   geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.4) +
#   #facet_wrap(~ species) +
#   xlim(c(-1, 1)) +
#   ylim(c(-1, 1)) +
#   coord_cartesian(xlim = c(0,0.75), ylim = c(0,0.75)) +
#   theme_classic()
# 
# ggsave(filename = "./output/nutrient_space.png", nutrient_space, device = "png")
# 
plant_labels = c("CT" = "Nicotiana tabacum",
                 "DD" = "Datura discolor", 
                 "DW" = "Datura wrightii", 
                 "PE" = "Capsicum annuum", 
                 "PP" = "Probiscidea parviflora", 
                 "TO" = "Lycopersicon esculentum", 
                 "WT" = "Nicotiana attenuata")


centroids = merge(plants_filtered, 
                  aggregate(cbind(mean.x = protein_percent, 
                                  mean.y = carb_percent
                                  ) ~ species, plants_filtered, mean),
                  by = "species")
sds = aggregate(cbind(sd.x = protein_percent,
                      sd.y = carb_percent) ~ species, plants_filtered, sd)
summed_sds = sds %>%
  mutate(summed_sds = sd.x + sd.y)

centroids = left_join(centroids, summed_sds, by = 'species')

clustering = centroids %>%
  ggplot(aes(x = protein_percent, y = carb_percent, color = species)) +
  #facet_wrap( ~ age) +
  geom_segment(aes(x = mean.x, y = mean.y, xend = protein_percent, yend = carb_percent), alpha = 0.3, size = 1) +
  geom_point(size = 2, alpha = 0.3) +
  geom_point(aes(x = mean.x, y = mean.y, size = summed_sds)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_classic() +
  scale_color_discrete(labels = c(expression(italic("Nicotiana tabacum")),
                                  expression(italic("Datura discolor")), 
                                  expression(italic("Datura wrightii")), 
                                  expression(italic("Capsicum annuum")), 
                                  expression(italic("Probiscidea parviflora")), 
                                  expression(italic("Lycopersicon esculentum")), 
                                  expression(italic("Nicotiana attenuata"))), name = "Species") +
  scale_size_continuous(name = "Summed SD") +
  xlab("Percent Protein") +
  ylab("Percent Carbohydrates") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text.align = 0)
  
ggsave(filename = "./output/Fig_2.tiff", clustering, device = "tiff", width = 14, height = 8, dpi = 150)
#Calculating p:c ratios

p_c = plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  mutate(p_c = protein_weight/carb_weight) %>%
  ggplot(aes(x = species, y = p_c, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = age), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2, size = 2) +
  theme_classic() +
  ylab("P:C") +
  xlab("") +
  scale_x_discrete(labels = c("Nicotiana tabacum",
                              "Datura discolor", 
                              "Datura wrightii", 
                              "Capsicum annuum", 
                              "Probiscidea parviflora", 
                              "Lycopersicon esculentum", 
                              "Nicotiana attenuata")) +
  theme(axis.text.x = element_text(face = "italic", size = 11, angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_discrete(name = "Leaf Age", 
                      labels = c("Old", "Young")) +
  ylim(c(0,6)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate(geom = "text", label = "a", x = 7, y = 6, size = 8)


#Total Nutrient Content
total_nutrient = plants_filtered %>%
  mutate(nutrient_percent = (carb_weight + protein_weight)/20000) %>%
  ggplot(aes(x = species, y = nutrient_percent*100, fill = age)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = age), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2, size = 2) +
  theme_classic() +
  ylab("Percent Nutrient Content") +
  xlab("") +
  scale_x_discrete(labels = c("Nicotiana tabacum",
                              "Datura discolor", 
                              "Datura wrightii", 
                              "Capsicum annuum", 
                              "Probiscidea parviflora", 
                              "Lycopersicon esculentum", 
                              "Nicotiana attenuata")) +
  theme(axis.text.x = element_text(face = "italic", size = 11, angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_discrete(name = "Leaf Age", 
                      labels = c("Old", "Young")) +
  ylim(c(0,100)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotate(geom = "text", label = "b", x = 7, y = 100, size = 8)

Fig_3 = ggarrange(p_c, total_nutrient, common.legend = TRUE)

ggsave(filename = "./output/Fig_3.tiff", Fig_3, device = "tiff", width = 14, height = 8, dpi = 150)

#Stats for figures above

fig_3_analysis = plants_filtered %>%
  filter(carb_percent < 0.75) %>%
  mutate(p_c = protein_weight/carb_weight,
         total_nutrient = (carb_weight + protein_weight)/20000)
  
#p_c model
lmm_pc = lmer(p_c ~ species + age + time_lag + (1|plant_id), 
           data = fig_3_analysis, REML = FALSE)
summary(lmm_pc)
anova(lmm_pc)

#Result is that there is a significant effect of age, but not species on protein. Not huge differences.

#Nutrient Content
lmm_nut = lmer(total_nutrient ~ species + age + time_lag + (1|plant_id), 
              data = fig_3_analysis, REML = FALSE)
summary(lmm_nut)
anova(lmm_nut)

fig_3_analysis %>%
  group_by(species) %>%
  summarize(mean_nut = mean(total_nutrient),
            sd_nut = sd(total_nutrient))
#Let's rethink the visualization. 
#Each point in the nutrient space plot above represents a single diet source for either a young or old leave on a plant. This is a measure of the slope of a rail in nutrient space. So can we plot a "rail-cloud" for each plant species?
#
library(lemon)
plant_labels = c("CT" = "Nicotiana tabacum",
                 "DD" = "Datura discolor", 
                 "DW" = "Datura wrightii", 
                 "PE" = "Capsicum annuum", 
                 "PP" = "Probiscidea parviflora", 
                 "TO" = "Lycopersicon esculentum", 
                "WT" = "Nicotiana attenuata")

nutrient_rails = plants_filtered %>%
  mutate(slope = carb_percent/protein_percent) %>%
  dplyr::select(age, species, protein_percent, carb_percent, slope) %>%
  ggplot(aes(x = protein_percent, y = carb_percent)) +
  geom_abline(aes(slope = slope, intercept = 0, lty = age), size = 0.25, alpha = 0.4) +
  facet_wrap(~species, labeller = as_labeller(plant_labels)) +
  theme_classic() +
  geom_abline(data = plants_filtered %>%
                mutate(slope = carb_percent/protein_percent) %>%
                dplyr::select(age, species, protein_percent,carb_percent, slope) %>%
                group_by(species, age) %>%
                summarize(median_slope = median(slope)),
              aes(slope = median_slope, intercept = 0, lty = age),
              size = 0.75) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  theme(strip.text = element_text(face = "italic")) +
  xlab("Protein (g)") +
  ylab("Carbohydrates (g)") +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 16), 
        axis.line = element_line()) +
  scale_linetype_discrete(name = "Age", labels = c("Old", "Young")) +
  lims(x = c(0,3), y = c(0,3)) +
  annotate(geom = "point", x= 1.587, y = 1.894, size = 4)

ggsave(filename = "./output/Fig_4.tiff", nutrient_rails, device = "tiff", dpi = 150)

#Summaries
plants_filtered %>%
  group_by(species) %>%
  summarize(mean_carb_percent = mean(carb_percent),
            sd_carb = sd(carb_percent),
            mean_protein_percent = mean(protein_percent),
            sd_protein = sd(protein_percent))
