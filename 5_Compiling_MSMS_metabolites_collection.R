rm(list = ls())

### Import libraries ----
library(tidyverse)
library(scales)
library(data.table)
library(ggplot2)
library(tools)
library(rstatix)

#### SET THRESHOLDS ----
# Tolerance for MS
abs_mass = 0.0025 # absolute mz tolerance for mass_diff
rel_mass = 25*10^-6 # relative mz tolerance (25 ppm)
abs_rt = 0.15 # absolute RT tolerance

#### DATA PATHS ----
# Set paths to working directory and data folders
wd_path <- getwd()
data_path <- paste0(wd_path, "/Metabolites_collection/4_Metabolite_lists")
setwd(wd_path)

#### Combine all metabolite lists into one ---
filenames <- list.files(path = data_path, pattern = "\\.csv$", full.names = TRUE)
print(filenames)

col_specs <- cols(
  Compound = col_character(),
  species = col_character(),
  metabolite = col_character(),
  MZ = col_double(),
  RT = col_double(),
  COSINE_SCORE = col_double(),
  F1_SCORE = col_double(),
  row_number = col_double())

#### Bind data ----
# Read individual files and create df with all metabolites and parent compounds
df_all <- filenames %>%
  lapply(read_csv, col_types = col_specs) %>%   
  bind_rows() %>% 
  select(-row_number) %>% 
  mutate(MZ = round(MZ,4),
         RT = round(RT, 2)) %>% 
  distinct()

#### MSMS score histograms ----
df_all %>%
  filter(F1_SCORE !=2) %>% 
  ggplot(aes(x = F1_SCORE)) +
  geom_histogram(binwidth = 0.01, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of F1_SCORE Values",
       x = "F1_SCORE",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/5_Histogram_F1_scores.png"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)

df_all %>%
  filter(COSINE_SCORE !=2) %>% 
  ggplot(aes(x = COSINE_SCORE)) +
  geom_histogram(binwidth = 0.01, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of COSINE_SCORE Values",
       x = "COSINE_SCORE",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/5_Histogram_cosine_scores.png"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)

#### Add parent data and check unique metabolites ----
# Read parent data; use H+ adduct instead of Na+ adduct for DexP
parent_data <- read.csv(paste0(wd_path, '/Reference_compounds_metadata.csv')) %>% 
  filter(PARENT_ID != "Dexamethasone-sodium-phosphate") %>% 
  mutate(PARENT_ID = gsub('Dexamethasone-sodium-phosphate_met', 'Dexamethasone-sodium-phosphate', PARENT_ID)) %>% 
  rename(MZ_PARENT = MZ,
         RT_PARENT = RT,
         MW_PARENT = MW,
         SMILES_PARENT = SMILES) %>% 
  filter(!str_detect(PARENT_ID, "_met")) %>% 
  select(-MS_DATA)

# Assign parent entries by MZ and RT tolerance
df_mass_diff <- left_join(df_all, parent_data, by = c("Compound" = "PARENT_ID")) %>% 
  select(2,1,7,5,6,3,4,11,9,10,8) %>% 
  mutate(delta_mass = MZ-MZ_PARENT,
         parent_match = ifelse((abs(delta_mass) <= abs_mass & abs(RT_PARENT - RT) <= abs_rt) 
                               | (abs(delta_mass) <= abs(rel_mass * MZ_PARENT) & abs(RT_PARENT - RT) <= abs_rt), 
                               "yes", NA)) %>% 
  select(-c("metabolite", "COSINE_SCORE", "F1_SCORE", "MW_PARENT")) %>%
  distinct()

df_mass_diff %>% filter(is.na(parent_match)) %>% nrow() #129 non-unique metabolites detected excluding parent hits
df_mass_diff %>% filter(is.na(parent_match)) %>% select(Compound) %>% distinct %>% nrow() #16 compounds
df_mass_diff %>% filter(is.na(parent_match)) %>% select(species) %>% distinct %>% nrow() #11 species

df_combinations <- df_mass_diff %>% 
  filter(is.na(parent_match)) %>% 
  select(Compound, species) %>% 
  distinct %>% 
  write_csv(paste0(wd_path, '/Metabolites_collection/5_Steroid_bacteria_combinations_with_MSMS_matches_filtered.csv'))

df_mass_diff$Compound <- gsub('Dexamethasone-sodium-phosphate', 'Dexamethasone-phosphate', df_mass_diff$Compound)
write_csv(df_mass_diff, paste0(wd_path, '/Metabolites_collection/5_Metabolites_filtered.csv'))

#### Mass difference histograms ----
df_mass_diff %>%
  filter(is.na(parent_match)) %>% 
  ggplot(aes(x = delta_mass)) +
  geom_histogram(binwidth = 1, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of delta_mass Values",
       x = "delta_mass",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/5_Mass_difference_histogram.png"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)

df_mass_diff %>%
  filter(is.na(parent_match)) %>% 
  filter(abs(delta_mass) <= 100) %>% 
  ggplot(aes(x = delta_mass)) +
  geom_histogram(binwidth = 1, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of delta_mass Values",
       x = "delta_mass",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/5_Mass_difference_histogram_abs100.png"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)

df_mass_diff %>%
  filter(is.na(parent_match)) %>% 
  filter(abs(delta_mass) <= 50) %>% 
  ggplot(aes(x = delta_mass)) +
  geom_histogram(binwidth = 1, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of delta_mass Values",
       x = "delta_mass",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/5_Mass_difference_histogram_abs50.png"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)
