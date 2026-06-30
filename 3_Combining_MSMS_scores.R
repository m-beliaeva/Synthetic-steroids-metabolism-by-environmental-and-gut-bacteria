##### This code reads data with cosine and F1 MSMS scores and create a new file #####

rm(list = ls())

### Import libraries ----
library(tidyverse)
library(scales)
library(data.table)

#### PATHS AND METADATA ----
wd_path <- getwd()
setwd(wd_path)

# Use data MergeION to get MSMS matches (putative metabolites)
match_path <- paste0(wd_path, "/MergeION/Output/")

# Import parent compounds
parents <- read_csv("Reference_compounds_metadata.csv")

# Get list of output files 
output_files <- list.files(
  path = match_path,
  pattern = "\\.csv$") 

output_files <- gsub("(_F1|_cosine)\\.csv$", "", output_files) %>% unique()

# Total steroid-bacteria combinations with MSMS matches
combination <- output_files %>%
  as.data.frame() %>%
  rename(files = 1) %>%
  filter(!grepl("_met", files)) %>% 
  mutate(files = gsub("_parent_similarity", "", files)) %>% 
  distinct() %>% 
  separate(files, into = c("compound", "species"), sep = "_") %>% 
  write_csv(paste0(wd_path, "/Metabolites_collection/3_Steroid_bacteria_combinations_with_MSMS_matches.csv"))

unique(combination$compound) #18 compounds
unique(combination$species) #11 species

#### HERE METACODE STARTS ####

for (i in output_files) {
  
  # Cosine similarity table
  msms_cosine_raw <- read_csv(paste0(match_path, i, "_cosine.csv")) %>% 
    select(c(1,2,9,22)) %>% 
    rename("COSINE_SCORE" = SCORE_MERGEION) %>% 
    separate(FILENAME, into = c("ID", "COMPOUND", "SPECIES", "TIME"), sep = "_") %>% 
    separate(TIME, into = "TIME", sep = "h") %>% 
    select(-c("ID"))
  
  # If any compound appears more than one time, keep only higher score value
  # Add neutral mass
  msms_cosine <- msms_cosine_raw %>% 
    mutate(TIME = as.numeric(TIME)) %>% 
    group_by(PEPMASS, RT, COMPOUND, SPECIES, TIME) %>%
    filter(row_number(desc(COSINE_SCORE)) == 1) %>%
    ungroup() %>%
    mutate(NEUTRAL_MASS = PEPMASS - 1.0073) %>% 
    select(3,4,5,7,1,2,6)    
  
  # F1 similarity table
  msms_F1_raw <- read_csv(paste0(match_path, i, "_F1.csv")) %>% 
    select(c(1,2,9,22)) %>% 
    rename("F1_SCORE" = SCORE_MERGEION) %>% 
    separate(FILENAME, into = c("ID", "COMPOUND", "SPECIES", "TIME"), sep = "_") %>% 
    separate(TIME, into = "TIME", sep = "h") %>% 
    select(-c("ID"))
  
  msms_F1 <- msms_F1_raw %>%
    mutate(TIME = as.numeric(TIME)) %>% 
    group_by(PEPMASS, RT, COMPOUND, SPECIES, TIME) %>%
    filter(row_number(desc(F1_SCORE)) == 1) %>%
    ungroup() %>%
    mutate(NEUTRAL_MASS = PEPMASS - 1.0073) %>% 
    select(3,4,5,7,1,2,6)  
  
  save_score <- left_join(msms_cosine, msms_F1, 
                          by = c("PEPMASS", "RT", "COMPOUND", 
                                 "SPECIES", "TIME","NEUTRAL_MASS"))
  
  write_csv(save_score, paste0(wd_path, "/Metabolites_collection/3_Combined_scores/", i, "_msms_score.csv"))
  
}

#### MSMS score full histograms ----
file_dir <- paste0(wd_path, "/Metabolites_collection/3_Combined_scores/")
files <- list.files(file_dir, pattern = "_msms_score\\.csv$", full.names = TRUE)
df_all <- bind_rows(lapply(files, read.csv))

df_all %>%
  filter(F1_SCORE !=2) %>% 
  ggplot(aes(x = F1_SCORE)) +
  geom_histogram(binwidth = 0.01, fill = "#FFD700", color = "black") +
  labs(title = "Histogram of F1_SCORE Values",
       x = "F1_SCORE",
       y = "Frequency") +
  theme_minimal()

ggsave(filename = paste0(wd_path, "/Metabolites_collection/3_Full_histogram_F1_scores.pdf"),
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

ggsave(filename = paste0(wd_path, "/Metabolites_collection/3_Full_histogram_cosine_scores.pdf"),
       width = 6,
       height = 3,
       unit = "in",
       dpi = 300)
