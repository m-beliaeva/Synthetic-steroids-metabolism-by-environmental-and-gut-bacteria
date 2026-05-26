#### Import libraries ----
library(tidyverse)
library(fs)
library(tools)

# Working directory
wd_path <- getwd()
setwd(wd_path)
# Data directory
data_path <- paste0(getwd(), "/Filtered_metabolites/")

# Get all include lists
all_csvs <- list.files(data_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

# Get all metabolites
df <- all_csvs %>%
  tibble(
    file_path = .,
    file_name = file_path_sans_ext(basename(.))) %>%
  mutate(data = map(file_path, read_csv)) %>%
  unnest(data) %>% 
  select(-1) %>% 
  separate(file_name, into = c("Name", "Condition"), sep = "_") %>% 
  mutate(Condition = str_replace(Condition, "^([A-Za-z])", "\\1.")) %>% 
  mutate(Name = gsub("67-Dehydro-ethynylestradiol", "67-Dehydro-ethinylestradiol", Name)) %>% 
  write_csv(paste0(wd_path, "/Total_metabolites_list.csv"))

