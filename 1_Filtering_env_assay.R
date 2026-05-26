### This code works from the parent directory, goes into every subfolder, 
### and creates lists for data-dependent MS/MS acquisition  (csv)
### It is important that data is homogeneous and the names are the same withing the folder

#### Import libraries ----
library(tidyverse)
library(fs)

# Data directory
data_path <- paste0(getwd(), "/Metaboanalyst_output_env/")
# Control directory
ctrl_path <- paste0(getwd(), "/Metaboanalyst_output_env/DMSO-control/")
# Output directories
met_path <- paste0(getwd(), "/Filtered_metabolites/")
lists_path <- paste0(getwd(), "/MSMS_include_lists/")
# Working directory
wd_path <- getwd()
setwd(wd_path)

#### FUNCTIONS ----

# This function makes a list of the sub directories withing the parent folder
get_sub_dirs <- function(path) {
  all_dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  sub_dirs <- dirname(all_dirs)
  leaf_dirs <- setdiff(all_dirs, sub_dirs)
  return(leaf_dirs)
}

# This function makes a list of the terminal directories withing the parent folder
get_leaf_dirs <- function(path) {
  all_dirs <- list.dirs(path, recursive = TRUE, full.names = TRUE)
  sub_dirs <- dirname(all_dirs)
  leaf_dirs <- setdiff(all_dirs, sub_dirs)
  return(leaf_dirs)
}

# As control and data are analysed in separate MZmine  batches, 
# same molecules can have slight deviations in RT and mz.
# Define tolerance values
mz_tolerance <- 0.001
RT_tolerance <- 0.1

#### PREPARE CONTROL DATA ----

### READ and MERGE DATA (DMSO) 
# Collect metabolites that appear over the time course in the control (DMSO) sample
# For each species separately

dmso_dirs <- get_sub_dirs(ctrl_path)

for (i in dmso_dirs) {
  folder <- i
  words <- strsplit(folder, "/")[[1]]
  species <- tail(words, n = 1)
  index <- substr(species, 1, 2)
  
  # Get metabolites from each time point and merge
  dmso_3vs0 <- read_csv(paste0(ctrl_path, species, "/volcano_dmso_3vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  name3 <- paste("dmso_3vs0", index, sep = "_")
  assign(name3, dmso_3vs0)
  
  dmso_9vs0 <- read_csv(paste0(ctrl_path, species, "/volcano_dmso_9vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  name9 <- paste("dmso_9vs0", index, sep = "_")
  assign(name9, dmso_9vs0)
  
  dmso_27vs0 <- read_csv(paste0(ctrl_path, species, "/volcano_dmso_27vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  name27 <- paste("dmso_27vs0", index, sep = "_")
  assign(name27, dmso_27vs0)
  
  dmso <- full_join(dmso_3vs0, dmso_9vs0, copy = TRUE) %>%
    full_join(., dmso_27vs0, copy = TRUE) %>% 
    distinct(mz, RT, .keep_all = TRUE)
  name <- paste("dmso", index, sep = "_")
  assign(name, dmso)
  
  # Collect ALL metabolites that are detected in the control (DMSO) sample
  dmso_all <- read_csv(paste0(ctrl_path, species, "/dmso_all_metabolites.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(c("mz","RT", "p.value", "-log10(p)")) %>% 
    mutate(across(1:4, as.numeric))
  name_all <- paste("dmso_all", index, sep = "_")
  assign(name_all, dmso_all)
  rm(dmso_all)
  rm(dmso_3vs0)
  rm(dmso_9vs0)
  rm(dmso_27vs0)
}

#### READ DATA IN EACH SUBFOLDER, FILTER AND SAVE ----
# Create a list of locations where data is located

data_dirs <- get_leaf_dirs(data_path)

# Remove current folder and DMSO folders
data_dirs <- data_dirs[-1]
indices <- grep("DMSO", data_dirs)
data_dirs <- data_dirs[-indices]

for (i in data_dirs) {
  folder <- i
  words <- strsplit(folder, "/")[[1]]
  compound <- tail(words, n = 2)[1]
  species <- tail(words, n = 1)
  index <- substr(species, 1, 2)
  if (index == "Ct") {
    dmso_3vs0 <- dmso_3vs0_Ct
    dmso_9vs0 <- dmso_9vs0_Ct
    dmso_27vs0 <- dmso_27vs0_Ct
    dmso_all <- dmso_all_Ct
  } else if (index == "Pa") {
    dmso_3vs0 <- dmso_3vs0_Pa
    dmso_9vs0 <- dmso_9vs0_Pa
    dmso_27vs0 <- dmso_27vs0_Pa
    dmso_all <- dmso_all_Pa
  } else if (index == "Re") {
    dmso_3vs0 <- dmso_3vs0_Re
    dmso_9vs0 <- dmso_9vs0_Re
    dmso_27vs0 <- dmso_27vs0_Re
    dmso_all <- dmso_all_Re
  } else if (index == "Sh") {
    dmso_3vs0 <- dmso_3vs0_Sh
    dmso_9vs0 <- dmso_9vs0_Sh
    dmso_27vs0 <- dmso_27vs0_Sh
    dmso_all <- dmso_all_Sh
  }
  
  ## Read and merge data  
  met_3vs0 <- read_csv(paste0(folder, "/volcano_3vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  
  met_9vs0 <- read_csv(paste0(folder, "/volcano_9vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  
  met_27vs0 <- read_csv(paste0(folder, "/volcano_27vs0.csv")) %>% 
    separate(col = "...1", into = c("id", "mz", "RT"), sep = "/", remove = TRUE) %>%
    separate(col = "mz", into = c("mz"), sep = "mz", remove = TRUE) %>%
    separate(col = "RT", into = c("RT"), sep = "min", remove = TRUE) %>%
    select(-id) %>% 
    mutate(across(1:6, as.numeric))
  
  met <- full_join(met_3vs0, met_9vs0, copy = TRUE) %>% 
    full_join(., met_27vs0, copy = TRUE) %>% 
    distinct(mz, RT, .keep_all = TRUE)
  
  # Select metabolites with positive FC and subtract DMSO control
  met_3vs0_pos <- met_3vs0 %>% 
    filter(log2(FC) > 0) 
  met_9vs0_pos <- met_9vs0 %>% 
    filter(log2(FC) > 0)
  met_27vs0_pos <- met_27vs0 %>% 
    filter(log2(FC) > 0)
  met_pos <- full_join(met_3vs0_pos, met_9vs0_pos, copy = TRUE) %>% 
    full_join(., met_27vs0_pos, copy = TRUE) %>%
    distinct(mz, RT, .keep_all = TRUE)
  
  dmso_3vs0_pos <- dmso_3vs0 %>% 
    filter(log2(FC) > 0) 
  dmso_9vs0_pos <- dmso_9vs0 %>% 
    filter(log2(FC) > 0)
  dmso_27vs0_pos <- dmso_27vs0 %>% 
    filter(log2(FC) > 0)
  dmso_pos <- full_join(dmso_3vs0_pos, dmso_9vs0_pos, copy = TRUE) %>% 
    full_join(., dmso_27vs0_pos, copy = TRUE) %>%
    distinct(mz, RT, .keep_all = TRUE)
  
  # Filter dmso metabolites out of met_pos based on mz and RT tolerances
  met_pos_wo_dmso <- met_pos %>%
    filter(!mz %in% (dmso_pos$mz - mz_tolerance) & !mz %in% (dmso_pos$mz + mz_tolerance),
           !RT %in% (dmso_pos$RT - RT_tolerance) & !RT %in% (dmso_pos$RT + RT_tolerance)) 
  
  # Filter negative values of log2(FC) and remove ALL metabolites that are found in DMSO samples
  met_3vs0_neg <- met_3vs0 %>% 
    filter(log2(FC) < 0)
  met_9vs0_neg <- met_9vs0 %>% 
    filter(log2(FC) < 0)
  met_27vs0_neg <- met_27vs0 %>% 
    filter(log2(FC) < 0)
  met_neg <- full_join(met_3vs0_neg, met_9vs0_neg, copy = TRUE) %>%
    full_join(., met_27vs0_neg, copy = TRUE) %>%
    distinct(mz, RT, .keep_all = TRUE)
  
  # Filter dmso metabolites out of met_neg based on mz and RT tolerances
  met_neg_wo_dmso <- met_neg %>%
    filter(!mz %in% (dmso_all$mz - mz_tolerance) & !mz %in% (dmso_all$mz + mz_tolerance),
           !RT %in% (dmso_all$RT - RT_tolerance) & !RT %in% (dmso_all$RT + RT_tolerance)) 
  
  # Merge positive and negative selection lists
  full_list <- full_join(met_pos_wo_dmso, met_neg_wo_dmso, copy = TRUE) %>% 
    distinct(mz, RT, .keep_all = TRUE)  %>% 
    filter(1 < RT & RT < 5.5) %>% 
    filter(mz < 900) %>% 
    write_csv(paste0(met_path, compound, "_", species, ".csv"))
  
  # Create an MS/MS inclusion list
  # Parameters may be changed depending on the MS/MS settings
  include_list <- full_list %>% 
    rename(`Prec. m/z` = mz, `Ret. Time (min)` = RT) %>% 
    mutate(`Delta m/z (ppm)` = 20,
           `Z` = 1,
           `Prec. Type` = "Preferred",
           `Delta Ret. Time (min)` = 0.4,
           `Iso. Width` = "Narrow (~1.3 m/z)",
           `Collision Energy` = 10) %>% 
    select(c(1,7,8,9,2,10,11,12)) %>% 
    write_csv(paste0(lists_path, compound, "_", species, ".csv"))
}
