rm(list = ls())

### Import libraries ----
library(tidyverse)
library(scales)
library(data.table)
library(ggplot2)
library(tools)
library(rstatix)
library(viridis)

# For plotting
scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))
}

### SET THRESHOLDS ----
# Tolerance for MS
abs_mass = 0.002 # absolute mz tolerance
rel_mass = 20*10^-6 # relative mz tolerance (20 ppm)
abs_rt = 0.15 # absolute RT tolerance
#rel_rt = 0.02 # relative RT tolerance

# MS baseline noise
ms_noise <- 1

### DATA PATHS ----
# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# Read metabolite data
metabolites <- read.csv(paste0(wd_path, "/Manually_detected_metabolites.csv"))

# Folder with untargeted MS data from MZMINE
mzmine_path_env <- paste0(dirname(wd_path), "/2_Untargeted_analysis/Environmental_bacteria/")
mzmine_path_gut <- paste0(dirname(wd_path), "/2_Untargeted_analysis/Gut_bacteria/")

### PREPARATION ----
# Make list of species
species_env <- c("Ct", "Re", "Pa", "Sh")
species_gut <- c("Bu", "Cs", "Er", "Ef", "Cc", "Lr", "Pv")

#### HERE METACODE STARTS ####

for (j in 1:nrow(metabolites)) {
  
  #### 1. Retrieve MSMS and untargeted data ----
  
  ##### 1.1. Extract information from metabolite entry ----
  compound_id <- metabolites$Compound[j]
  bact_name <- metabolites$species[j]
  bact <- gsub("\\.", "", metabolites$species[j])
  index <- substr(bact, 1, 2)
  
  ##### 1.2. Get untargeted data ----
  # Import data here
  # Depending on env or gut bacteria, select folder
  if (index %in% species_env) {
    
    # We skip non-existing data 
    if (!file.exists(paste0(mzmine_path_env, compound_id, "/", bact, "_metaboanalyst.csv"))) {
      message("File does not exist:", paste0(mzmine_path_env, compound_id, "/", bact, "_metaboanalyst.csv"))
      next
    }
    
    mzmine_raw <- read_csv(paste0(mzmine_path_env, compound_id, "/", bact, "_metaboanalyst.csv"))
  } else if (index %in% species_gut) {
    
    # We skip non-existing data 
    if (!file.exists(paste0(mzmine_path_gut, compound_id, "/", bact, "_metaboanalyst.csv"))) {
      message("File does not exist:", paste0(mzmine_path_gut, compound_id, "/", bact, "_metaboanalyst.csv"))
      next
    }
    
    mzmine_raw <- read_csv(paste0(mzmine_path_gut, compound_id, "/", bact, "_metaboanalyst.csv"))
  }
  
  # Fix column names if needed
  names(mzmine_raw) <- gsub("Detamethasone", "Betamethasone", names(mzmine_raw))
  names(mzmine_raw) <- gsub("Declomethasone", "Beclomethasone", names(mzmine_raw))
  names(mzmine_raw) <- gsub("_dipropionate", "-dipropionate", names(mzmine_raw))
  names(mzmine_raw) <- gsub("_acetate", "-acetate", names(mzmine_raw))
  names(mzmine_raw) <- gsub("_sodium", "-sodium", names(mzmine_raw))
  names(mzmine_raw) <- gsub("_phosphate", "-phosphate", names(mzmine_raw))
  names(mzmine_raw) <- gsub("R.aeruginosa", "P.aeruginosa", names(mzmine_raw))
  names(mzmine_raw) <- gsub("PH4", "pH4", names(mzmine_raw))
  names(mzmine_raw) <- gsub("PH7", "pH7", names(mzmine_raw))
  
  mzmine_reshaped <- mzmine_raw %>%
    filter(Filename != "time") %>% 
    separate(col = "Filename", into = c("mzmine_id", "MZ", "RT"), sep = "/", remove = TRUE) %>%
    mutate(MZ = str_replace(MZ, "mz", "")) %>% 
    mutate(RT = str_replace(RT, "min", "")) %>% 
    mutate_all(as.numeric) %>% 
    pivot_longer(cols = colnames(.)[4:length(colnames(.))]) %>% 
    mutate("timepoint" = as.factor(as.numeric(str_extract(name, "(?<=\\D)(\\d+)(?=h)")))) %>% 
    select(c(4, 6, 1, 2, 3, 5)) %>% 
    filter(RT > 1, RT < 5.5) %>% 
    separate(col = name, into = c("id", "Compound", "species", "time", "rep"), sep = "_") %>% 
    separate(col = rep, into = c("rep"), sep = ".mzML") %>% 
    select(-id, -time, -rep, -mzmine_id)
  
  ##### 1.3. Retrieve areas from untargeted data for metabolites ----
  retrieved_data <- list()
  ion_data_mass <- metabolites$MZ[j]
  ion_data_rt <- metabolites$RT[j]
  
  metabolite_area <- mzmine_reshaped %>% 
    mutate(Mass_diff = MZ - ion_data_mass, 
           RT_diff = RT - ion_data_rt) %>% 
    filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*ion_data_mass,
           abs(RT_diff) <= abs_rt)
  
  parent_data_mass <- metabolites$MZ_PARENT[j]
  parent_data_rt <- metabolites$RT_PARENT[j]
  
  if (nrow(metabolite_area)==0) {
    next
  }
  
  # Get parent data and merge
  parent_area <- mzmine_reshaped %>% 
    mutate(Mass_diff = MZ - parent_data_mass, 
           RT_diff = RT - parent_data_rt) %>% 
    filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*parent_data_mass,
           abs(RT_diff) <= abs_rt)
  
  df1 <- rbind(metabolite_area, parent_area)
  
  ## Do final adjustments to df1 and keep in df2
  # Impute ms_noise as a ms area background instead of 0 where area not detected
  df2 <- df1 %>% 
    select(-c(Mass_diff, RT_diff)) %>% 
    mutate(value = ifelse(is.na(value), ms_noise, value)) %>% 
    mutate(met = paste0(MZ, "_", RT))
  
  # Check point
  # How many entries for each compound for each timepoint in df2 (must be 4 ~ 4 replicates)
  rep_count <- df2 %>%
    group_by(MZ, RT, timepoint) %>%
    summarize(entries = n())
  
  any(rep_count$entries != 4) # (must be FALSE: all 4 replicates)
  
  #### 2. Retrieve pH control data ----
  # Add controls from  MZmine output (Metaboanalyst format) 
  
  ##### 2.1. Get untargeted mzmine data ----
  # Then find metabolites from dfs (what we are interested in) in the pH control
  if (index %in% species_env) {
    mzmine_raw_pH <- read_csv(paste0(mzmine_path_env, compound_id, "/", "pH_metaboanalyst.csv"))
  } else if (index %in% species_gut) {
    mzmine_raw_pH <- read_csv(paste0(mzmine_path_gut, compound_id, "/", "pH_metaboanalyst.csv"))
  }
  
  # Fix column names if needed
  names(mzmine_raw_pH) <- gsub("Detamethasone", "Betamethasone", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("Declomethasone", "Beclomethasone", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("_dipropionate", "-dipropionate", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("_acetate", "-acetate", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("_sodium", "-sodium", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("_phosphate", "-phosphate", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("R.aeruginosa", "P.aeruginosa", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("PH4", "pH4", names(mzmine_raw_pH))
  names(mzmine_raw_pH) <- gsub("PH7", "pH7", names(mzmine_raw_pH))
  
  # Split the first column and reshape
  mzmine_reshaped_pH <- mzmine_raw_pH %>%
    filter(Filename != "time") %>% 
    separate(col = "Filename", into = c("mzmine_id", "MZ", "RT"), sep = "/", remove = TRUE) %>%
    mutate(MZ = str_replace(MZ, "mz", "")) %>% 
    mutate(RT = str_replace(RT, "min", "")) %>% 
    mutate_all(as.numeric) %>% 
    pivot_longer(cols = colnames(.)[4:length(colnames(.))]) %>% 
    mutate("timepoint" = as.factor(as.numeric(str_extract(name, "(?<=\\D)(\\d+)(?=h)")))) %>% 
    select(c(4, 6, 1, 2, 3, 5)) %>% 
    filter(RT > 1, RT < 5.5) %>% 
    separate(col = name, into = c("id", "Compound", "species", "time", "rep"), sep = "_") %>% 
    separate(col = rep, into = c("rep"), sep = ".mzML") %>% 
    select(-id, -time, -rep, -mzmine_id)
  
  ##### 3.2. Search over MZ and RT in pH data ----
  # Impute MZ and RT as in parent compounds
  df3 <- df2 %>% 
    select(Compound, species, MZ, RT) %>% 
    distinct()
    
  pH <- data.frame()
  df_pH <- list()
  df_pH[[1]] <- pH
  
  for(n in 1:nrow(df3)){
    ion_data_mass <- df3$MZ[n]
    ion_data_rt <- df3$RT[n]
    
    pH <- mzmine_reshaped_pH %>% 
      mutate(Mass_diff = MZ - ion_data_mass, 
             RT_diff = RT - ion_data_rt) %>% 
      filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*ion_data_mass,
             abs(RT_diff) <= abs_rt) 
    
    # If metabolites are not found in pH data
    # Impute MZ and RT as in parent compounds and ms_noise in area
    if (nrow(pH) == 0) {
      pH <- data.frame(
        MZ = rep(ion_data_mass, 8), 
        RT = rep(ion_data_rt, 8), 
        value = rep(ms_noise, 8), 
        Compound = rep(compound_id, 8), 
        species = rep(c("pH4", "pH7"), each = 4),
        Mass_diff = rep(0, 8), 
        RT_diff = rep(0, 8), 
        timepoint = rep(c(0, 3, 9, 27), 2) 
      )
    } 
    
    # If there are more matches found in pH data, get the data with closest RT
    else if (nrow(pH) > 32) {
      pH <- pH %>% 
        filter(RT == min(abs(RT))) %>% 
        mutate(MZ = ion_data_mass, RT = ion_data_rt)           
    }
    
    else {
      pH <- pH %>% 
        mutate(MZ = ion_data_mass, RT = ion_data_rt) 
    }
    
    df_pH[[n+1]] <- pH  
    
  }
  
  # Collect pH data together and calculate mean values
  pH_retrieved <- rbindlist(df_pH, fill = TRUE)
  
  if (ncol(pH_retrieved) > 0) {
    pH_retrieved <- pH_retrieved %>% 
      select(-Mass_diff, -RT_diff) %>% 
      mutate(met = paste0(MZ, "_", RT))
  }
  
  df4 <- rbind(df2, pH_retrieved) 
  
  #### 3. Plot metabolites with controls ----
  df5 <- df4 %>% 
    group_by(Compound, species, met, timepoint) %>%
    summarise(mean = mean(value), sd = sd(value)) %>% 
    mutate(timepoint = as.numeric(as.character(timepoint)))
  
  plot <- df5 %>% 
    group_by(Compound, species, met, timepoint) %>%
    ggplot(mapping = aes(x = timepoint, y = mean, ymin = 0,
                         group = interaction(met, species))) +
    geom_point(mapping = aes(shape = species), size = 2) +
    geom_line(mapping = aes(linetype = species), color = "black", size = 0.75) +
    geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), size = 0.5, width = 0.5) +
    scale_x_continuous(labels = label_number(1)) +
    scale_y_continuous(label = scientific_10) +
    facet_wrap(~ met, scales = "free_y", ncol = 1) +
    theme_bw() +
    labs (x= "time, h", y= "peak area") +
    ggtitle(paste0(compound_id, "_", bact_name)) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 8),
          strip.text.y = element_text(size = 8),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position= "bottom",
          legend.direction = "vertical") +
    guides(
      color = guide_legend(ncol = 1),  
      shape = guide_legend(ncol = 1),  
      linetype = guide_legend(ncol = 1)) 
  
  ggsave(filename = paste0(wd_path, "/Metabolites_collection/4_Metabolite_plots_curated/", compound_id, "_", bact, "_", j, ".png"),
         plot = plot,
         width = 3,
         height = 4.5,
         unit = "in",
         dpi = 300)
}


