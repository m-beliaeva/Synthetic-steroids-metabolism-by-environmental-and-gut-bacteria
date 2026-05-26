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
area_thresh <- 500

# MSMS thresholds
F1_thresh <- 0.1
cosine_thresh <- 0.5

# Thresholds for produced and reduced metabolites
FC_thresh_prod <- 3 #new metabolite should appear 8x in comparison with background
FC_thresh_red <- -1 #50% of present metabolite (i.e. parent compound) should disappear 

# Adjusted pval threshold
padj_thresh <- 0.05

### DATA PATHS ----
# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# Folder with msms scores from MERGEION
msms_path <- paste0(wd_path, "/Metabolites_collection/3_Combined_scores/")

# Folder with untargeted MS data from MZMINE
mzmine_path_env <- paste0(dirname(wd_path), "/2_Untargeted_analysis/Environmental_bacteria/")
mzmine_path_gut <- paste0(dirname(wd_path), "/2_Untargeted_analysis/Gut_bacteria/")

### PREPARATION ----
# Get the list of parent compounds with RT and MZ
parents <- read_csv(paste0(wd_path, "/Reference_compounds_metadata.csv"))

# Make list of species
species_env <- c("Ct", "Re", "Pa", "Sh")
species_gut <- c("Bu", "Cs", "Er", "Ef", "Cc", "Lr", "Pv")

filenames <- list.files(path = msms_path, pattern = "\\.csv$", full.names = TRUE)
filenames <- file_path_sans_ext(basename(filenames))

print(filenames)

#### HERE METACODE STARTS ####

for (j in 1:length(filenames)) {
  filename <- filenames[j]
  words <- strsplit(filename, "_")[[1]]
  compound <- words[1]
  
  if (words[2] == "met") {
    compound_id <- paste0(compound, "_", words[2])
  } else {
    compound_id <- compound
  }
  
  bact <- words[length(words) - 4]
  index <- substr(bact, 1, 2)
  
  ms_data <- parents %>% filter(PARENT_ID == compound_id) %>% pull(MS_DATA)
  
  # Get bact name in proper format with dot
  bact_name <- sub("([A-Z])", "\\1.", bact)
  
  #### 1. Retrieve MSMS and untargeted data ----
  
  ##### 1.1. Get a list of msms matches ----
  msms_score <- read_csv(paste0(msms_path, "/", filename, ".csv")) %>% 
    rename(MZ = PEPMASS) %>% 
    distinct()
  
  msms_hits <- msms_score %>% filter(F1_SCORE >= F1_thresh | COSINE_SCORE >= cosine_thresh) #%>% 
  
  ##### 1.2. Get untargeted data ----
  
  # Import data here
  # Depending on env or gut bacteria, select folder
  if (index %in% species_env) {
    
    # We skip non-existing data 
    if (!file.exists(paste0(mzmine_path_env, ms_data, "/", bact, "_metaboanalyst.csv"))) {
      message("File does not exist:", paste0(mzmine_path_env, ms_data, "/", bact, "_metaboanalyst.csv"))
      
      next
    }
    
    mzmine_raw <- read_csv(paste0(mzmine_path_env, ms_data, "/", bact, "_metaboanalyst.csv"))
  } else if (index %in% species_gut) {
    
    # We skip non-existing data 
    if (!file.exists(paste0(mzmine_path_gut, ms_data, "/", bact, "_metaboanalyst.csv"))) {
      message("File does not exist:", paste0(mzmine_path_gut, ms_data, "/", bact, "_metaboanalyst.csv"))
      
      next
    }
    
    mzmine_raw <- read_csv(paste0(mzmine_path_gut, ms_data, "/", bact, "_metaboanalyst.csv"))
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
  
  ##### 1.3. Retrieve areas from untargeted data for msms_hits ----
  retrieved_data <- list()
  
  for(i in 1:nrow(msms_hits)){
    rm(tmp)
    ion_data <- msms_hits[i,]
    ion_data_mass <- ion_data %>% pull(MZ)
    ion_data_rt <- ion_data %>% pull(RT)
    cosine_score <- ion_data %>% pull(COSINE_SCORE)
    f1_score <- ion_data %>% pull(F1_SCORE)
    
    metabolite_area <- mzmine_reshaped %>% 
      mutate(Mass_diff = MZ - ion_data_mass, 
             RT_diff = RT - ion_data_rt) %>% 
      filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*ion_data_mass,
             abs(RT_diff) <= abs_rt)
    
    # If more than one matching RT/MZ combination found
    # Keep only those with higher area
    if (nrow(metabolite_area) > 16) {
      tmp <- metabolite_area %>% 
        group_by(Compound, species, MZ, RT, timepoint) %>% 
        mutate(value = replace_na(value, 0)) %>% 
        summarize(mean_area = mean(value)) %>% 
        filter(mean_area > area_thresh) 
      
      if (nrow(tmp) > 0) {
        mz_selected <- tmp %>% pull(MZ)
        rt_selected <- tmp %>% pull(RT)
        
        metabolite_area <- metabolite_area %>% 
          filter(MZ == mz_selected,
                 RT == rt_selected)
      }
      
    } 
    
    # Add MSMS score to retrieved ions 
    # And keep the same MZ and RT as in MSMS hits list
    metabolite_area <- metabolite_area %>% 
      mutate(COSINE_SCORE = cosine_score, 
             F1_SCORE = f1_score) %>% 
      mutate(MZ = ion_data_mass, RT = ion_data_rt)
    
    retrieved_data[[i+1]] <- metabolite_area  
    
  }
  
  # Get collected data into df1
  df1 <- rbindlist(retrieved_data, fill = TRUE) 
  anyNA(df1$COSINE_SCORE) 
  
  ## Check is there is parent data in df1
  # Get parent data
  parent_data <- parents %>% filter(PARENT_ID == compound_id) %>% 
    rename(PEPMASS = "MZ")
  
  parent_data_mass <- parent_data %>% pull(PEPMASS)
  parent_data_rt <- parent_data %>% pull(RT)
  
  # Screen df1 for parent data
  parent_area <- df1 %>% 
    mutate(Mass_diff = MZ - parent_data_mass, 
           RT_diff = RT - parent_data_rt) %>% 
    filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*parent_data_mass,
           abs(RT_diff) <= abs_rt)
  
  # If parent_area is empty (msms_hits did not capture parent compound) 
  # retrieve parent data separately
  # and impute COSINE_SCORE and F1_SCORE == 2 as an indication that its included manually
  if (nrow(parent_area) == 0) {
    parent_area <- mzmine_reshaped %>% 
      mutate(Mass_diff = MZ - parent_data_mass, 
             RT_diff = RT - parent_data_rt) %>% 
      filter(abs(Mass_diff) <= abs_mass | abs(Mass_diff) <= rel_mass*parent_data_mass,
             abs(RT_diff) <= abs_rt) %>% 
      mutate(COSINE_SCORE = 2, 
             F1_SCORE = 2)
    
    df1 <- rbind(df1, parent_area)
  }
  
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
  
  anyNA(df2$COSINE_SCORE)
  anyNA(df2$value)
  
  # Get metabolites_list
  metabolites_list1 <- df2 %>% 
    arrange(desc(F1_SCORE)) %>% 
    distinct(MZ, RT, COSINE_SCORE, F1_SCORE, met)
  
  #### 2. Calculate FC and pvals for retrieved data on msms hits ----
  
  ##### 2.1. Compute pval and padj ----
  
  df_anova <- data.frame()
  df_anova_list <- list()
  df_anova_list[[1]] <- df_anova
  
  for (m in unique(df2$met)) {
    tmp <- df2 %>% 
      filter(met == m)
    
    test <- tmp %>% 
      anova_test(value ~ timepoint)
    
    pval <- get_anova_table(test)$p
    
    df_anova <- tmp %>% 
      mutate(pval = pval)
    
    df_anova_list[[length(df_anova_list)+1]] <- df_anova
  }
  
  df_anova_pval <- bind_rows(df_anova_list)
  
  df_anova_padj <- df_anova_pval %>% 
    select(c(Compound, species, met, pval)) %>% 
    distinct()
  
  df_anova_padj$padj <- p.adjust(df_anova_padj$pval, method = "fdr")
  
  df3 <- full_join(df_anova_pval, df_anova_padj, 
                   by = c("Compound", "species", "met","pval"))
  
  ##### 2.2. Compute FC for each timepoint ----
  df4 <- df3 %>% 
    group_by(Compound, species, MZ, RT, met, COSINE_SCORE, F1_SCORE, timepoint, pval, padj) %>%
    summarise(mean = mean(value), sd = sd(value)) %>% 
    pivot_wider(names_from = timepoint, values_from = c(mean, sd)) %>% 
    mutate(FC_3vs0 = log2(mean_3/mean_0), 
           FC_9vs0 = log2(mean_9/mean_0),
           FC_27vs0 = log2(mean_27/mean_0)) %>% 
    ungroup()
  
  ##### 2.3. Collect produced and reduced metabolites by FC_thresh ----
  df4_filtered <- df4 %>%
    filter(FC_3vs0 > FC_thresh_prod | FC_9vs0 > FC_thresh_prod | FC_27vs0 > FC_thresh_prod |
             FC_3vs0 < FC_thresh_red | FC_9vs0 < FC_thresh_red | FC_27vs0 < FC_thresh_red) %>% 
    filter(padj <= padj_thresh)
  
  # List of compounds retrieved from mzmine and following pattern 
  metabolites_list2 <- df4_filtered %>%
    arrange(desc(F1_SCORE)) %>% 
    distinct(MZ, RT, COSINE_SCORE, F1_SCORE, met) 
  
  # Get raw data back for metabolites_list2
  df3_filtered <- df3 %>% filter(met %in% metabolites_list2$met) %>% 
    select(-pval, -padj)
  
  if (nrow(df3_filtered) == 0) {
    next
  }
  
  #### 3. Retrieve pH control data ----
  # Only for metabolites in metabolites_list2
  # Add controls from  MZmine output (Metaboanalyst format) 
  
  ##### 3.1. Get untargeted mzmine data ----
  # Then find metabolites from dfs (what we are interested in) in the pH control
  if (index %in% species_env) {
    mzmine_raw_pH <- read_csv(paste0(mzmine_path_env, ms_data, "/", "pH_metaboanalyst.csv"))
  } else if (index %in% species_gut) {
    mzmine_raw_pH <- read_csv(paste0(mzmine_path_gut, ms_data, "/", "pH_metaboanalyst.csv"))
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
  
  pH <- data.frame()
  df_pH <- list()
  df_pH[[1]] <- pH
  
  if (nrow(metabolites_list2) > 0) {
    for(n in 1:nrow(metabolites_list2)){
      
      ion_data <- metabolites_list2[n,]
      ion_data_mass <- ion_data %>% pull(MZ)
      ion_data_rt <- ion_data %>% pull(RT)
      cosine_score <- ion_data %>% pull(COSINE_SCORE)
      f1_score <- ion_data %>% pull(F1_SCORE)
      
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
          COSINE_SCORE = rep(cosine_score, 8), 
          F1_SCORE = rep(f1_score, 8), 
          value = rep(ms_noise, 8), 
          Compound = rep(compound, 8), 
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
          mutate(COSINE_SCORE = cosine_score, 
                 F1_SCORE = f1_score) %>% 
          mutate(MZ = ion_data_mass, RT = ion_data_rt)           
      }
      
      else {
        pH <- pH %>% mutate(COSINE_SCORE = cosine_score, 
                            F1_SCORE = f1_score) %>% 
          mutate(MZ = ion_data_mass, RT = ion_data_rt) 
      }
      
      df_pH[[n+1]] <- pH  
      
    }
  }
  
  # Collect pH data together and calculate mean values
  pH_retrieved <- rbindlist(df_pH, fill = TRUE)
  
  # Check point
  anyNA(pH_retrieved$COSINE_SCORE)
  anyNA(pH_retrieved$F1_SCORE)
  
  if (ncol(pH_retrieved) > 0) {
    pH_retrieved <- pH_retrieved %>% 
      select(-Mass_diff, -RT_diff) %>% 
      mutate(met = paste0(MZ, "_", RT))
  }
  
  #### 4. Analyze species vs pH controls ----
  
  ##### 4.1. Bind pH and species data together ----
  # If metabolites are not found (value = 0), ms_noise
  
  df5 <- rbind(df3_filtered, pH_retrieved) %>% 
    mutate_at(vars(value), ~ifelse(is.na(.), ms_noise, .)) 
  
  # Check point
  # For which metabolites we retrieved pH controls
  species_counts <- df5 %>%
    group_by(met) %>%
    summarize(unique_species = n_distinct(species))
  any(species_counts$unique_species != 3) # must be 1-3 (max) for bact, pH7 and pH4
  
  anyNA(df5$F1_SCORE)
  anyNA(df5$value)
  
  # Fix names if needed
  df5$MZ <- as.factor(df5$MZ)
  df5$RT <- as.factor(df5$RT)
  df5$species <- as.factor(df5$species)
  df5$Compound <- gsub("Detamethasone", "Betamethasone", df5$Compound)
  df5$Compound <- gsub("Declomethasone", "Beclomethasone", df5$Compound)
  df5$Compound <- gsub("_dipropionate", "-dipropionate", df5$Compound)
  df5$Compound <- gsub("_acetate", "-acetate", df5$Compound)
  df5$Compound <- gsub("_sodium", "-sodium", df5$Compound)
  df5$Compound <- gsub("_phosphate", "-phosphate", df5$Compound)
  df5$species <- gsub("R.aeruginosa", "P.aeruginosa", df5$species)
  df5$species <- gsub("PH4", "pH4", df5$species)
  df5$species <- gsub("PH7", "pH7", df5$species)
  
  ##### 4.2. Calculate pvals between pH controls and species across timepoints ----
  # Two-way interactions between condition pH7 and time for each metabolite
  met_list <- unique(df5$met)
  
  # Save results here
  df_anova_pH <- data.frame()
  df_anova_list_pH <- list()
  df_anova_list_pH[[1]] <- df_anova_pH
  
  for (l in met_list) {
    tmp <- df5 %>% 
      filter(met == l, species %in% c(bact_name, "pH7"))
    
    test <- tmp %>% 
      anova_test(value ~ species*timepoint)
    
    pval <- get_anova_table(test)$p[3]
    
    df_anova_pH <- tmp %>% 
      mutate(pval = pval) %>% 
      filter(species != "pH7")
    
    df_anova_list_pH[[length(df_anova_list_pH)+1]] <- df_anova_pH
  }
  
  df_anova_pval_pH <- bind_rows(df_anova_list_pH)
  
  # Adjust p-values with fdr
  df_anova_padj_pH <- df_anova_pval_pH %>% 
    select(c(Compound, species, met, pval)) %>% 
    distinct()
  
  df_anova_padj_pH$padj <- p.adjust(df_anova_padj_pH$pval, method = "fdr")
  
  # df6 contains metabolites which are increasing or decreasing significantly in species presence
  df6 <- full_join(df_anova_pval_pH, df_anova_padj_pH, 
                   by = c("Compound", "species", "met", "pval"))
  
  ##### 4.3. Get significant metabolites candidates ----
  # df6_filtered also filters for metabolites that are significantly different from pH control
  df6_filtered <- df6 %>% filter(padj <= padj_thresh)
  
  if (nrow(df6_filtered) == 0) {
    next
  }
  
  # Now retrieve pH controls for metabolites from df6_filtered
  metabolites_list3 <- df6_filtered %>%
    arrange(desc(F1_SCORE)) %>% 
    distinct(MZ, RT, COSINE_SCORE, F1_SCORE, met)
  
  df7 <- df5 %>% 
    filter(met %in% metabolites_list3$met)
  
  # Create summary df
  df8 <- df7 %>% 
    mutate(RT = as.numeric(as.character(RT)),
           MZ = as.numeric(as.character(MZ))) %>% 
    group_by(Compound, species, MZ, RT, met, COSINE_SCORE, F1_SCORE, timepoint) %>%
    summarise(mean = mean(value), sd = sd(value)) %>% 
    pivot_wider(names_from = timepoint, values_from = c(mean, sd)) %>% 
    mutate(FC_3vs0 = log2(mean_3/mean_0), 
           FC_9vs0 = log2(mean_9/mean_0),
           FC_27vs0 = log2(mean_27/mean_0)) %>% 
    mutate(sd_0 = ifelse(is.na(sd_0), ms_noise, sd_0),
           sd_3 = ifelse(is.na(sd_3), ms_noise, sd_3),
           sd_9 = ifelse(is.na(sd_9), ms_noise, sd_9),
           sd_27 = ifelse(is.na(sd_27), ms_noise, sd_27)) %>% 
    ungroup()
  
  
  #### 5. Filter data ----
  
  ##### 5.1. Filter out all values where overall intensity < area_thresh ----
  to_exclude <- df8 %>% 
    filter(species == bact_name,
           mean_0 < area_thresh & mean_3 < area_thresh & mean_9 < area_thresh & mean_27 < area_thresh)
  
  df9 <- df8 %>% 
    filter(!met %in% c(to_exclude$met))
  
  unique(df9$met)
  
  ##### 5.2. Get a list of in-source fragments  ----
  # For fragments RT must be matching, therefore use lower RT tolerance (abs_rt/2)
  # Get close RT pairs and filter to exclude pairs with identical MZ values
  rt_pairs <- expand.grid(met1 = df9$met, 
                          met2 = df9$met) %>%
    separate(met1, into = c("MZ1", "RT1"), sep = "_") %>% 
    separate(met2, into = c("MZ2", "RT2"), sep = "_") %>% 
    mutate_all(as.numeric) %>% 
    filter(RT1 <= RT2 & MZ1 != MZ2 & abs(RT1 - RT2) <= abs_rt/2) %>% 
    distinct()
  
  # Now compare patterns in each pair by correlation
  # If correlation coefficient >= 0.7, one of ion is likely a fragments another one
  # Label ions that will be kept (smaller ion is fragment, bigger - no)
  
  fragm <- data.frame()
  fragments_list <- list()
  fragments_list[[1]] <- fragm
  
  if (nrow(rt_pairs) > 0) {
    for (p in 1:nrow(rt_pairs)) {
      RT1 <- rt_pairs$RT1[p]
      RT2 <- rt_pairs$RT2[p]
      
      fragm_df <- df9 %>% 
        filter(species == bact_name) %>% 
        filter(RT == RT1 | RT == RT2) %>% 
        mutate(delta_mz = max(MZ)-min(MZ)) %>% 
        select(-contains("FC"), -COSINE_SCORE, -F1_SCORE) %>% 
        group_by(MZ, RT, met) %>% 
        pivot_longer(cols = matches("^mean_|^sd_"),  # Select columns starting with "mean_" or "sd_"
                     names_to = c(".value", "timepoint"), 
                     names_sep = "_") %>% 
        group_by(Compound, species, MZ, RT, met) %>% 
        mutate(norm_area = mean/max(mean), norm_sd = sd/max(mean))
      
      met1 <- fragm_df %>% filter(RT == RT1) %>% pull(met) %>% unique()
      met2 <- fragm_df %>% filter(RT == RT2) %>% pull(met) %>% unique()
      
      correlation <- cor.test(fragm_df$norm_area[fragm_df$met == met1],
                              fragm_df$norm_area[fragm_df$met == met2])
      
      cor_value <- correlation$estimate[[1]]
      
      # Plot and save pairs where correlation value >= 0.7
      if (cor_value >= 0.7) {
        frag_plot <- fragm_df %>%
          group_by(timepoint) %>% 
          mutate(diff_intensity = max(mean)-mean) %>% 
          ungroup() %>% 
          group_by(met) %>% 
          mutate(max_intensity_met = ifelse(mean(diff_intensity) == 0, met, NA)) %>% 
          ggplot(aes(x = as.numeric(timepoint), y = norm_area, linetype = met)) +
          geom_errorbar(aes(ymin = norm_area - norm_sd, ymax = norm_area + norm_sd), width = 0.1) +
          geom_point() +
          geom_line() +
          labs(x = "time, h", y = "normalized area") +
          scale_linetype_discrete(name = "metabolites") +
          theme_minimal()
        
        frag_plot <- frag_plot +
          annotate("text", x = max(as.numeric(fragm_df$timepoint)), y = max(fragm_df$norm_area), 
                   label = paste("Correlation:", round(cor_value, 4)), hjust = 1)
        
        plot(frag_plot)
        
        frag_plot2 <- fragm_df %>%
          ggplot(aes(x = as.numeric(timepoint), y = mean, linetype = met)) +
          geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
          geom_point() +
          geom_line() +
          labs(x = "time", y = "area") +
          scale_linetype_discrete(name = "metabolites") +
          theme_minimal()
        
        plot(frag_plot2)
        
        fragm <- fragm_df %>% 
          ungroup() %>% 
          select(1,2,3,4,5,6,7,8,9) %>% 
          distinct() %>% 
          mutate(correlation_value = cor_value,
                 fragment = ifelse(MZ == min(MZ), "yes", "no"))
        
        fragments_list[[p+1]] <- fragm
        
        ggsave(filename = paste0(wd_path, "/Metabolites_collection/4_Insource_fragments/", compound_id, "_", bact, "_", p, ".png"),
               plot = frag_plot2,
               width = 6,
               height = 3,
               unit = "in",
               dpi = 300)
      } 
    }
    
    fragments <- rbindlist(fragments_list) 
    if (nrow(fragments)>0){
      write_csv(fragments, paste0(wd_path, "/Metabolites_collection/4_Insource_fragments/", compound_id, "_", bact, ".csv"))
    }
  }
  
  if (!exists("fragments")) {
    fragments <- data.frame(
      Compound = character(),
      species = character(),
      MZ = numeric(),
      RT = numeric(),
      met = character(),
      delta_mz = numeric(),
      correlation_value = numeric(),
      fragment = character())
  }
  
  ##### 5.3. Filter out in-source fragments  ----
  # Get ions labelled as fragments
  if (length(fragments$fragment) > 0) {
    frag_to_exclude <- fragments %>% filter(fragment == "yes")
  } else {
    frag_to_exclude <- data.frame(
      Compound = character(),
      species = character(),
      MZ = numeric(),
      RT = numeric(),
      met = character(),
      delta_mz = numeric(),
      correlation_value = numeric(),
      fragment = character())
  }
  
  # Remove them from df9 and save in df10
  df10 <- df9 %>% filter(!met %in% frag_to_exclude$met) 
  
  write_csv(df10, paste0(wd_path, "/Metabolites_collection/4_Metabolite_data/", compound_id, "_", bact, ".csv"))
  
  metabolites_list_final <- df10 %>%
    filter(species == bact_name) %>% 
    select(-met) %>% 
    arrange(desc(F1_SCORE)) %>% 
    distinct(Compound, species, MZ, RT, COSINE_SCORE, F1_SCORE) %>% 
    mutate(row_number = row_number()) %>% 
    mutate(metabolite = paste0(Compound, "_met", row_number))
  
  if (nrow(metabolites_list_final) > 0) {
    write_csv(metabolites_list_final, paste0(wd_path, "/Metabolites_collection/4_Metabolite_lists/", compound_id, "_", bact, ".csv"))
  }
  
  rm(fragments)
  
  #### 6. Plot metabolites with controls ----
  
  df11 <- df10 %>% 
    select(-contains("FC")) %>% 
    pivot_longer(cols = matches("^mean_|^sd_"),  # Select columns starting with "mean_" or "sd_"
                 names_to = c(".value", "timepoint"), 
                 names_sep = "_") %>% 
    arrange(desc(F1_SCORE))
  
  if (nrow(metabolites_list_final) == 0) {
    message("No hits found: ", paste0(compound, " + ", bact))
    
    next
  }
  
  df11$MZ <- as.factor(df11$MZ)
  df11$RT <- as.factor(df11$RT)
  df11$timepoint <- as.numeric(df11$timepoint)
  
  conditions <- unique(df11$species)
  df11$species <- factor(df11$species, levels = c(bact_name, "pH4", "pH7"))
  
  plot <- df11 %>% 
    unite(score, COSINE_SCORE, F1_SCORE, sep = "/", remove = TRUE) %>%
    group_by(Compound, species, met, timepoint, score) %>%
    ggplot(mapping = aes(x = timepoint, y = mean, ymin = 0,
                         group = interaction(met, species, score))) +
    geom_point(mapping = aes(shape = species), size = 2) +
    geom_line(mapping = aes(linetype = species), color = "black", size = 0.75) +
    geom_errorbar(mapping = aes(ymin = mean - sd, ymax = mean + sd), size = 0.5, width = 0.5) +
    scale_x_continuous(labels = label_number(1)) +
    scale_y_continuous(label = scientific_10) +
    facet_wrap(~ score + met, scales = "free_y", ncol = 1) +
    theme_bw() +
    labs (x= "time, h", y= "peak area") +
    ggtitle(paste0(compound, "_", bact_name)) +
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
  
  ggsave(filename = paste0(wd_path, "/Metabolites_collection/4_Metabolite_plots/", compound_id, "_", bact, ".png"),
         plot = plot,
         width = 3,
         height = 1.5*(length(unique(df11$met))+1),
         unit = "in",
         dpi = 300)
}

