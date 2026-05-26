rm(list = ls())

#### Import libraries ----
library(tidyverse)
library(scales)
library(tools)
library(rstatix)
library(pheatmap)
library(rJava)
library(MergeION) 
library(RChemMass)
library(MSnbase)
library(ChemmineR)
library(rcdk)
library(data.table) 

#### SET THRESHOLDS ----
# Tolerance for MS
abs_mass = 0.002 # absolute mz tolerance
rel_mass = 2*10^-6 # relative mz tolerance (20 ppm)
abs_rt = 0.15*60 # absolute RT tolerance in seconds 
#rel_rt = 0.02 # relative RT tolerance

#### READ DATA ----

# Set paths to working directory and mgf folder
wd_path <- getwd()
mgf_path <- paste0(wd_path, "/MergeION/mgf/")
setwd(wd_path)

# Read data, add manually detected metabolites 
metabolites_manual <- read_csv(paste0(wd_path, "/Manually_detected_metabolites.csv")) %>% 
  mutate(delta_mass = round(MZ-MZ_PARENT,4),
         identified_by = "manual_curation") %>% 
  distinct()

metabolites <- read_csv(paste0(wd_path, "/Metabolites_collection/5_Metabolites_filtered.csv")) %>%
  mutate(identified_by = "MSMS") %>% #168 rows
  bind_rows(metabolites_manual) #210 rows

#### CREATE FILTERED MFG ---- 
files <- list.files(mgf_path, pattern = "\\.mgf$", full.names = TRUE)

unique_interactions <- metabolites %>% select(species, Compound) %>% distinct() #70

for (n in 1:nrow(unique_interactions)) {
  tmp <- unique_interactions %>% slice(n)
  species <- tmp$species %>% gsub("\\.", "", .)
  compound <- tmp$Compound
  met_list <- metabolites %>%  filter(Compound == compound, species == species)
  MZ_met <- met_list %>% pull(MZ)
  RT_met <- met_list %>% pull(RT)*60 # In seconds

  # Read mgf file of certain compound-species combination
  combination <- paste0(compound, "_", species)
  combination <- str_replace(combination, "Dexamethasone-phosphate", "Dexamethasone-sodium-phosphate")
  filename <- files[grepl(combination, files)]
  
  if (length(filename)==0) {
    next
  }
  
  # Prepare mgf to the correct RT format
  filename_fixed <- "tmp.mgf"
  mgf_lines <- readLines(filename)
  mgf_lines <- gsub("^RT=([0-9.]+)", "RTINSECONDS=\\1", mgf_lines) 
  mgf_lines <- sapply(mgf_lines, function(line) {
    if (grepl("^RTINSECONDS=", line)) {
      val <- as.numeric(sub("RTINSECONDS=", "", line))
      line <- paste0("RTINSECONDS=", val * 60)
    }
    line
  }, USE.NAMES = FALSE)
  
  writeLines(mgf_lines, filename_fixed)
  mgf_raw <- readMgfData(filename_fixed)

  # Extract m/z and RT data from the mgf file
  mz_values <- precursorMz(mgf_raw)
  rt_values <- rtime(mgf_raw)
  
  filtered_indices <- integer(0)
  filtered_mgf_data <- new("MSnExp")

  # Loop through each metabolite in met_list
  for (i in seq_along(MZ_met)) {
    mz_range <- c(MZ_met[i] - abs_mass, MZ_met[i] + abs_mass)  # Define m/z range
    rt_range <- c(RT_met[i] - abs_rt, RT_met[i] + abs_rt)  # Define RT range
    
    # Find indices that match the criteria
    matched_indices <- which(mz_values >= mz_range[1] & mz_values <= mz_range[2] &
                               rt_values >= rt_range[1] & rt_values <= rt_range[2])
    
    # Append the matched indices to the filtered_indices list
    filtered_indices <- unique(c(filtered_indices, matched_indices))
  }
  
  filtered_mgf_data <- mgf_raw[filtered_indices]
  new_filename <- basename(filename)
  writeMgfData(filtered_mgf_data, paste0(paste0(wd_path, "/Metabolites_collection/6_Filtered_mgf/", new_filename)))
}

#### SUGGEST ATOMIC DELTA FORMULA ----
# Generate parent formula first 
df_formula <- metabolites %>%
  mutate(FORMULA_PARENT = sapply(SMILES_PARENT, MolFormFromSmiles.rcdk)) %>% 
  distinct()

# Find difference between parent and metabolite in atoms
metabolites_list <- list()
df <- data.frame()

for (m in 1:nrow(df_formula)) {
  d_mass <- df_formula$delta_mass[m]
  parent_frml <- df_formula$FORMULA_PARENT[m]
  
  if (d_mass <= 0) {
    elements_list <- list(c("C", 0, 30), c("H", 0, 50), c("O", 0, 10))
  } else {
    elements_list <- list(c("C", 0, 30), c("H", 0, 50), c("P", 0, 2),
                          c("N", 0, 5), c("O", 0, 10), c("S", 0, 5))
  }  
  
  if (grepl("P", parent_frml) & grepl("F", parent_frml)) {
    elements_list <- append(elements_list, list(c("P", 0, 1))) %>% 
      append(list(c("F", 0, 2)))
  } else if (grepl("F", parent_frml)) {
    elements_list <- append(elements_list, list(c("F", 0, 2)))
  } else if (grepl("Cl", parent_frml)) {
    elements_list <- append(elements_list, list(c("Cl", 0, 2)))
  } else if (grepl("P", parent_frml)) {
    elements_list <- append(elements_list, list(c("P", 0, 1)))
  } else {
    elements_list <- elements_list
  }
  
  tmp <- generate.formula(abs(d_mass),
                          window = abs_mass*2,  # allow wider mass tolerance
                          elements = elements_list,
                          charge = 0)
  if (length(tmp) > 0){
    tmp_string <- c()
    for (k in 1:length(tmp)) {
      tmp_string[k] <- tmp[[k]]@string 
    }
    
    tmp_string_comb <- paste(unlist(tmp_string), collapse = ", ")
    
    df <- df_formula %>% 
      filter(FORMULA_PARENT == parent_frml & delta_mass == d_mass) %>% 
      mutate(atomic_difference = c(tmp_string_comb))
  } else {
    df <- df_formula %>% 
      filter(FORMULA_PARENT == parent_frml & delta_mass == d_mass) %>% 
      mutate(atomic_difference = NA)
  }
  
  metabolites_list[[m]] <- df
}

df_atomic_diff <- do.call(rbind, metabolites_list) %>% distinct()

#### ADD REACTIONS METADATA ----
reactions <- read_csv(paste0(wd_path, "/Reaction_mass_differences.csv"))
reaction_list <- list()
met_react <- data.frame()

for (j in 1:nrow(reactions)) {
  react <- reactions$reaction[j]
  mass_diff <- reactions$mass_difference[j]
  
  met_react <- df_atomic_diff %>% 
    mutate(dev = abs(delta_mass-mass_diff)) %>% 
    filter(dev <= abs_mass*2 | dev <= abs(rel_mass*mass_diff)) %>% # allow wider mass tolerance
    mutate(reaction = react) %>% 
    select(-dev)
  
  reaction_list[[j+1]] <- met_react
}

met_reactions <- rbindlist(reaction_list, fill = TRUE)

df_reactions <- left_join(df_atomic_diff, met_reactions, by=c("species", "Compound", "MZ", "RT", "MZ_PARENT", "RT_PARENT", "SMILES_PARENT", 
                                                              "delta_mass", "parent_match", "identified_by", "FORMULA_PARENT", "atomic_difference"))

write_csv(df_reactions, paste0(wd_path, "/Metabolites_collection/6_Metabolites_filtered_with_reactions.csv")) 