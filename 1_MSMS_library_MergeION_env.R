rm(list = ls())

## THIS CODE CALCULATES SIMILARITY-TO-PARENT-COMPOUND LISTS FROM DDA MSMS ##

#### INSTALL PACKAGES AND LIBRIARIES ----
##### Packages and dependencies ----
# # Install BiocManager if it has not been installed previously:
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# 
# # Check Java
# Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/jdk-20.jdk/Contents/Home")
# install.packages("rJava")
# 
# # Install dependencies 
# BiocManager::install("RMassBank")
# 
# # Install MergeION
# install.packages("remotes")
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") 
# remotes::install_github("daqana/dqshiny", force = T)
# remotes::install_github("schymane/RChemMass", force = T)
# BiocManager::install("daniellyz/MergeION2")

##### Load libraries ----
library(rJava)
library(MergeION) 
library(RChemMass)
library(tidyverse)

#### PATHS AND METADATA ----
# Set path for metacode
# This path contains subfolders with compounds' names
wd_path <- getwd()
setwd(wd_path)

# Set path to include lists used for MSMS acquisition
metadata_path <- paste0(wd_path, "/MSMS_include_lists/")
# Set path to MSMS data folder
msms_path <- "/Volumes/Mariia/MAB001_Steroid_metabolism/ms_data_001C/Environmental_bacteria/converted_data/"

# Import parent compounds
parents <- read_csv("Reference_compounds_metadata.csv")

# Get list of compound-species combinations 
msms_files <- list.files(
  path = msms_path,
  pattern = "\\.mzML$", 
  full.names = FALSE) %>% 
  str_remove("\\.mzML$") %>%              
  str_remove("_\\d{2}h") %>%
  str_remove("_10eV") %>% 
  str_remove("MAB001_") %>% 
  unique()

#### METACODE ####
for (i in msms_files) {
  index <- tail(strsplit(i, "_")[[1]], n=1)
  
  if (index == "Ct") {
    bact <- "Ctestosteroni"
  } else if (index == "Pa") {
    bact <- "Paeruginosa"
  } else if (index == "Re") {
    bact <- "Rerythropolis"
  } else if (index == "Sh") {
    bact <- "Sherbicidovorans"
  }
  
  # Retrieve compound information 
  id <- head(strsplit(i, "_")[[1]], n=1)
  
  # Retrieve MSMS metadata (include lists) 
  if (!file.exists(paste0(metadata_path, id, "_", bact, ".csv"))) {
    message("File does not exist:", paste0(metadata_path, id, "_", bact, ".csv"))
    
    next
  }
  
  include_list <- read_csv(paste0(metadata_path, id, "_", bact, ".csv"))
  metadata_file <-  include_list %>% 
    rename(PEPMASS = "Prec. m/z",	RT = "Ret. Time (min)") %>% 
    mutate(ID = paste0(id, "_", row_number())) %>% 
    select(c(ID, PEPMASS, RT))
  
  # Get list of MSMS files
  lcms_files <- list.files(path = msms_path, pattern = ".mzML") %>% 
    str_subset(index) %>%
    str_subset(id)
  
  ##### Create local MSMS library with data ----
  # It is important to do from the location of mzML files!
  setwd(msms_path)
  
  # Set parameters for library creation
  # Here you can change mz and RT tolerance and noise (baseline) level
  input_library = NULL 
  processing.algorithm = "Default"
  params.search = list(mz_search = 0.002, ppm_search = 10, rt_search = 6, rt_gap = 10)
  params.ms.preprocessing = list(normalized = TRUE, baseline = 200, relative = 0, max_peaks = 200, recalibration = 0)
  params.consensus = list(consensus = TRUE, consensus_method = "consensus", consensus_window = 0.005)
  library1 = library_generator(input_library = input_library, 
                               lcms_files = lcms_files, 
                               metadata_file = metadata_file,
                               polarity = "Positive", mslevel = 2, add.adduct = FALSE, 
                               processing.algorithm = processing.algorithm, 
                               params.search = params.search,
                               params.ms.preprocessing = params.ms.preprocessing, 
                               params.consensus = params.consensus)
  
  # Export library
  library_writer(library1, con = paste0(wd_path, "/MergeION/mgf/", id, "_", bact, ".mgf"), type ="complete")
  library_writer(library1, con = paste0(wd_path, "/MergeION/Rdata/", id, "_", bact, ".Rdata"), type ="complete")
  
  #### Do spectral comparison ----
  # Get parent data and related metabolites for the reference MSMS
  id_core <- head(strsplit(id, "-")[[1]], n=1)
  
  if (nchar(id) == nchar(id_core)) {
    parents_id <- parents$PARENT_ID[parents$PARENT_ID == id]
  } else {
    parents_id <- parents$PARENT_ID[grepl(id_core, parents$PARENT_ID)]
  }
  
  # Loop over reference compounds
  for (compound in parents_id) {
    parent_mz <- parents$MZ[parents$PARENT_ID == compound]
    parent_rt <- parents$RT[parents$PARENT_ID == compound]
    
    # Extracting spectral records of compound of interest based on its precursor mass and retention time 
    query_result = library_query(input_library = library1, query_expression = c(paste0("PEPMASS", "=", parent_mz), paste0("RT", "=", parent_rt)))
    library_reporter(query_result) 
    print(query_result$complete$metadata) 
    
    # Search parent analogs in created library (no external library used) 
    query_sp = query_result$consensus$sp[[1]] # Query_sp is a two-column matrix m/z, intensity of the parent compound consensus spectrum
    query_mz = parent_mz # Precursor mass 
    
    # Using F1 algorithm
    params.query.sp = list(prec_mz = query_mz, use_prec = FALSE, polarity = "Positive", method = "F1", min_frag_match = 3, min_score = 0, reaction_type = "Metabolic")
    search_F1 = library_query(input_library = library1, 
                              query_spectrum = query_sp, 
                              params.query.sp = params.query.sp)
    
    head(search_F1$consensus$metadata)
    write_csv(search_F1$consensus$metadata, paste0(wd_path, "/MergeION/Output/", compound, "_", bact, "_parent_similarity_F1.csv"))
    
    # Using Cosine algorithm
    params.query.sp = list(prec_mz = query_mz, use_prec = FALSE, polarity = "Positive", method = "Cosine", min_frag_match = 3, min_score = 0, reaction_type = "Metabolic")
    search_cosine = library_query(input_library = library1, 
                                  query_spectrum = query_sp, 
                                  params.query.sp = params.query.sp)
    
    head(search_cosine$consensus$metadata)
    write_csv(search_cosine$consensus$metadata, paste0(wd_path, "/MergeION/Output/", compound, "_", bact, "_parent_similarity_cosine.csv"))
  }
  setwd(wd_path)
}
