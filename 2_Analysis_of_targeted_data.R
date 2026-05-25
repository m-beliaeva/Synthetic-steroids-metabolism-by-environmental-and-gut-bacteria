rm(list = ls())

#### Import libraries ----
library(tidyverse)
library(scales)
library(tools)
library(rstatix)
library(viridis)
library(ComplexHeatmap)
library(circlize)

#### READ AND PREPARE DATA ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# Read data
data <- read_csv(paste0(wd_path, "/Input/MS1_raw_data.csv")) 

# Check data
any(is.na(data)) #FALSE
unique(data$name)
unique(data$condition)

# Set condition as factor
data$condition <- factor(data$condition, 
                             levels = c("C.testosteroni", "P.aeruginosa",
                                        "R.erythropolis", "S.herbicidovorans",
                                        "B.uniformis", "B.longum", "C.comes",
                                        "C.scindens", "E.faecalis", "E.rectale",
                                        "L.reuteri", "P.vulgatus",
                                        "pH4", "pH7"))

species_env <- c("C.testosteroni", "P.aeruginosa",
                 "R.erythropolis", "S.herbicidovorans")

species_gut <- c("B.uniformis", "B.longum", "C.comes",
                 "C.scindens", "E.faecalis", "E.rectale",
                 "L.reuteri", "P.vulgatus")

# Separate by assay
data_env <- data %>% filter(assay == "env") %>% select(-assay)
data_gut <- data %>% filter(assay == "gut") %>% select(-assay)

#### CALCULATE PVALS ----
data <- data %>% 
  mutate(time = as.factor(time),
         name = as.factor(name),
         position = as.factor(position)) 

class(data$time)
class(data$name)

drug_list <- unique(data$name)

# Test whether each drug is significantly changing in pH7 control
# Calculate p-values based on the time-course only
df_anova_pH7 <- data.frame()
df_anova_list_pH7 <- list()
df_anova_list_pH7[[1]] <- df_anova_pH7

for (eg in unique(data$assay)) {
  for (drug in drug_list) {
    tmp <- data %>% 
      filter(assay == eg, name == drug, condition == "pH7")
    
    test <- tmp %>% 
      anova_test(area ~ time)
    
    pval <- get_anova_table(test) %>% pull(p)
    
    df_anova_pH7 <- tmp %>% 
      mutate(pval = pval)
    
    df_anova_list_pH7[[length(df_anova_list_pH7)+1]] <- df_anova_pH7
  }
}

pH7_data_anova <- bind_rows(df_anova_list_pH7)

# Adjust p-values with fdr
pH7_data_padj <- pH7_data_anova %>% 
  select(-c(time, position, area)) %>% 
  distinct()

pH7_data_padj$padj <- p.adjust(pH7_data_padj$pval, method = "fdr")

pH7 <- full_join(pH7_data_anova, pH7_data_padj, 
                 by = c("name", "condition", "assay", "pval"))

# Two-way interactions between condition and time
condition_list <- unique(data$condition)
condition_list <- condition_list[!grepl("pH7", condition_list)] # pH7 is reference

# Save results here
df_anova <- data.frame()
df_anova_list <- list()
df_anova_list[[1]] <- df_anova

for (eg in unique(data$assay)) {
  for (bact in condition_list) {
    for (drug in drug_list) {
      tmp <- data %>% 
        filter(assay == eg, name == drug, condition %in% c(bact, "pH7"))
      
      if (length(unique(tmp$condition)) > 1) {
        test <- tmp %>% 
          anova_test(area ~ condition*time) # add paired here
        
        pval <- get_anova_table(test)$p[3] #%>% pull(p)
        
        df_anova <- tmp %>% 
          filter(condition != "pH7") %>% 
          mutate(pval = pval)
        
        df_anova_list[[length(df_anova_list)+1]] <- df_anova
      } else {
        message("Skipping ANOVA: missing combination ", drug, "+", bact, 
                " (not ", eg, " bacteria)")
      }
    }
  }
}

data_anova <- bind_rows(df_anova_list)

# Adjust p-values with fdr
data_padj <- data_anova %>% 
  select(-c(time, position, area)) %>% 
  distinct()

data_padj$padj <- p.adjust(data_padj$pval, method = "fdr")

data_stats <- full_join(data_anova, data_padj, 
                        by = c("name", "condition", "assay", "pval"))

# Merge data 
data_stats <- rbind(data_stats, pH7)

#### CHECK COMPOUNDS THAT ARE DEGRADED AT t0 AND IMPUTE FROM CONTROL ----
# Impute mean of pH7 at t0 if the compound is degraded very fast (area(t0)==0)
t0_check1 <- data_stats %>% filter(time == 0 & area == 0) #only some compounds with C.testosteroni

for (i in 1:nrow(t0_check1)) {
  to_fix <- slice(t0_check1, i)
  name_check <- as.character(to_fix$name)
  condition_check <- as.character(to_fix$condition)
  assay_check <- as.character(to_fix$assay)
  
  pH7_t0 <- data_stats %>% filter(name == name_check,
                                  condition == "pH7",
                                  time == 0,
                                  assay == assay_check) %>% 
    group_by(name, condition, time, assay) %>% 
    mutate(mean = mean(area)) %>% 
    select(-area) %>% distinct()
  
  to_impute <- pH7_t0$mean
  
  data_stats <- data_stats %>%
    mutate(area = ifelse(name == name_check 
                         & time == 0 
                         & condition == condition_check 
                         & assay == assay_check, 
                         to_impute, area))
  
}


t0_check2 <- data_stats %>% filter(time == 0 & is.na(area)) #empty
t0_check3 <- data_stats %>% filter(time == 0 & area == 0) #empty

#### CALCULATE DEGRADATION ----
data_mean <- data_stats %>% 
  select(-position) %>% 
  group_by(name, condition, time, assay, pval, padj) %>% 
  summarise(mean = mean(area), sd = sd(area)) 

data_mean_wide <- data_mean %>% 
  pivot_wider(names_from = time, values_from = c(mean, sd))

# Calculate % degradation for each timepoint
data_degr <- data_mean_wide %>% 
  group_by(name, condition, assay) %>% 
  mutate(degr_3 = (100-(mean_3/mean_0*100)), 
         degr_9 = (100-(mean_9/mean_0*100)), 
         degr_27 = (100-(mean_27/mean_0*100))) %>% 
  mutate(degr_3 = if_else(degr_3 < 0, 0, degr_3),
         degr_9 = if_else(degr_9 < 0, 0, degr_9),
         degr_27 = if_else(degr_27 < 0, 0, degr_27)) %>% 
  distinct()

write_csv(data_degr, paste0(wd_path, "/Output/Steroid_degradation_statistics.csv"))

# How many compounds degraded across environmental and gut species
data_degr %>% filter(!condition%in%c("pH4", "pH7")) %>% filter(assay == "env") %>% nrow() #88
data_degr %>% filter(!condition%in%c("pH4", "pH7")) %>% filter(assay == "env") %>% filter(degr_27 >= 50) %>% filter(padj <= 0.05) %>% nrow() #35
data_degr %>% filter(!condition%in%c("pH4", "pH7")) %>% filter(assay == "gut") %>% nrow() #176
data_degr %>% filter(!condition%in%c("pH4", "pH7")) %>% filter(assay == "gut") %>% filter(degr_27 >= 50) %>% filter(padj <= 0.05) %>% nrow() #24


compounds <- data_degr %>% 
  pull(name) %>% 
  unique() #22

degraded_compounds <- data_degr %>% 
  filter(!condition%in%c("pH4", "pH7")) %>% 
  filter(degr_27 >= 50) %>% filter(padj <= 0.05) %>% 
  pull(name) %>% 
  unique() #17

non_degraded_compounds <- setdiff(compounds, degraded_compounds) 

#### PLOT HEATMAPS ----
data_degr$name <- gsub("12-Dihydro-dexamethasone", "1,2-Dihydro-dexamethasone", data_degr$name)
data_degr$name <- gsub("67-Dehydro-ethinylestradiol", "6,7-Dehydro-ethinylestradiol", data_degr$name)
data_degr$name <- gsub("Dexamethasone-sodium-phosphate", "Dexamethasone-phosphate", data_degr$name)


##### Add order -----
compound_order <- c("Finasteride", "Altrenogest", "Cyproterone-acetate", "Chlormadinone-acetate", "Betamethasone-dipropionate", 
                    "Beclomethasone-dipropionate", "Fludrocortisone-acetate", "Triamcinolone-acetonide", 
                    "Dexamethasone-phosphate", "1,2-Dihydro-dexamethasone", "Flumethasone", "Triamcinolone",
                    "Dexamethasone", "Betamethasone", "Exemestane", "Estrone", "6,7-Dehydro-ethinylestradiol",
                    "17-Alpha-ethinylestradiol", "Inostral", "Testosterone", "Levonorgestrel", "Tibolone")

species_order <- c("S.herbicidovorans", "C.testosteroni", "P.aeruginosa",  "R.erythropolis", 
                   "control NB (pH4)", "control NB (pH7)", 
                   "B.longum", "B.uniformis", "P.vulgatus", "E.faecalis", "L.reuteri" , "C.scindens", "E.rectale", "C.comes", 
                   "control GMM (pH4)", "control GMM (pH7)")

##### Long data format ----
data_degr_long <- data_degr %>% 
  pivot_longer(cols = starts_with("mean_") | starts_with("sd_") | starts_with("degr_"),
               names_to = c(".value", "time"),
               names_pattern = "(mean|sd|degr)_(.*)") %>% 
  mutate(condition = case_when(
    condition == "pH4" & assay == "env" ~ "control NB (pH4)",
    condition == "pH7" & assay == "env" ~ "control NB (pH7)",
    condition == "pH4" & assay == "gut" ~ "control GMM (pH4)",
    condition == "pH7" & assay == "gut" ~ "control GMM (pH7)",
    TRUE ~ condition)) %>% 
  mutate(time = as.numeric(time),
         name = factor(name, levels = compound_order),
         condition = factor(condition, levels = species_order)) %>% 
  replace_na(list(degr = 0))

# Split into env vs gut
row_groups <- factor(ifelse(species_order %in% c("S.herbicidovorans", "C.testosteroni", "P.aeruginosa", 
                                                 "R.erythropolis", "control NB (pH4)", "control NB (pH7)"), 
                            "Env. species", "Gut species"))

#### Heatmaps ----
timepoints <- c(3,9,27)

for (t in timepoints) {
  heatmap <- data_degr_long %>% 
    filter(time == t) %>% 
    select(name, condition, degr) %>%
    pivot_wider(names_from = name, values_from = degr) %>%
    column_to_rownames(var = "condition")
  
  heatmap <- heatmap[species_order, , drop = FALSE]
  heatmap <- heatmap[, compound_order]
  
  
  pdf_file <- paste0(wd_path, "/Plots/heatmap_", t, "h.pdf")
  pdf(pdf_file, width = 6, height = 6)
  
  ht <-  Heatmap(
    as.matrix(t(heatmap)),  
    name = paste0("Degradation (%) at t = ", t, "h"),
    column_order = species_order,  
    column_split = row_groups,   
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),  
    col = viridis(256),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = gpar(fontface = "italic"),  
    row_names_side = "left",
    border = TRUE,
    column_gap = unit(5, "mm"),  
    rect_gp = gpar(col = "white", lwd = 1.5),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(2, "cm"),  
      legend_width = unit(3, "cm"),
      direction = "horizontal",
      title_position = "lefttop"))
  
  draw(ht, heatmap_legend_side = "top")
  
  dev.off()
  
}
