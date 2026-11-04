rm(list = ls())

#### Import libraries ----
library(tidyverse)
library(ggembl)
library(data.table)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggalluvial)
library(stringr)

#### SET WD ----

# Set paths to working directory and mgf folder
wd_path <- getwd()
setwd(wd_path)

#### METABOLITE DATA -----
# Get list of curated metabolites
metabolites <- read_csv(paste0(wd_path, "/Metabolites_collection/6_Metabolites_filtered_with_reactions_curated_upd.csv")) %>% 
  mutate(parent_match = NA) # 166 curated biotransformations detected excluding parent hits

metabolites$Compound <- gsub("^1,?2-Dihydro-?dexamethasone$",
                             "1,2-Dihydro-dexamethasone", metabolites$Compound)
metabolites$Compound <- gsub("^6,?7-Dehydro-?eth[iy]nylestradiol$",
                             "6,7-Dehydro-ethinylestradiol", metabolites$Compound)

# Tolerance for MZ_RT consensus
abs_mass = 0.001 # absolute mz tolerance
rel_mass = 20*10^-6 # relative mz tolerance (20 ppm)
# allows for between-runs shift, however, merges some isomers (to be fixed manually)
abs_rt = 0.12 # absolute rt tolerance

# Cluster metabolites based on MZ_RT combinations by tolerance
metabolites_clustered <- metabolites %>%
  arrange(MZ, RT) %>%
  mutate(cluster_id = NA_integer_)

cluster <- 1
metabolites_clustered$cluster_id[1] <- cluster

# Reference values for current cluster
cluster_mz <- metabolites_clustered$MZ[1]
cluster_rt <- metabolites_clustered$RT[1]

for (i in 2:nrow(metabolites_clustered)) {
  
  mz_curr <- metabolites_clustered$MZ[i]
  rt_curr <- metabolites_clustered$RT[i]
  
  mz_tol <- max(abs_mass, rel_mass * cluster_mz)
  
  if (abs(mz_curr - cluster_mz) <= mz_tol &&
      abs(rt_curr - cluster_rt) <= abs_rt) {
    
    metabolites_clustered$cluster_id[i] <- cluster
    
  } else {
    
    cluster <- cluster + 1
    metabolites_clustered$cluster_id[i] <- cluster
    
    # Start new cluster
    cluster_mz <- mz_curr
    cluster_rt <- rt_curr
  }
}

# --- manual correction of clusters that merge distinct compounds ----------
# Three m/z windows chain across the RT tolerance and merge compounds that
# are chromatographically distinct: the 17-ketone epimers (m/z 333.19,
# dexamethasone- vs betamethasone-17-ketone) and the monopropionate
# positional isomers (m/z 449.23 and m/z 465.20). Rows in these windows are
# re-assigned to a sub-cluster by the curated compound name rather than by
# retention time, because between-run RT drift is larger than the separation
# between the two isomers.

manual_splits <- tibble::tribble(
  ~mz_low,  ~mz_high, ~label,
  333.180,  333.190,  "m333_17ketone",
  393.200,  393.210,  "m393_parent",
  449.225,  449.240,  "m449_monopropionate",
  465.195,  465.210,  "m465_monopropionate")

metabolites_clustered <- metabolites_clustered %>%
  mutate(cluster_id = as.character(cluster_id))

for (i in seq_len(nrow(manual_splits))) {
  s   <- manual_splits[i, ]
  hit <- metabolites_clustered$MZ >= s$mz_low &
    metabolites_clustered$MZ <= s$mz_high
  # every row in these windows must carry a curated name
  stopifnot(all(!is.na(metabolites_clustered$METABOLITE_PROPOSED[hit]) &
                  nzchar(trimws(metabolites_clustered$METABOLITE_PROPOSED[hit]))))
  metabolites_clustered$cluster_id[hit] <- paste0(
    s$label, "_", trimws(metabolites_clustered$METABOLITE_PROPOSED[hit]))
}

# how many rows each manual window captured, as a check
metabolites_clustered %>%
  filter(str_starts(cluster_id, "m333|m393|m449|m465")) %>%
  count(cluster_id)

cluster_consensus <- metabolites_clustered %>%
  group_by(cluster_id) %>%
  summarise(
    consensus_MZ = mean(MZ, na.rm = TRUE),
    consensus_RT = mean(RT, na.rm = TRUE),
    .groups = "drop")

metabolites <- metabolites_clustered %>%
  left_join(cluster_consensus, by = "cluster_id") %>%
  mutate(
    consensus_MZ_RT = paste0(round(consensus_MZ, 4), "_", round(consensus_RT, 1))
  ) %>%
  write_csv(paste0(wd_path,
                   "/Metabolites_collection/6_Metabolites_filtered_with_reactions_curated_upd_clustered.csv")) %>%
  select(-consensus_MZ, -consensus_RT)   

stopifnot(n_distinct(metabolites$cluster_id) ==
            n_distinct(metabolites$consensus_MZ_RT))

##### Count metabolites ----
# Funtion to count
count_by_id <- function(df, desc, filter_expr, count_expr) {
  
  count_expr <- enquo(count_expr)
  filter_expr <- enquo(filter_expr)
  
  tibble(
    description = desc,
    MSMS = df %>%
      filter(!!filter_expr, identified_by == "MSMS") %>%
      summarise(n = n_distinct(!!count_expr)) %>%
      pull(n),
    
    manual_curation = df %>%
      filter(!!filter_expr, identified_by == "manual_curation") %>%
      summarise(n = n_distinct(!!count_expr)) %>%
      pull(n),
    
    total = df %>%
      filter(!!filter_expr) %>%
      summarise(n = n_distinct(!!count_expr)) %>%
      pull(n))}

# Summary table 
summary_table <- bind_rows(
  
  # Total biotransformation number
  count_by_id(
    metabolites,
    "Total number of biotransfomrations",
    is.na(parent_match),
    paste(species, Compound, consensus_MZ_RT)),
  
  # Biotransformations with assigned reaction
  count_by_id(
    metabolites,
    "Biotransformations with assigned reaction",
    is.na(parent_match) & !is.na(reaction),
    paste(species, Compound, consensus_MZ_RT)),
  
  # Biotransformations with initial metabolites
  count_by_id(
    metabolites,
    "Biotransformations with initial metabolites",
    is.na(parent_match) & str_detect(metabolite_type, "initial"),
    paste(species, Compound, consensus_MZ_RT)),
  
  # Biotransformations with intermediate metabolites
  count_by_id(
    metabolites,
    "Biotransformations with intermediate metabolites",
    is.na(parent_match) & str_detect(metabolite_type, "intermediate"),
    paste(species, Compound, consensus_MZ_RT)),
  
  # Biotransformations with accumulated metabolites
  count_by_id(
    metabolites,
    "Biotransformations with accumulated metabolites",
    is.na(parent_match) & str_detect(metabolite_type, "accumulated"),
    paste(species, Compound, consensus_MZ_RT)),
  
  # Unique metabolites related to compound (Unique Compound–MZ–RT)
  count_by_id(
    metabolites,
    "Unique Compound-MZ-RT combinations",
    is.na(parent_match),
    paste(Compound, consensus_MZ_RT)),
  
  # Unique metabolites (MZ–RT)
  count_by_id(
    metabolites,
    "Unique MZ-RT combinations",
    is.na(parent_match),
    paste(consensus_MZ_RT)),
  
  # Unique metabolites (MZ–RT)
  count_by_id(
    metabolites,
    "Unique MZ-RT with COSMIC_STRUCTURE_SCORE >= 80",
    is.na(parent_match) & COSMIC_STRUCTURE_SCORE >= 80,
    paste(consensus_MZ_RT)),
  
  # Associated parent compounds
  count_by_id(
    metabolites,
    "Associated parent compounds",
    is.na(parent_match),
    Compound),
  
  # Associated species
  count_by_id(
    metabolites,
    "Associated species",
    is.na(parent_match),
    species),
  
  # Unique compound–species combinations
  count_by_id(
    metabolites,
    "Unique compound-species combinations",
    is.na(parent_match),
    paste(Compound, species)),
  
  # Compound-species combinations with assigned reactions
  count_by_id(
    metabolites,
    "Unique compound-species combinations with assigned reactions",
    is.na(parent_match) & !is.na(reaction),
    paste(Compound, species)
  )
)

# Add annotations count
level_order <- c("1", "2b", "3", "5")

metabolites_best <- metabolites %>%
  filter(is.na(parent_match)) %>%
  mutate(level_rank = match(ANNOTATION_LEVEL, level_order)) %>%
  group_by(cluster_id) %>%
  mutate(best_level = level_order[min(level_rank, na.rm = TRUE)]) %>%
  ungroup()

stopifnot(n_distinct(metabolites$cluster_id) ==
            n_distinct(metabolites$consensus_MZ_RT))

# per biotransformation (species x compound x feature): no ambiguity
level_rows <- map(level_order, function(lv) {
  count_by_id(
    metabolites,
    paste0("Biotransformations at annotation level ", lv),
    is.na(parent_match) & ANNOTATION_LEVEL == lv,
    paste(species, Compound, consensus_MZ_RT))
})

# per unique metabolite (MZ-RT): assign each feature its best level first
metabolites_best <- metabolites %>%
  filter(is.na(parent_match)) %>%
  mutate(level_rank = match(ANNOTATION_LEVEL, level_order)) %>%
  group_by(consensus_MZ_RT) %>%
  mutate(best_level = level_order[min(level_rank, na.rm = TRUE)]) %>%
  ungroup()

level_rows_unique <- map(level_order, function(lv) {
  count_by_id(
    metabolites_best,
    paste0("Unique MZ-RT at annotation level ", lv),
    best_level == lv,
    paste(consensus_MZ_RT))
})

summary_table <- bind_rows(
  summary_table,
  bind_rows(level_rows),
  bind_rows(level_rows_unique))

write_csv(summary_table, paste0(wd_path, "/Metabolites_collection/7_Summary_of_unique_metabolites_upd.csv"))

#### BARPLOT METABOLITES PER STEROID ----
df_plot <- metabolites %>%
  filter(is.na(parent_match)) %>% 
  mutate(MetaboliteID = consensus_MZ_RT)

# Count unique metabolites per Compound and reaction status
metabolite_counts <- df_plot %>%
  group_by(Compound, reaction_status = ifelse(is.na(reaction), "No", "Yes")) %>%
  summarise(unique_metabolites = n_distinct(MetaboliteID), .groups = "drop") %>% 
  write_csv(paste0(wd_path, "/Metabolites_collection/7_Unique_metabolites_per_steroid.csv"))

# Add order
compound_order <- metabolite_counts %>%
  group_by(Compound) %>%
  summarise(total_metabolites = sum(unique_metabolites)) %>%
  arrange(desc(total_metabolites)) %>%
  pull(Compound)

metabolite_counts$Compound <- factor(metabolite_counts$Compound, levels = compound_order)

# Make stacked barplot
plot <- metabolite_counts %>%
  ggplot(aes(x = Compound, y = unique_metabolites, fill = reaction_status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  theme_minimal(base_size = 14) +
  labs(x = "Compound",
       y = "# of unique metabolites",
       fill = "Reaction assigned") +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

print(plot)

ggsave(filename = paste0(wd_path, "/Metabolites_collection/7_Count_by_steroid.pdf"),
       plot = plot,
       width = length(unique(metabolite_counts$Compound))/2+2,
       height = 6,
       unit = "in",
       dpi = 300)

#### BARPLOT METABOLITES PER SPECIES ----
df_plot <- metabolites %>%
  filter(is.na(parent_match)) %>% 
  mutate(MetaboliteID = consensus_MZ_RT)

# Count unique metabolites per species and reaction status
metabolite_counts <- df_plot %>%
  group_by(species, reaction_status = ifelse(is.na(reaction), "No", "Yes")) %>%
  summarise(unique_metabolites = n_distinct(MetaboliteID), .groups = "drop") %>% 
  write_csv(paste0(wd_path, "/Metabolites_collection/7_Unique_metabolites_per_species.csv"))

# Add order
species_order <- metabolite_counts %>%
  group_by(species) %>%
  summarise(total_metabolites = sum(unique_metabolites)) %>%
  arrange(desc(total_metabolites)) %>%
  pull(species)

metabolite_counts$species <- factor(metabolite_counts$species, levels = species_order)

# Make stacked barplot
plot <- metabolite_counts %>%
  ggplot(aes(x = species, y = unique_metabolites, fill = reaction_status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  theme_minimal(base_size = 14) +
  labs(x = "Species",
       y = "# of unique metabolites",
       fill = "Reaction assigned") +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

print(plot)

ggsave(filename = paste0(wd_path, "/Metabolites_collection/7_Count_by_species.pdf"),
       plot = plot,
       width = length(unique(metabolite_counts$species))/2+2,
       height = 6,
       unit = "in",
       dpi = 300)

#### PLOT METABOLITE MAP ----

# Get parent data
data_degr <- read_csv(paste0(dirname(wd_path), "/1_Targeted_analysis/Output/Steroid_degradation_statistics.csv"))
data_degr$name <- gsub("12-Dihydro-dexamethasone", "1,2-Dihydro-dexamethasone", data_degr$name)
data_degr$name <- gsub("67-Dehydro-ethinylestradiol", "6,7-Dehydro-ethinylestradiol", data_degr$name)
data_degr$name <- gsub('Dexamethasone-sodium-phosphate', 'Dexamethasone-phosphate', data_degr$name)

# Add order of species and compounds
compound_order <- c("Finasteride", "Altrenogest", "Cyproterone-acetate", "Chlormadinone-acetate", "Betamethasone-dipropionate", 
                    "Beclomethasone-dipropionate", "Fludrocortisone-acetate", "Triamcinolone-acetonide", 
                    "Dexamethasone-phosphate", "1,2-Dihydro-dexamethasone", "Flumethasone", "Triamcinolone",
                    "Dexamethasone", "Betamethasone", "Exemestane", "Estrone", "6,7-Dehydro-ethinylestradiol",
                    "17-Alpha-ethinylestradiol", "Inostral", "Testosterone", "Levonorgestrel", "Tibolone")

species_order <- c("S.herbicidovorans", "C.testosteroni", "P.aeruginosa",  "R.erythropolis", 
                   "control NB (pH4)", "control NB (pH7)", 
                   "B.longum", "B.uniformis", "P.vulgatus", "E.faecalis", "L.reuteri" , "C.scindens", "E.rectale", "C.comes", 
                   "control GMM (pH4)", "control GMM (pH7)")

# Long data format
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

# Prepare and combine metabolite and parent data
metabolites_df <- metabolites %>% 
  separate(reaction, into = c("reaction1", "reaction2", "reaction3"), sep = " & ") %>% 
  pivot_longer(cols = c("reaction1", "reaction2", "reaction3"), values_to = "reaction") %>% 
  select(-c(name, SMILES_PARENT, FORMULA_PARENT, MZ_PARENT, RT_PARENT, MZ, RT, delta_mass, atomic_difference, 
            SIRIUS_FORMULA, SIRIUS_FORMULA_SCORE, CSIFingerID_STRUCTURE, COSMIC_STRUCTURE_SCORE, METABOLITE_PROPOSED,
            identified_by, Metabolite, ANNOTATION_LEVEL, LEVEL_EVIDENCE, metabolite_type, SIRIUS_PARAMETERS, 
            parent_match, consensus_MZ_RT, cluster_id)) %>% 
  distinct() %>% 
  drop_na(reaction)

df <- data_degr_long %>%
  filter(time == 27 & !grepl("pH", condition, ignore.case = TRUE)) %>%
  rename(species = 'condition', Compound = 'name') %>% 
  select(Compound, species, degr) %>%
  full_join(metabolites_df) %>%
  mutate(#reaction = factor(reaction, levels = reaction_order),
    Compound = factor(Compound, levels = compound_order),
    species = factor(species, levels = species_order))

any(is.na(df$degr))
any(is.na(df$Compound))

df_full <- df %>%
  mutate(interaction1 = as.factor(paste(Compound, species, sep = " + ")),
         interaction2 = as.factor(paste(species, Compound, sep = " + ")))

# Ordered by species
heatmap <- df_full %>% 
  arrange(species) %>% 
  select(interaction2, reaction) %>%
  drop_na() %>% 
  distinct() %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = interaction2, values_from = value, values_fill = 0) %>% 
  column_to_rownames(var = "reaction")

legend_labels <- c("metabolite detected", "NA")

heatmap_anno <- df_full %>% 
  drop_na(reaction) %>% 
  arrange(species) %>% 
  select(interaction2, degr) %>% 
  distinct()

degradation <- HeatmapAnnotation(`Degradation (%) at t = 27h` = 
                                   anno_barplot(heatmap_anno$degr, gp = gpar(fill = "grey30")), 
                                 annotation_name_gp = gpar(fontsize = 14, fontface = "bold"))

column_groups <- gsub("\\s*\\+.*", "", heatmap_anno$interaction2)
column_groups <- factor(column_groups, levels = species_order)

pdf_file <- paste0(wd_path, "/Metabolites_collection/7_Metabolite_map_by_species.pdf")
pdf(pdf_file, width = 12, height = 7.7)

ht <- Heatmap(as.matrix(heatmap),
              col = c("grey90","#440154FF"), 
              top_annotation = degradation,
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              show_row_names = TRUE, 
              show_column_names = TRUE, 
              column_split = column_groups,
              column_title = NULL,
              border = TRUE,
              row_gap = unit(5, "mm"),
              rect_gp = gpar(col = "white", lwd = 1.5),
              #column_names_rot = 45,
              row_names_gp = gpar(fontsize = 10, hjust = 1, vjust = 1),
              column_names_gp = gpar(fontsize = 10, hjust = 1, vjust = 1),
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 14, fontface = "bold"),
                legend_height = unit(4, "cm"),
                legend_width = unit(4, "cm"),
                title = "",
                labels = legend_labels))

draw(ht, heatmap_legend_side = "top", padding = unit(c(60, 20, 5, 10), "mm"))

dev.off()

# Ordered by compound
heatmap <- df_full %>% 
  arrange(Compound) %>% 
  select(interaction1, reaction) %>%
  drop_na() %>% 
  distinct() %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = interaction1, values_from = value, values_fill = 0) %>% 
  column_to_rownames(var = "reaction")

heatmap_anno <- df_full %>% 
  drop_na(reaction) %>% 
  arrange(Compound) %>% 
  select(interaction1, degr) %>% 
  distinct()

degradation <- HeatmapAnnotation(`Degradation (%) at t = 27h` = 
                                   anno_barplot(heatmap_anno$degr, gp = gpar(fill = "grey30")), 
                                 annotation_name_gp = gpar(fontsize = 14, fontface = "bold"))

column_groups <- gsub("\\s*\\+.*", "", heatmap_anno$interaction1)
column_groups <- factor(column_groups, levels = compound_order)

pdf_file <- paste0(wd_path, "/Metabolites_collection/7_Metabolite_map_by_compound.pdf")
pdf(pdf_file, width = 12, height = 7.7)

ht <- Heatmap(as.matrix(heatmap),
              col = c("grey90","#440154FF"), 
              top_annotation = degradation,
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              show_row_names = TRUE, 
              show_column_names = TRUE, 
              column_split = column_groups,
              column_title = NULL,
              border = TRUE,
              row_gap = unit(5, "mm"),
              rect_gp = gpar(col = "white", lwd = 1.5),
              #column_names_rot = 45,
              row_names_gp = gpar(fontsize = 10, hjust = 1, vjust = 1),
              column_names_gp = gpar(fontsize = 10, hjust = 1, vjust = 1),
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 14, fontface = "bold"),
                legend_height = unit(3, "cm"),
                legend_width = unit(4, "cm"),
                title = "",
                labels = legend_labels))

draw(ht, heatmap_legend_side = "top", padding = unit(c(60, 20, 5, 10), "mm"))

dev.off()

#### S.herbicidovorans data for GoF ----
sherbi <- metabolites %>%  
  filter(species == "S.herbicidovorans") %>% 
  mutate(RT = RT + 0.1) %>% # for the ms method 
  mutate(metabolite = paste0(Compound, "_", MZ, "_", RT)) %>% 
  select(species, Compound, metabolite, MZ, RT) %>% 
  arrange(MZ) %>% 
  distinct() %>% 
  write_csv(paste0(wd_path, "/Metabolites_collection/7_GoF_metabolite_list_Sherbicidovorans.csv"))

