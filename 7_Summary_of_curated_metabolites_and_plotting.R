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

#### PARENT DATA ----
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

#### METABOLITE DATA -----
# Get list of curated metabolites
metabolites <- read_csv(paste0(wd_path, "/Metabolites_collection/6_Metabolites_filtered_with_reactions_curated.csv")) %>% 
  filter(is.na(parent_match)) %>% 
  filter(metabolite_type != "in-source defluorination") # 167 curated non-unique metabolites detected excluding parent hits

metabolites$Compound <- gsub("12-Dihydro-dexamethasone", "1,2-Dihydro-dexamethasone", metabolites$Compound)
metabolites$Compound <- gsub("67-Dehydro-ethinylestradiol", "6,7-Dehydro-ethinylestradiol", metabolites$Compound)
metabolites$Compound <- gsub("67-Dehydro-ethynylestradiol", "6,7-Dehydro-ethinylestradiol", metabolites$Compound)

# Tolerance for MZ_RT consensus
abs_mass = 0.001 # absolute mz tolerance
rel_mass = 1*10^-6 # relative mz tolerance (10 ppm)
abs_rt = 0.1 # absolute rt tolerance

# Cluster metabolites based on MZ_RT combinations by tolerance
metabolites_clustered <- metabolites %>%
  arrange(MZ, RT) %>%
  mutate(cluster_id = NA_integer_)

cluster <- 1
metabolites_clustered$cluster_id[1] <- cluster

for (i in 2:nrow(metabolites_clustered)) {
  
  mz_prev <- metabolites_clustered$MZ[i - 1]
  mz_curr <- metabolites_clustered$MZ[i]
  
  rt_prev <- metabolites_clustered$RT[i - 1]
  rt_curr <- metabolites_clustered$RT[i]
  
  mz_tol <- max(abs_mass, rel_mass * mz_prev)
  
  if (abs(mz_curr - mz_prev) <= mz_tol &&
      abs(rt_curr - rt_prev) <= abs_rt) {
    
    metabolites_clustered$cluster_id[i] <- cluster
    
  } else {
    cluster <- cluster + 1
    metabolites_clustered$cluster_id[i] <- cluster
  }
}

cluster_consensus <- metabolites_clustered %>%
  group_by(cluster_id) %>%
  summarise(
    consensus_MZ = mean(MZ, na.rm = TRUE),
    consensus_RT = mean(RT, na.rm = TRUE),
    .groups = "drop"
  )

metabolites <- metabolites_clustered %>%
  left_join(cluster_consensus, by = "cluster_id") %>%
  mutate(
    consensus_MZ_RT = paste0(
      round(consensus_MZ, 4), "_",
      round(consensus_RT, 1)
    )
  ) %>% 
  select(-consensus_MZ, -consensus_RT, -cluster_id)

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
  
  # Unique metabolites (Compound–MZ–RT)
  count_by_id(
    metabolites,
    "Unique metabolites (Compound-MZ-RT)",
    is.na(parent_match),
    paste(Compound, consensus_MZ_RT)),
  
  # Unique compounds
  count_by_id(
    metabolites,
    "Parent compounds",
    is.na(parent_match),
    Compound),
  
  # Unique species
  count_by_id(
    metabolites,
    "Species",
    is.na(parent_match),
    species),
  
  # Compound–species combinations
  count_by_id(
    metabolites,
    "Compound-species combinations",
    is.na(parent_match),
    paste(Compound, species)
  ),
  
  # Assigned reactions (Compound–reaction)
  count_by_id(
    metabolites,
    "Assigned reactions and their combinations (Compound-reaction)",
    !is.na(reaction),
    paste(Compound, reaction)
  ),
  
  # Compound-species combinations with assigned reactions
  count_by_id(
    metabolites,
    "Compound-species combinations with assigned reactions",
    is.na(parent_match) & !is.na(reaction),
    paste(Compound, species)
  )
)

# Add confidently identified structures
cosmic80_row <- tibble(
  description = "Unique metabolites with COSMIC_STRUCTURE_SCORE >= 80",
  MSMS = metabolites %>%
    filter(is.na(parent_match), identified_by == "MSMS", COSMIC_STRUCTURE_SCORE >= 80) %>%
    select(Compound, consensus_MZ_RT) %>%
    nrow(),
  
  manual_curation = NA,
  
  total = metabolites %>%
    filter(
      is.na(parent_match),
      identified_by == "MSMS",
      COSMIC_STRUCTURE_SCORE >= 80) %>%
    select(Compound, consensus_MZ_RT) %>%
    nrow())

summary_table <- bind_rows(summary_table, cosmic80_row)

# Add total biotransformation number
total <- tibble(
  description = "Total biotransfomrations detected",
  MSMS = metabolites %>%
    filter(is.na(parent_match), identified_by == "MSMS") %>%
    nrow(),
  
  manual_curation = metabolites %>%
    filter(is.na(parent_match), identified_by == "manual_curation") %>%
    nrow(),
  
  total = metabolites %>%
    filter(is.na(parent_match), identified_by == "MSMS" | identified_by == "manual_curation") %>%
    nrow())

summary_table <- bind_rows(summary_table, total)

# Add total biotransformation number with assigned reactions
total_reactions <- tibble(
  description = "Biotransfomrations with assigned reactions",
  MSMS = metabolites %>%
    filter(is.na(parent_match), !is.na(reaction), identified_by == "MSMS", ) %>%
    nrow(),
  
  manual_curation = metabolites %>%
    filter(is.na(parent_match), !is.na(reaction), identified_by == "manual_curation") %>%
    nrow(),
  
  total = metabolites %>%
    filter(is.na(parent_match), !is.na(reaction), identified_by == "MSMS" | identified_by == "manual_curation") %>%
    nrow())

summary_table <- bind_rows(summary_table, total_reactions)

write_csv(summary_table, paste0(wd_path, "/Metabolites_collection/7_Summary_of_unique_metabolites.csv"))

#### BARPLOT METABOLITES PER STEROID ----
df_plot <- metabolites %>%
  filter(is.na(parent_match)) %>% 
  mutate(MetaboliteID = consensus_MZ_RT)

# Count unique metabolites per Compound and reaction status
metabolite_counts <- df_plot %>%
  group_by(Compound, reaction_status = ifelse(is.na(reaction), "No", "Yes")) %>%
  summarise(unique_metabolites = n_distinct(MetaboliteID), .groups = "drop")%>% 
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
  scale_fill_viridis_d(option = "D") +
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
  scale_fill_viridis_d(option = "D") +
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
# Prepare and combine metabolite and parent data
metabolites_df <- metabolites %>% 
  separate(reaction, into = c("reaction1", "reaction2", "reaction3"), sep = " & ") %>% 
  pivot_longer(cols = c("reaction1", "reaction2", "reaction3"), values_to = "reaction") %>% 
  select(-c(name, SMILES_PARENT, FORMULA_PARENT, MZ_PARENT, RT_PARENT, MZ, RT, delta_mass, atomic_difference, 
            SIRIUS_FORMULA, SIRIUS_FORMULA_SCORE, CSIFingerID_STRUCTURE, COSMIC_STRUCTURE_SCORE, METABOLITE,
            identified_by, metabolite_type, SIRIUS_PARAMETERS, parent_match, SIRIUS_VERSION, consensus_MZ_RT)) %>% 
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
  mutate(RT = RT + 0.1) %>% #for the ms method 
  mutate(metabolite = paste0(Compound, "_", MZ, "_", RT)) %>% 
  select(species, Compound, metabolite, MZ, RT) %>% 
  arrange(MZ) %>% 
  distinct() %>% 
  write_csv(paste0(wd_path, "/Metabolites_collection/7_GoF_metabolite_list_Sherbicidovorans.csv"))

