rm(list = ls())

#### Import libraries ----
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales) 
library(readr)
library(rstatix)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))}

#### SET WD ----

# Set paths to working directory
wd_path <- getwd()
setwd(wd_path)
data_path <- paste0(wd_path, "/Input/")

#### Extract IDs ----
thr = 1e-6

id_A0A086P4J8 <- read_csv(paste0(data_path, "A0A086P4J8.csv")) 

id_A0A086P4J8 %>% pull(`Target Name`) %>% write_lines(paste0(data_path, "A0A086P4J8_IDs.txt")) # save all IDs

id_A0A086P4J8 <- id_A0A086P4J8 %>% 
  filter(`E-value` <= thr) %>% 
  pull(`Target Name`)
  
id_B0NC66 <- read_csv(paste0(data_path, "B0NC66.csv")) 

id_B0NC66 %>% pull(`Target Name`) %>% write_lines(paste0(data_path, "B0NC66_IDs.txt"))

id_B0NC66 <- id_B0NC66 %>% 
  filter(`E-value` <= thr) %>% 
  pull(`Target Name`)
  
#### READ DATA ----
df_A0A086P4J8 <- read_tsv(
  paste0(data_path, "A0A086P4J8_biomes.tsv"), 
  col_names = c("MGYP", "count", "biome")) %>% 
  filter(MGYP%in%id_A0A086P4J8) %>% 
  mutate(protein = "A0A086P4J8") # biomes associated with each sequence and respective protein count

df_B0NC66 <- read_tsv(
  paste0(data_path, "B0NC66_biomes.tsv"),
  col_names = c("MGYP", "count", "biome")) %>% 
  filter(MGYP%in%id_B0NC66) %>% 
  mutate(protein = "B0NC66")

df_combined_raw <- bind_rows(df_A0A086P4J8, df_B0NC66)

df_biome_count <- read_tsv(
  paste0(data_path, "mgy_biome_counts.tsv"),
  col_names = c("biome_count", "biome")) # observation counts of each biome

#### MERGE AND UNIFY DATA ----
df_combined <- df_combined_raw %>%
  left_join(df_biome_count) %>% 
  separate(biome, into = c("root", "type", "level1", "level2", "level3"), sep = ":") %>% 
  mutate(biome_group = case_when(
    level1 %in% c("Human", "Mammals") ~ "Gut-associated host",
    level1 %in% c("Fish", "Birds", "Insecta", "Mollusca", "Amphibia") ~ "Other host-associated",
    level1 %in% c("Aquatic", "Algae") ~ "Aquatic environment",
    level1 %in% c("Terrestrial", "Plants") ~ "Terrestrial environment",
    level1 %in% c("Wastewater", "Solid waste", "Bioremediation") ~ "Waste / polluted environment",
    level1 %in% c("Biogas plant", "Bioreactor", "Food production",
                  "Built environment", "Lab enrichment",
                  "Lab Synthesis", "Microbial", "Modeled") ~ "Engineered / laboratory",
    TRUE ~ "Other / unknown")) %>% 
  mutate(biome_group = factor(biome_group, levels = c("Gut-associated host", "Other host-associated", 
                                                      "Aquatic environment", "Terrestrial environment",
                                                      "Waste / polluted environment", "Engineered / laboratory",
                                                      "Other / unknown"))) 

biome_totals_grouped <- df_combined %>%
  distinct(biome_group, biome_count) %>%
  group_by(biome_group) %>%
  summarise(biome_total_count = sum(biome_count), .groups = "drop")

#### RELATIVE PREVALENCE ----
# Where do homologs of these two proteins tend to occur?
plot_df <- df_combined %>%
  group_by(protein, biome_group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>% 
  group_by(protein) %>%
  mutate(percent = 100 * total_count / sum(total_count)) %>%
  ungroup()

plot <- ggplot(plot_df, 
               aes(x = biome_group, 
                   y = percent, 
                   fill = biome_group)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ protein, ncol = 1) +
  scale_fill_viridis_d(option = "D", begin = 0.05, end = 0.9) +
  labs(
    x = NULL,
    y = "% homolog observations across environmets",
    fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    strip.text   = element_text(size = 24, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")

plot

ggsave(plot,
       filename = file.path(wd_path, "Plots", "Percent_of_homolog_observations.pdf"),
       width = 14,
       height = 10,
       units = "in")

plot <- ggplot(
  plot_df,
  aes(
    x = protein,
    y = percent,
    fill = biome_group
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  scale_fill_viridis_d(option = "D", begin = 0.05, end = 0.9) +
  labs(
    x = NULL,
    y = "% homolog observations across environments",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

plot

ggsave(plot,
       filename = file.path(wd_path, "Plots", "Percent_of_homolog_observations_stack.pdf"),
       width = 10,
       height = 10,
       units = "in")

#### ABSOLUTE PREVALENCE ----
plot_df_abs <- df_combined %>%
  group_by(protein, biome_group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  left_join(biome_totals_grouped, by = "biome_group") %>%
  mutate(enrichment = total_count / biome_total_count)

plot_abs <- ggplot(plot_df_abs, aes(
  x = biome_group, enrichment,
  y = enrichment,
  fill = biome_group)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ protein, ncol = 1) +
  scale_fill_viridis_d(option = "D", begin = 0.05, end = 0.9) +
  labs(
    x = NULL,
    y = "Normalized enrichment",
    fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    strip.text   = element_text(size = 24, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")

plot_abs

ggsave(plot_abs,
       filename = file.path(wd_path, "Plots", "Normalized_enrichment.pdf"),
       width = 14,
       height = 10,
       units = "in")

plot_abs <- ggplot(
  plot_df_abs,
  aes(
    x = protein,
    y = enrichment,
    fill = biome_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  facet_wrap(~ biome_group, scales = "free_y") +
  scale_fill_viridis_d(option = "D", begin = 0.05, end = 0.9) +
  scale_y_continuous(labels = scientific_10) +
  labs(
    x = NULL,
    y = "Normalized enrichment",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    strip.text   = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")

plot_abs

ggsave(plot_abs,
       filename = file.path(wd_path, "Plots", "Normalized_enrichment_facet.pdf"),
       width = 14,
       height = 10,
       units = "in")

heatmap <- ggplot(
  plot_df_abs,
  aes(
    x = protein,
    y = biome_group,
    fill = enrichment
  )
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = signif(enrichment, 2)),
    size = 5,
    color = "black"
  ) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log10",
    labels = scientific_10
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Normalized\nenrichment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_blank()
  )

heatmap

ggsave(heatmap,
       filename = file.path(wd_path, "Plots", "Normalized_enrichment_heatmap.pdf"),
       width = 10,
       height = 10,
       units = "in")
