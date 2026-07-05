rm(list = ls())

#### Import libraries ----
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales) 
library(rstatix)

#### THEME PUBLICATION ----
scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))}

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            
            # Legend settings
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1.3, "cm"),
            legend.key.height = unit(0.4, "cm"),
            legend.margin = unit(1, "cm"),
            legend.title = element_blank(),
            
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#### SET WD ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)
data_path <- paste0(wd_path, "/Input/")

### READ DATA ----
data <- read.csv(paste0(data_path, 'Additional_compounds.csv')) 

data %>% filter(!str_detect(metabolite, "parent")) %>% select(mz, rt) %>% distinct() %>% nrow() # Unique metabolites 

### CALCULATE PVALS ----
data <- data %>% 
  mutate(time = as.factor(time),
         name = as.factor(name),
         metabolite = as.factor(metabolite)) 

met_list <- unique(data$metabolite)
data_ctrl <- data %>% 
  filter(condition == 'ctrl')

# Two-way interactions between condition and time
df_anova <- data.frame()
df_anova_list <- list()
df_anova_list[[1]] <- df_anova

for (met in met_list) {
  tmp <- data %>% 
    filter(metabolite == met)
  
  test <- tmp %>% 
    anova_test(area ~ condition*time) # add paired here
  
  pval <- get_anova_table(test)$p[3] #%>% pull(p)
  
  df_anova <- tmp %>% 
    filter(condition != "ctrl") %>% 
    mutate(pval = pval)
  
  df_anova_list[[length(df_anova_list)+1]] <- df_anova
}

data_anova <- bind_rows(df_anova_list)

# Adjust p-values with fdr
data_padj <- data_anova %>% 
  select(-c(time, area)) %>% 
  distinct()

data_padj$padj <- p.adjust(data_padj$pval, method = "fdr")

data_stats <- full_join(data_anova, data_padj, 
                        by = c("name", "condition", "metabolite", "mz", "rt", "pval")) %>% 
  bind_rows(data_ctrl)

data_stats$condition <- gsub('Sh', 'S.herbicidovorans', data_stats$condition)

# Add reactions
reactions <- read.csv(paste0(data_path, 'Metabolite_list.csv')) %>% 
  select(metabolite, delta_mass, reaction)

data_long <- data_stats %>%
  left_join(reactions) %>% 
  mutate(condition = factor(condition, levels = c("S.herbicidovorans", "ctrl"))) %>%
  mutate(time = as.numeric(as.character(time))) %>% 
  group_by(name, metabolite, mz, rt, condition, time, reaction, pval, padj) %>%
  summarise(mean_area = mean(area), sd = sd(area), .groups = "drop") %>% 
  write_csv(paste0(wd_path, '/Output/1_Additional_compounds_analysis.csv')) 

metabolites <- data_long %>% 
  filter(reaction != 'parent') %>% 
  select(mz, rt) %>% 
  distinct() %>% 
  nrow() # 32

### PLOT ALL METABOLITES ----
for (j in unique(data_long$name)) {
  
  plot <- data_long %>%
    filter(name == j)
  
  plot %>%
    ggplot(aes(x = time, y = mean_area,
               linetype = condition,
               group = interaction(metabolite, condition))) +
    geom_ribbon(aes(ymin = mean_area - sd,
                    ymax = mean_area + sd),
                fill = "#440154", alpha = 0.2, color = NA) +
    geom_line(linewidth = 1, color = "#440154") +
    scale_linetype_manual(values = c(
      "S.herbicidovorans" = "solid",
      "ctrl" = "dotted")) +
    scale_x_continuous(labels = label_number(1)) +
    scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
    facet_wrap(~ metabolite, scales = "free_y", ncol = 1) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Time, h",
      y = "Peak Area",
      linetype = "Sample") +
    theme(
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),
      axis.text   = element_text(size = 18, color = "black"),
      strip.text  = element_text(size = 20, face = "bold"),
      legend.title = element_blank(),
      legend.text  = element_text(size = 24),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))
  
  ggsave(filename = paste0(wd_path, "/Plots/1_", j, ".pdf"),
         width = 10,
         height = 10*(nrow(plot)/80)+1,
         units = "in",
         dpi = 300,
         device = cairo_pdf)
}

### PLOT ALL METABOLITES (SAME PLOT, SCALED) ----
for (j in unique(data_long$name)) {
  
  plot <- data_long %>%
    filter(name == j) %>% 
    group_by(name,metabolite,mz,rt) %>%
    mutate(max_area = max(mean_area)) %>% 
    ungroup() %>% 
    mutate(norm_area = mean_area/max_area, norm_sd = sd/max_area)
  
  plot %>%
    ggplot(aes(
      x = time,
      y = norm_area,
      color = metabolite,
      fill = metabolite,
      linetype = condition,
      group = interaction(metabolite, condition))) +
    geom_ribbon(
      aes(ymin = norm_area - norm_sd, ymax = norm_area + norm_sd), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(
      option = "D",
      direction = -1,
      end = 0.85) +
    scale_fill_viridis_d(
      option = "D",
      direction = -1,
      end = 0.85) +
    scale_linetype_manual(values = c(
      "S.herbicidovorans" = "solid",
      "ctrl" = "dotted")) +
    scale_x_continuous(labels = label_number(1)) +
    scale_y_continuous(labels = scientific_10) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_minimal(base_size = 14) +
    labs(x = "Time, h", y = "Peak Area", color = "Metabolite", fill = "Metabolite", linetype = "Sample") +
    theme(
      axis.title.x = element_text(size = 20, face = "bold", color = "black"),
      axis.title.y = element_text(size = 20, face = "bold", color = "black"),
      axis.text = element_text(size = 18, color = "black"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 24),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.5))
  
  ggsave(filename = paste0(wd_path, "/Plots/1_", j, "_scaled.pdf"),
         width = 24,
         height = 10,
         units = "in",
         dpi = 300,
         device = cairo_pdf)
}

### PLOT PARENT COMPOUNDS ----
parents <- data_long %>% 
  filter(reaction == 'parent')

parents_plot <- parents %>%
  ggplot(aes(x = time, y = mean_area,
             linetype = condition,
             group = interaction(metabolite, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd,
                  ymax = mean_area + sd),
              fill = "#440154", alpha = 0.2, color = NA) +
  geom_line(linewidth = 1, color = "#440154") +
  scale_linetype_manual(values = c(
    "S.herbicidovorans" = "solid",
    "ctrl" = "dotted")) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
  facet_wrap(vars(name), scales = "free", ncol = 4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(parents_plot,
       filename = paste0(wd_path, "/Plots/1_All_parent_compounds.pdf"),
       width = 20,
       height = 9,
       unit = "in",
       dpi = 200,
       device = cairo_pdf)

### PLOT SIDE-CHAIN CLEAVED METABOLITES ----
sidechain <- data_long %>%
  filter(metabolite != "Corticosterone_303.1955_2.6") %>% # Select metabolite with higher intensity
  filter(metabolite != "Fludrocortisone_321.1861_2.58") %>% # Select metabolite with higher intensity
  filter(!name %in% c("12-Didehydrocorticosterone", "Fludrocortisone-acetate", "Prednisolone", "Dihydrocortisol")) %>% # Select compounds no keto metabolite
  filter(reaction %in% c(
    "parent",
    "side-chain cleavage* [-C2H4O2]",
    "side-chain cleavage [-C2H4O2]")) %>%
  mutate(metabolite_type = case_when(
    reaction == "parent" ~ "Parent",
    TRUE ~ "17-Ketone"))

sidechain_plot <- sidechain %>%
  ggplot(aes(x = time, y = mean_area, color = metabolite_type, 
             fill = metabolite_type, linetype = condition,
             group = interaction(metabolite_type, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd,
                  ymax = mean_area + sd),
              alpha = 0.2,
              color = NA) +
  geom_line(size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_color_viridis_d(option = "D", direction = -1, end = 0.85) +
  scale_fill_viridis_d(option = "D", direction = -1, end = 0.85) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(sidechain_plot,
       filename = paste0(wd_path, "/Plots/1_Ketone_metabolites.pdf"),
       width = 20,
       height = 9,
       unit = "in",
       dpi = 200,
       device = cairo_pdf)

sidechain_plot2 <- sidechain %>%
  filter(reaction != "parent") %>% 
  ggplot(aes(x = time, y = mean_area, color = metabolite_type, 
             fill = metabolite_type, linetype = condition,
             group = interaction(metabolite_type, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd,
                  ymax = mean_area + sd),
              alpha = 0.2,
              color = NA) +
  geom_line(size = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_color_viridis_d(option = "D", direction = -1, end = 0.85) +
  scale_fill_viridis_d(option = "D", direction = -1, end = 0.85) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(sidechain_plot2,
       filename = paste0(wd_path, "/Plots/1_Ketone_metabolites_only.pdf"),
       width = 20,
       height = 5.8,
       unit = "in",
       dpi = 200,
       device = cairo_pdf)

### FLUDROCORTISONE ----
plot <- data_long %>%
  filter(name == "Fludrocortisone") %>% 
  mutate(metabolite = factor(metabolite)) 

plot %>%
  ggplot(aes(
    x = time,
    y = mean_area,
    color = metabolite,
    fill = metabolite,
    linetype = condition,
    group = interaction(metabolite, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd, ymax = mean_area + sd), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(
    option = "D",
    direction = -1,
    end = 0.85) +
  scale_fill_viridis_d(
    option = "D",
    direction = -1,
    end = 0.85) +
  scale_linetype_manual(values = c(
    "S.herbicidovorans" = "solid",
    "ctrl" = "dotted")) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 14) +
  labs(x = "Time, h", y = "Peak Area", color = "Metabolite", fill = "Metabolite", linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.text = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/1_Fludrocortisone_SI.pdf"),
       width = 15,
       height = 7,
       units = "in",
       dpi = 300,
       device = cairo_pdf)

### HYDROCORTISONE ----
plot <- data_long %>%
  filter(name == "Hydrocortisone") %>% 
  mutate(metabolite = factor(metabolite)) 

plot %>%
  ggplot(aes(
    x = time,
    y = mean_area,
    color = metabolite,
    fill = metabolite,
    linetype = condition,
    group = interaction(metabolite, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd, ymax = mean_area + sd), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(
    option = "D",
    direction = -1,
    end = 0.85) +
  scale_fill_viridis_d(
    option = "D",
    direction = -1,
    end = 0.85) +
  scale_linetype_manual(values = c(
    "S.herbicidovorans" = "solid",
    "ctrl" = "dotted")) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 14) +
  labs(x = "Time, h", y = "Peak Area", color = "Metabolite", fill = "Metabolite", linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.text = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/1_Hydrocortisone_SI.pdf"),
       width = 15,
       height = 2.8,
       units = "in",
       dpi = 300,
       device = cairo_pdf)

### CORTICOSTERONE ----
plot <- data_long %>%
  filter(name == "Corticosterone") 
  # filter(metabolite %in% c(
  #   "Corticosterone-parent_347.2217_2.85",
  #   "Corticosterone_303.1955_2.8",
  #   "Corticosterone_363.2166_2.93"
  # )) %>%
  # mutate(
  #   metabolite = factor(
  #     metabolite,
  #     levels = c(
  #       "Corticosterone-parent_347.2217_2.85",
  #       "Corticosterone_363.2166_2.93",
  #       "Corticosterone_303.1955_2.8"))) %>% 
  # filter(time <= 8)

plot %>%
  ggplot(aes(
    x = time,
    y = mean_area,
    color = metabolite,
    fill = metabolite,
    linetype = condition,
    group = interaction(metabolite, condition))) +
  geom_ribbon(aes(ymin = mean_area - sd, ymax = mean_area + sd), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(
    option = "D",
    direction = 1,
    end = 0.85) +
  scale_fill_viridis_d(
    option = "D",
    direction = 1,
    end = 0.85) +
  scale_linetype_manual(values = c(
    "S.herbicidovorans" = "solid",
    "ctrl" = "dotted")) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = label_number(1)) +
  scale_y_continuous(labels = scientific_10) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 14) +
  labs(x = "Time, h", y = "Peak Area", color = "Metabolite", fill = "Metabolite", linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.text = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/1_Corticosterone_SI.pdf"),
       width = 15,
       height = 6,
       units = "in",
       dpi = 300,
       device = cairo_pdf)
