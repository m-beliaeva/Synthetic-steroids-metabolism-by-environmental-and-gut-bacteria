rm(list = ls())

#### Import libraries ----
library(tidyverse)
library(scales)
library(viridis)

scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))}

#### SET PATH ----
wd_path <- getwd()
setwd(wd_path)

#### Beclomethasone dipropionate ----
df <- read_csv(paste0(wd_path, "/Input/Sherbicidovorans_Beclomethasone-dipropionate.csv")) %>% 
  rename('Drug' = name,
         'Condition' = condition) %>% 
  group_by(time, Condition, metabolite) %>% 
  summarise(mean_area = mean(area), sd_area = sd(area), .groups = "drop") %>% 
  mutate(metabolite = factor(metabolite, 
                             levels = c("Beclomethasone dipropionate", 
                                        "Beclomethasone-17-propionate", 
                                        "Beclomethasone-21-propionate",
                                        "Beclomethasone")))
df %>% 
  ggplot(aes(x = time, y = mean_area,
             color = metabolite,
             fill = metabolite,
             linetype = Condition,
             group = interaction(metabolite, Condition))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_linetype_manual(values = c(
    "S.herbicidovorans" = "solid",
    "sterile control" = "dotted")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
  #facet_wrap(~ metabolite, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/Sherbicidovorans.pdf"),
       width = 12,
       height = 4.5,
       unit = "in",
       dpi = 300,  
       device = cairo_pdf)

### Fludrocortisone acetate ----
df <- read_csv(paste0(wd_path, "/Input/Fludrocortisone-acetate.csv")) %>% 
  rename('Drug' = name,
         'Condition' = condition) %>% 
  group_by(Condition, time, metabolite) %>%
  summarize(mean_area = mean(area),
            sd_area = sd(area), .groups = 'drop') %>% 
  mutate(metabolite = factor(metabolite, 
                             levels = c("Fludrocortisone-acetate", "Fludrocortisone")))

df %>% 
  ggplot(aes(x = time, y = mean_area,
             color = metabolite,
             fill = metabolite,
             linetype = Condition,
             group = interaction(metabolite, Condition))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_linetype_manual(values = c(
    "B.uniformis" = "solid",
    "sterile control" = "dotted")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
  #facet_wrap(~ metabolite, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/Buniformis.pdf"),
       width = 12,
       height = 4.5,
       units = "in",
       dpi = 300,
       device = cairo_pdf)

### Dexmethasone phosphate ----
df <- read_csv(paste0(wd_path, "/Input/Dexmethasone-phosphate.csv")) %>% 
  rename('Drug' = name,
         'Condition' = condition) %>% 
  group_by(Condition, time, metabolite) %>%
  summarize(mean_area = mean(area),
            sd_area = sd(area), .groups = 'drop') %>% 
  mutate(metabolite = factor(metabolite, 
                             levels = c("Dexamethasone-phosphate", "Dexamethasone")))
df %>% 
  ggplot(aes(x = time, y = mean_area,
             color = metabolite,
             fill = metabolite,
             linetype = Condition,
             group = interaction(metabolite, Condition))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_linetype_manual(values = c(
    "P.vulgatus" = "solid",
    "sterile control" = "dotted")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
  #facet_wrap(~ metabolite, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/Pvulgatus.pdf"),
       width = 12,
       height = 4.5,
       units = "in",
       dpi = 300,
       device = cairo_pdf)

### Cscindens metabolites ----
df <- read_csv(paste0(wd_path, "/Input/Cscindens.csv")) %>% 
  rename('Drug' = name,
         'Condition' = condition) %>% 
  separate(col = metabolite, into = c("Compound"), sep = "-", remove = FALSE) %>% 
  mutate(Compound = gsub("12", "12-Dihydro-dexamethasone", Compound)) %>% 
  group_by(Condition, time, metabolite, Compound) %>%
  summarize(mean_area = mean(area),
            sd_area = sd(area), .groups = 'drop') %>% 
  mutate(type = ifelse(grepl("ketone", metabolite, ignore.case = TRUE),
                       "ketone",
                       "parent")) %>% 
  mutate(type = factor(type, levels = c('parent', 'ketone')))

df %>% 
  ggplot(aes(x = time, y = mean_area,
             color = type,
             fill = type,
             linetype = Condition,
             group = interaction(metabolite, Condition))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_linetype_manual(values = c(
    "C.scindens" = "solid",
    "sterile control" = "dotted")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10, limits = c(0, NA)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  facet_wrap(~ Compound, scales = "free_y", nrow = 5) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

ggsave(filename = paste0(wd_path, "/Plots/Cscindens_v2.pdf"),
       width = 10,
       height = 20,
       units = "in",
       dpi = 300,
       device = cairo_pdf)

