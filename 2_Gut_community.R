rm(list = ls())

#### Import libraries ----
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales) 

scientific_10 <- function(x) {   
  parse(text = gsub("e\\+*", " %*% 10^", scientific_format()(x)))}

#### SET WD ----

# Set paths to working directory
wd_path <- getwd()
setwd(wd_path)
data_path <- paste0(wd_path, "/Input/")

#### DATA ----
data <- read_csv(paste0(data_path, "Gut_community.csv"))

data$condition <- factor(data$condition, levels = c("Comm1","Comm2","Comm3","Comm4", "GMM"))

#### Dexamethasone-phosphate ----
DP <- data %>% 
  filter(compound == 'Dexamethasone-phosphate') %>% 
  group_by(condition, time, metabolite) %>%  
  summarize(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE),
            .groups   = "drop")

DP$metabolite <- factor(DP$metabolite, 
                        levels =c("Dexamethasone-phosphate", 
                                  "Dexamethasone", 
                                  "Dexamethasone-17-Ketone"))

com_conditions <- unique(DP$condition[DP$condition != "GMM"])

DP_plot <- bind_rows(
  DP %>% filter(condition != "GMM"),
  DP %>%
    filter(condition == "GMM") %>%
    select(-condition) %>%
    crossing(condition = com_conditions)) %>%
  mutate(sample = ifelse(condition == "GMM", "sterile control", "community"))

DP_plot <- bind_rows(
  DP %>% filter(condition != "GMM") %>% mutate(sample = "community"),
  DP %>%
    filter(condition == "GMM") %>%
    select(-condition) %>%
    crossing(condition = com_conditions) %>%
    mutate(sample = "sterile control"))

p <- DP_plot %>%
  ggplot(aes(x = time, y = mean_area,
             color = condition,
             fill = condition,
             linetype = sample,
             group = interaction(condition, sample))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c(
    "community" = "solid",
    "sterile control" = "dashed")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 1) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

plot(p)

ggsave(p, filename = file.path(wd_path, "Plots", "2_Dexamethasone_phosphate_communities.pdf"),
       width    = 8,
       height   = 10,
       units    = "in")

#### Fludrocortisone-acetate ----
FA <- data %>% 
  filter(compound == 'Fludrocortisone-acetate') %>% 
  group_by(condition, time, metabolite) %>%  
  summarize(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE),
            .groups   = "drop") %>% 
  filter(metabolite != "Fludrocortisone_319.1704_2.72") %>% 
  filter(metabolite != "Fludrocortisone_381.2069_2.3")

FA$metabolite <- gsub("Fludrocortisone_321.1861_2.82", "Fludrocortisone-17-Ketone", FA$metabolite)
FA$metabolite <- gsub("Fludrocortisone_383.2228_2.3", "Dihydro-fludrocortisone", FA$metabolite)
FA$metabolite <- gsub("Fludrocortisone_323.2017_2.79", "Dihydro-fludrocortisone-17-Ketone", FA$metabolite)

FA$metabolite <- factor(FA$metabolite, 
                        levels =c("Fludrocortisone-acetate", 
                                  "Fludrocortisone", "Dihydro-fludrocortisone", 
                                  "Fludrocortisone-17-Ketone", "Dihydro-fludrocortisone-17-Ketone"))

com_conditions <- unique(FA$condition[FA$condition != "GMM"])

FA_plot <- bind_rows(
  FA %>% filter(condition != "GMM"),
  FA %>%
    filter(condition == "GMM") %>%
    select(-condition) %>%
    crossing(condition = com_conditions)) %>%
  mutate(sample = ifelse(condition == "GMM", "sterile control", "community")) 

FA_plot <- bind_rows(
  FA %>% filter(condition != "GMM") %>% mutate(sample = "community"),
  FA %>%
    filter(condition == "GMM") %>%
    select(-condition) %>%
    crossing(condition = com_conditions) %>%
    mutate(sample = "sterile control")) %>% 
  filter(metabolite != "Dihydro-fludrocortisone-17-Ketone")
  

p <- FA_plot %>%
  ggplot(aes(x = time, y = mean_area,
             color = condition,
             fill = condition,
             linetype = sample,
             group = interaction(condition, sample))) +
  geom_ribbon(aes(ymin = mean_area - sd_area,
                  ymax = mean_area + sd_area),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c(
    "community" = "solid",
    "sterile control" = "dashed")) +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0, end = 0.85) +
  scale_y_continuous(labels = scientific_10) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time, h",
    y = "Peak Area",
    color = "Metabolite",
    linetype = "Sample") +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 4) +
  theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text   = element_text(size = 18, color = "black"),
    strip.text  = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 24),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

plot(p)

ggsave(p, filename = file.path(wd_path, "Plots", "2_Fludrocortisone_acetate_communities.pdf"),
       width    = 24,
       height   = 4,
       units    = "in")

