rm(list = ls())

#### Import libraries ----
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales) 

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
            legend.title = element_text(face="italic"),
            
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#### SET WD ----

# Set paths to working directory and mgf folder
wd_path <- getwd()
setwd(wd_path)
data_path <- paste0(wd_path, "/Input/")

#### Fludrocortisone acetate ----
FA_raw <- read_csv(paste0(data_path, "Fludrocortisone-acetate.csv")) %>% 
  mutate(step = as.factor(ifelse(nchar(condition) >= 3 & nchar(condition) <= 4, 1, 
                                 ifelse(nchar(condition) >= 7 & nchar(condition) <= 9, 2, 
                                        ifelse(nchar(condition) >= 11 & nchar(condition) <= 14, 3, 0))))) %>% 
  filter(step == 2) %>% 
  select(-c(step, compound, replicate)) %>%
  filter(time != 'filtrate') %>% 
  pivot_longer(cols = c(3:6),
               names_to = "metabolite",
               values_to = "area") %>% 
  filter(metabolite != "Fludrocortisone-319.1704")

FA_raw$condition <- gsub('Csci', 'C.scindens', FA_raw$condition)
FA_raw$condition <- gsub('Buni', 'B.uniformis', FA_raw$condition)

FA_raw$metabolite <- factor(FA_raw$metabolite, 
                            levels =c("Fludrocortisone-acetate", 
                                      "Fludrocortisone-381.2072", 
                                      "Fludrocortisone-321.1860"))   

FA <- FA_raw %>% 
  group_by(condition, time, metabolite) %>%  
  mutate(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE),
            max_area  = max(area, na.rm = TRUE),
            .groups   = "drop") %>% 
  mutate(norm_area    = ifelse(max_area == 0, 0, mean_area / max_area),
         sd_norm_area = ifelse(max_area == 0, 0, sd_area / max_area)) 

FA <- FA %>% 
  separate(condition, into = c("condition1","condition2"), sep = "-")

FA$condition1 <- factor(FA$condition1, 
                        levels =c("GMM",  'C.scindens', 'B.uniformis'))  

FA$condition2 <- factor(FA$condition2, 
                        levels =c("GMM",  'C.scindens', 'B.uniformis'))  

FA %>% 
  ungroup() %>%    
  mutate(time = factor(time)) %>%  
  ggplot(mapping = aes(x = time, y = mean_area, fill = metabolite, group = metabolite)) +  
  geom_bar(stat = "identity", position = position_dodge(1)) +   
  geom_errorbar(aes(ymin = mean_area - sd_area, ymax = mean_area + sd_area),                  
                position = position_dodge(1), width = 0.2) +  
  geom_jitter(aes(x = time, y = area, fill = metabolite), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
              size = 1, alpha = 0.8, inherit.aes = FALSE) +
  labs(x = "Time, h", y = "Peak Area") +        
  scale_fill_viridis_d() +   
  scale_y_continuous(label = scientific_10) +
  theme_Publication() +   
  facet_grid(condition2 ~ condition1, scales = 'free') +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),         
        #legend.position = "right",         
        legend.direction = "vertical",         
        legend.key.size= unit(0.5, "cm"),         
        legend.key.height = unit(0.5, "cm"),         
        legend.margin = unit(1, "cm"),         
        legend.title = element_text(face="italic"))

ggsave(filename = file.path(wd_path, "Plots", "1_Fludrocortisone_acetate.pdf"),
       device   = cairo_pdf, 
       width    = 4,
       height   = 5,
       units    = "in")

#### Dexamethasone phosphate ----
DP_raw <- read_csv(paste0(data_path, "Dexamethasone-phosphate.csv")) %>% 
  mutate(step = as.factor(ifelse(nchar(condition) >= 3 & nchar(condition) <= 4, 1, 
                                 ifelse(nchar(condition) >= 7 & nchar(condition) <= 9, 2, 
                                        ifelse(nchar(condition) >= 11 & nchar(condition) <= 14, 3, 0))))) %>% 
  filter(step == 2) %>% 
  select(-c(step, compound, replicate)) %>%
  filter(time != 'filtrate') %>% 
  pivot_longer(cols = c(3:5),
               names_to = "metabolite",
               values_to = "area")

DP_raw$condition <- gsub('Csci', 'C.scindens', DP_raw$condition)
DP_raw$condition <- gsub('Pvul', 'P.vulgatus', DP_raw$condition)

DP_raw$metabolite <- factor(DP_raw$metabolite, 
                              levels =c("Dexamethasone-phosphate-473.1735", 
                                        "Dexamethasone-393.2072", 
                                        "Dexamethasone-333.1860"))

DP <- DP_raw %>% 
  group_by(condition, time, metabolite) %>%  
  mutate(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE),
            max_area  = max(area, na.rm = TRUE),
            .groups   = "drop") %>% 
  mutate(norm_area    = ifelse(max_area == 0, 0, mean_area / max_area),
         sd_norm_area = ifelse(max_area == 0, 0, sd_area / max_area)) 

DP <- DP %>% 
  separate(condition, into = c("condition1","condition2"), sep = "-")

DP$condition1 <- factor(DP$condition1, 
                          levels =c("GMM",  'C.scindens', 'P.vulgatus'))  

DP$condition2 <- factor(DP$condition2, 
                          levels =c("GMM",  'C.scindens', 'P.vulgatus')) 

DP %>% 
  ungroup() %>%    
  mutate(time = factor(time)) %>%  
  ggplot(mapping = aes(x = time, y = mean_area, fill = metabolite, group = metabolite)) +  
  geom_bar(stat = "identity", position = position_dodge(1)) +   
  geom_errorbar(aes(ymin = mean_area - sd_area, ymax = mean_area + sd_area),                  
                position = position_dodge(1), width = 0.2) +  
  geom_jitter(aes(x = time, y = area, fill = metabolite), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
              size = 1, alpha = 0.8, inherit.aes = FALSE) +
  labs(x = "Time, h", y = "Peak Area") +        
  scale_fill_viridis_d() + 
  scale_y_continuous(label = scientific_10) +
  theme_Publication() +   
  facet_grid(condition2 ~ condition1, scales = 'free') +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),         
        #legend.position = "right",         
        legend.direction = "vertical",         
        legend.key.size= unit(0.5, "cm"),         
        legend.key.height = unit(0.5, "cm"),         
        legend.margin = unit(1, "cm"),         
        legend.title = element_text(face="italic"))

ggsave(filename = file.path(wd_path, "Plots", "1_Dexamethasone_phosphate.pdf"),
       device   = cairo_pdf, 
       width    = 4,
       height   = 5,
       units    = "in")
