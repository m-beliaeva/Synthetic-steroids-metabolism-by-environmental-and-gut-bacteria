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

# Set paths to working directory
wd_path <- getwd()
setwd(wd_path)
data_path <- paste0(wd_path, "/Input/")

#### DATA ----
data <- read_csv(paste0(data_path, "Gut_community_Sherbicidovorans.csv")) %>% 
  separate(condition, into = c("condition1", "condition2"), sep = "-") %>% 
  mutate(condition2 = gsub("Sh", "S.herbicidovorans", condition2)) %>% 
  filter(!position%in%c("B5", "B6", "B7", "B8")) # Keep 4 control replicates

data$condition1 <- factor(data$condition1, 
                          levels = c( "GMM", "Comm1","Comm2","Comm3","Comm4"))  
                          
data$condition2 <- factor(data$condition2, 
                          levels =c( "NB", "S.herbicidovorans")) 

#### Fludrocortisone-acetate ----
FA <- data %>% 
  filter(compound == 'Fludrocortisone-acetate') %>% 
  select(-position) %>% 
  group_by(condition1, condition2, time, metabolite) %>%  
  mutate(mean_area = mean(area, na.rm = TRUE),
            sd_area   = sd(area, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(metabolite != "Fludrocortisone_321.1861_2.4") %>%
  filter(metabolite != "Fludrocortisone_381.2072_2.4") %>% 
  #filter(metabolite != "Fludrocortisone_383.2228_2.4") %>% 
  filter(metabolite != "Fludrocortisone_323.2017_2.2") #%>% 
  #filter(condition1 != "GMM")

FA$metabolite <- gsub("Fludrocortisone_17-Ketone", "Fludrocortisone-17-Ketone", FA$metabolite)
FA$metabolite <- gsub("Fludrocortisone_319.1704_2.9", "1,2-Dehydrofludrocortisone-17-Ketone", FA$metabolite)
FA$metabolite <- gsub("Fludrocortisone_383.2228_2.4", "Dihydro-fludrocortisone", FA$metabolite)
# FA$metabolite <- gsub("Fludrocortisone_323.2017_2.2", "Dihydro-fludrocortisone-17-Ketone", FA$metabolite)

FA$metabolite <- factor(FA$metabolite,
                        levels =c("Fludrocortisone-acetate",
                                  "Fludrocortisone",
                                  "Dihydro-fludrocortisone",
                                  "Fludrocortisone-17-Ketone", 
                                  "1,2-Dehydrofludrocortisone-17-Ketone"))
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

ggsave(filename = file.path(wd_path, "Plots", "3_Fludrocortisone_acetate.pdf"),
       device   = cairo_pdf, 
       width    = 12,
       height   = 6,
       units    = "in")
