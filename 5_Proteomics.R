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
            legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#### SET WD ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# Threshold for log2FC
for (thr in c(1,2)) {
  
  ### READ AND PLOT DATA
  NB <- read_csv(paste0(wd_path, '/Input/Proteomics_NB.csv')) %>% 
    mutate(logFC_adj = ifelse(grepl("^ctrl-", comparison), -logFC, logFC),
           hit_flag = abs(logFC_adj) > thr,
           negLogP = -log10(adj.P.Val))
  
  # Loop over each comparison and save a PDF
  for (c in unique(NB$comparison)) {
    df <- NB %>% filter(comparison == c)
    
    xlim <- range(df$logFC_adj) * 1.05
    ylim <- range(df$negLogP) * 1.05
    
    p <- ggplot(df, aes(x = logFC_adj, y = negLogP)) +
      
      geom_point(color = "lightgrey", size = 2) +
      geom_point(data = df %>% filter(hit_flag), color = "red3", size = 2.5) +
      geom_text_repel(data = df %>% filter(hit_flag),
                      aes(label = Gene), size = 3, max.overlaps = 50) +
      
      coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
      
      # Top border
      annotate("segment", x = xlim[1], xend = xlim[2],
               y = ylim[2], yend = ylim[2], size = 1) +
      
      # Right border
      annotate("segment", x = xlim[2], xend = xlim[2],
               y = ylim[1], yend = ylim[2], size = 1) +
      
      labs(title = paste('NB:', c), x = "log2 Fold Change", y = "-log10(adj. p-value)") +
      theme_Publication()
    
    ggsave(
      filename = paste0(wd_path, '/Plots/Proteomics/NB_volcano_thr_',thr, '_', c, ".pdf"),
      plot = p,
      width = 3,
      height = 3,
      device = cairo_pdf
    )
  }
}

NB_hits <- NB %>%
  filter(adj.P.Val <= 0.01,
         abs(logFC_adj) >= 1) %>%
  select(comparison, Entry.Name, Gene, Organism, Protein.Description) %>%
  group_by(Entry.Name, Gene, Organism, Protein.Description) %>%
  mutate(sample_n = n()) %>%
  ungroup()

write_csv(NB_hits, paste0(wd_path, "/Output/5_Proteomics_NB_hits.csv"))

