rm(list = ls())

#### Import libraries ----
library(ChemmineR)
library(tidyverse)
library(ggdendro)

#### READ AND PREPARE DATA ----

# Set paths to working directory and data folders
wd_path <- getwd()
setwd(wd_path)

# Import SMILES
smiset <- read.SMIset(paste0(wd_path, "/Input/SMILES.smi"))

# Convert smiles into SDF 
sdfset <- smiles2sdf(smiset)

# Generate fingerprints from Open Babel
fpset <- fingerprintOB(sdfset, "FP2") ## FP2 indexes linear fragments up to 7 atoms.

# Create Atom Pair fingerprint
apset <- sdf2ap(sdfset)

# Clustering
clusters <- cmp.cluster(db = apset, cutoff = c(0.7, 0.8, 0.9), quiet = TRUE, save.distances = "distmat.rda")

#### GENERATE DENDROGRAMS ----

# Hierarchical clustering with hclust
load("distmat.rda")
hc <- hclust(as.dist(distmat), method = "single") 
hc[["labels"]] <- cid(apset) # Assign correct item labels 

# Create dendrograms
dendrogram <- ggdendrogram(hc, rotate = TRUE, size = 2) 

ggsave(filename = paste0(wd_path, "/Plots/Dendrogram_complete.pdf"),
       plot = dendrogram,
       bg = "transparent",
       width = 5,
       height = 5,
       unit = "in",
       dpi = 300)

dhc <- as.dendrogram(hc)
data <- dendro_data(dhc, type = "rectangle")

dendrogram_plain <- ggplot(segment(data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 2) + 
  theme_dendro() +
  coord_flip() + 
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA)) 

ggsave(filename = paste0(wd_path, "/Plots/Dendrogram_without_labels.pdf"),
       plot = dendrogram_plain,
       bg = "transparent",
       width = 5,
       height = 30,
       unit = "in",
       dpi = 300)
