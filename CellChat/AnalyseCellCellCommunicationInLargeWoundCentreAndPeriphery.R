### Analysis of the differences in signalling pathways between fibrotic and regenerative scars, as analysed in Gay et al. (2020).
### We first reproduce the clustering analysis reported in the paper, then 

# Load the relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)

### We first need to cluster the scRNA-seq data, as CellChat requires these labels.
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact