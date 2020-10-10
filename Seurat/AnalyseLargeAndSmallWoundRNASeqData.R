############################################################################################################
### Integration and clustering analysis of unwounded, large wound centre and periphery, and small wound
### murine skin tissue scRNA-seq from Abbasi et al. (2020) (https://doi.org/10.1016/j.stem.2020.07.008).
### There's 7 conditions that need to be integrated:
###
### Uninjured
### Day 14 post-wounding: small wound, large wound centre, large wound periphery
### Day 18 post-wounding: small wound, large wound centre, large wound peripherygoogle
############################################################################################################

### Load the relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

### Load the data

# The various experimental conditions we will consider
# conditions <- c("Uninjured", "SmallWound_PWD8", "SmallWound_PWD14", "LargeWoundCentre_PWD14", "LargeWoundPeriphery_PWD14", "LargeWoundCentre_PWD18", "LargeWoundPeriphery_PWD18")
conditions <- c("Uninjured", "SmallWound_PWD8", "SmallWound_PWD14", "LargeWoundCentre_PWD14", "LargeWoundPeriphery_PWD14")

# Read in each dataset
abbasi2020.list <- NULL
num_conditions <- length(conditions)

for (i in 1:num_conditions)
{
  dataDirectory = paste("~/Documents/CellCellCommunicationModelling/Data/Abbasi2020/", conditions[[i]], "/", sep="")
  data <- Read10X(data.dir=dataDirectory)
  
  # Initialise the Seurat object (with the raw, non-normalised data)
  object <- CreateSeuratObject(counts = data, project=conditions[[i]], min.cells = 3, min.features = 200)
  object@meta.data$sample <- conditions[[i]]
  abbasi2020.list[[i]] <- object
}

# We need to merge the list into a single Seurat object for Harmony
abbasi2020.merged <- merge(abbasi2020.list[[1]], y = c(abbasi2020.list[2:num_conditions]))

# Do the standard pre-processing
abbasi2020.merged[["percent.mt"]] <- PercentageFeatureSet(abbasi2020.merged, pattern="^mt-") # Filter out low-quality cells
VlnPlot(abbasi2020.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Plot the feature count, RNA count, and percentage of mitochondrial genes per cell

abbasi2020.merged <- subset(abbasi2020.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10) # Subset the data
abbasi2020.merged <- NormalizeData(abbasi2020.merged, verbose = FALSE)
abbasi2020.merged <- FindVariableFeatures(abbasi2020.merged, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
abbasi2020.merged <- ScaleData(abbasi2020.merged, verbose = FALSE)
abbasi2020.merged <- RunPCA(pc.genes = abbasi2020.merged@var.genes, abbasi2020.merged, npcs = 20, verbose = FALSE)

# Compare the PCAs across the samples to see how much integration we need to do
p1 <- DimPlot(object = abbasi2020.merged, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = abbasi2020.merged, features = "PC_1", group.by = "sample", pt.size = .1)
p1 + p2

abbasi2020.merged <- JackStraw(abbasi2020.merged, num.replicate = 100)
abbasi2020.merged <- ScoreJackStraw(abbasi2020.merged, dims = 1:20)

JackStrawPlot(abbasi2020.merged, dims = 1:20)
ElbowPlot(abbasi2020.merged)

# Run Harmony to project the data 
abbasi2020.merged <- RunHarmony(abbasi2020.merged, "sample", max.iter.harmony = 10, plot_convergence = TRUE)
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = abbasi2020.merged, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = abbasi2020.merged, features = "harmony_1", group.by = "sample", pt.size = .1)

# A more quantitative way of determining dimensionality, as explained explained here https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
pct <- abbasi2020.merged[["harmony"]]@stdev / sum(abbasi2020.merged[["harmony"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct) # Calculate cumulative percents for each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
num_pcs_merged <- min(co1, co2)

# Run the UMAP (all dimensionalities seem to be significant)
abbasi2020.merged <- RunUMAP(abbasi2020.merged, reduction = "harmony", dims = 1:num_pcs_merged)

# Cluster and find the neighbours
abbasi2020.merged <- FindNeighbors(abbasi2020.merged, reduction = "harmony", dims = 1:num_pcs_merged) 
abbasi2020.merged <- FindClusters(abbasi2020.merged, resolution = seq(0.3, 0.7, 0.05)) # We check which resolution is best

### We now try and figure out what the optimal resolution is (they all seem to be fine, let's go with 0.3)
Idents(abbasi2020.merged) <- "RNA_snn_res.0.3"

abbasi2020.merged.markers <- FindAllMarkers(abbasi2020.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

abbasi2020.merged.top10.markers <- abbasi2020.merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
abbasi2020.merged.top50.markers <- abbasi2020.merged.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(abbasi2020.merged, features = abbasi2020.merged.top10.markers$gene) + NoLegend()

# Let's relabel the clusters now
abbasi2020.merged.old.idents <- levels(abbasi2020.merged)
# We use the markers identified in Abbasi et al. (2020) and Guerrero-Juarez (2019) to help us identify what's what.
# Here, I will put more specifically what I think the cell types are
# 0 - Fibroblast-1 (periphery), 1 - Fibroblasts  -2 (centre), 2 - Fibroblasts-3 (periphery) 3 - Immune cell (could be macrophage or DC)
# 4 - Myofibroblast, 5 - Epidermal (has both basal and spinous), 6 - Endothelial 1, 7 - Fibroblasts 4
# 8 - Endothelial cell 2 9 - Fibroblast 5  10 - Schwann, 11 - Immune 2
abbasi2020.merged.new.idents <- c("Fibroblast 1", "Fibroblast 2", "Fibroblast 3", "Immune 1", "Myofibroblast", "Epidermal",
                                  "Endothelial 1", "Fibroblast 4", "Endothelial 2", "Fibroblast 5", "Schwann", "Immune 2")
names(abbasi2020.merged.new.idents) <- abbasi2020.merged.old.idents

# Rename the cell clusters and save this data
abbasi2020.merged <- RenameIdents(abbasi2020.merged, abbasi2020.merged.new.idents) # Rename the clusters

# Let's re-cluster different subsets to decide how certain groups can be marked
abbasi2020.fibroblasts <- subset(abbasi2020.merged, idents = c("Fibroblast 1", "Fibroblast 2",  "Fibroblast 3", "Fibroblast 4", "Fibroblast 5", "Myofibroblast"))
abbasi2020.immune <- subset(abbasi2020.merged, idents = c("Immune 1", "Immune 2"))
abbasi2020.epidermal <- subset(abbasi2020.merged, idents = c("Epidermal"))
abbasi2020.endothelial <- subset(abbasi2020.merged, idents = c("Endothelial 1", "Endothelial 2"))

### Let's recluster these subsets, as we know there should be better cell types
### First the epidermal, as we know there should be both basal and suprabasal cells
abbasi2020.epidermal <- RunPCA(abbasi2020.epidermal, npcs = 20, verbose = FALSE)
abbasi2020.epidermal <- RunHarmony(abbasi2020.epidermal, "sample", plot_convergence = TRUE)

pct <- abbasi2020.epidermal[["harmony"]]@stdev / sum(abbasi2020.epidermal[["harmony"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct) # Calculate cumulative percents for each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
num_pcs_epidermal <- min(co1, co2)

abbasi2020.epidermal <- RunUMAP(abbasi2020.epidermal, reduction = "harmony", dims = 1:num_pcs_epidermal)
abbasi2020.epidermal <- FindNeighbors(abbasi2020.epidermal, reduction = "harmony", dims = 1:num_pcs_merged) 
abbasi2020.epidermal <- FindClusters(abbasi2020.epidermal, resolution = seq(0.2, 0.5, 0.05)) # We check which resolution is best
Idents(abbasi2020.epidermal) <- "RNA_snn_res.0.4" # 0.4 seems to produce the best resolution (we were having trouble splitting certain basal from suprabasal cells)

abbasi2020.epidermal.markers <- FindAllMarkers(abbasi2020.epidermal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

abbasi2020.epidermal.top10.markers <- abbasi2020.epidermal.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
abbasi2020.epidermal.top50.markers <- abbasi2020.epidermal.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(abbasi2020.epidermal, features = abbasi2020.epidermal.top10.markers$gene) + NoLegend()

# We relabel the epidermal labels now
abbasi2020.epidermal.old.idents <- levels(abbasi2020.epidermal)
# We use the markers identified in Abbasi et al. (2020), Haensel (2020), and Guerrero-Juarez (2019) to help us identify what's what.
# Here, I will put more specifically what I think the cell types are
# 0 - Epidermal basal 1, 1 - Epidermal basal 2, 2 - Epidermal spinous 1, 3 - Epidermal spinous 2, 
# 4 - Epidermal basal 3, 5 - Epidermal basal 5
abbasi2020.epidermal.new.idents <- c("Epidermal basal 1", "Epidermal basal 2", "Epidermal spinous 1", "Epidermal spinous 2", "Epidermal basal 3", "Epidermal basal 4")
names(abbasi2020.epidermal.new.idents) <- abbasi2020.epidermal.old.idents

abbasi2020.epidermal <- RenameIdents(abbasi2020.epidermal, abbasi2020.epidermal.new.idents) # Rename the epidermal clusters

### Now cluster the immune cells
abbasi2020.immune <- RunPCA(abbasi2020.immune, npcs = 20, verbose = FALSE)
abbasi2020.immune <- RunHarmony(abbasi2020.immune, "sample", plot_convergence = TRUE)

pct <- abbasi2020.immune[["harmony"]]@stdev / sum(abbasi2020.immune[["harmony"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct) # Calculate cumulative percents for each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
num_pcs_immune <- min(co1, co2)

abbasi2020.immune <- RunUMAP(abbasi2020.immune, reduction = "harmony", dims = 1:num_pcs_immune)
abbasi2020.immune <- FindNeighbors(abbasi2020.immune, reduction = "harmony", dims = 1:num_pcs_immune) 
abbasi2020.immune <- FindClusters(abbasi2020.immune, resolution = seq(0.2, 0.5, 0.05)) # We check which resolution is best

Idents(abbasi2020.immune) <- "RNA_snn_res.0.2" # 0.2 seems to produce the best resolution

abbasi2020.immune.markers <- FindAllMarkers(abbasi2020.immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

abbasi2020.immune.top10.markers <- abbasi2020.immune.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
abbasi2020.immune.top50.markers <- abbasi2020.immune.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(abbasi2020.immune, features = abbasi2020.immune.top10.markers$gene) + NoLegend()

# We relabel the immune labels now
abbasi2020.immune.old.idents <- levels(abbasi2020.immune)
# We use the markers identified in Abbasi et al. (2020), Haensel (2020), and Guerrero-Juarez (2019) to help us identify what's what.
# Here, I will put more specifically what I think the cell types are
# 0 -  Macrophage 1 1 - Macrophage 2 2 - Dendritic cells
# 3 - Macrophage 3 4 - T cell 1 5 - T cell 2
abbasi2020.immune.new.idents <- c("Macrophage 1", "Macrophage 2","Dendritic cell", "Macrophage 3", "T cell", "Other immune cell")
names(abbasi2020.immune.new.idents) <- abbasi2020.immune.old.idents

abbasi2020.immune <- RenameIdents(abbasi2020.immune, abbasi2020.immune.new.idents) # Rename the immune clusters

### Look at the reclustering the fibroblasts
abbasi2020.fibroblasts <- RunPCA(abbasi2020.fibroblasts, npcs = 20, verbose = FALSE)
abbasi2020.fibroblasts <- RunHarmony(abbasi2020.fibroblasts, "sample", plot_convergence = TRUE)

pct <- abbasi2020.fibroblasts[["harmony"]]@stdev / sum(abbasi2020.fibroblasts[["harmony"]]@stdev) * 100 # Determine percent of variation associated with each PC
cumu <- cumsum(pct) # Calculate cumulative percents for each PC
co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
num_pcs_fibroblasts <- min(co1, co2)

abbasi2020.fibroblasts <- RunUMAP(abbasi2020.fibroblasts, reduction = "harmony", dims = 1:num_pcs_fibroblasts)
abbasi2020.fibroblasts <- FindNeighbors(abbasi2020.fibroblasts, reduction = "harmony", dims = 1:num_pcs_fibroblasts) 
abbasi2020.fibroblasts <- FindClusters(abbasi2020.fibroblasts, resolution = seq(0.2, 0.5, 0.05)) # We check which resolution is best

Idents(abbasi2020.fibroblasts) <- "RNA_snn_res.0.2" # 0.2 seems to produce the best resolution

abbasi2020.fibroblasts.markers <- FindAllMarkers(abbasi2020.fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

abbasi2020.fibroblasts.top10.markers <- abbasi2020.fibroblasts.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
abbasi2020.fibroblasts.top50.markers <- abbasi2020.fibroblasts.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(abbasi2020.fibroblasts, features = abbasi2020.fibroblasts.top10.markers$gene) + NoLegend()


# Check which clusters line up with Abbasi (2020)
# featurePlot <- FeaturePlot(abbasi2020.merged, features = c("Dlk1", "Sca1", "Mest", "Crabp1", "Fabp5", "Runx1"))
swd8Plot <- FeaturePlot(abbasi2020.merged, features = c("Saa3", "Tnc", "Gm26917")) # Small wound day 8
swd14Plot <- FeaturePlot(abbasi2020.merged, features = c("Mgp", "Dlk1", "Fbln7", "Thbs4")) # Small wound day 14
lwu14Plot <- FeaturePlot(abbasi2020.merged, features = c("Fabp5", "Crabp1", "Prss35")) # Large wound day 14 upper fibroblast markers
lwd14Plot <- FeaturePlot(abbasi2020.merged, features = c("Dlk1", "Sca1", "Mest")) # Large wound day 14 upper fibroblast markers
uwPlot <- FeaturePlot(abbasi2020.merged, features = c("Tnmd", "Ptn", "Igfbp2", "Crabp1")) # Uninjured
umapPlot <- DimPlot(abbasi2020.merged, reduction = "umap", group.by = "sub.cluster", label = TRUE)
### From this, Fibroblast 1, Fibroblast 3, Fibroblasts 4 line up with "lower dermis fibroblasts"

# Let's add the sub-cluster labels
abbasi2020.merged$sub.cluster <- as.character(Idents(abbasi2020.merged))
abbasi2020.merged$sub.cluster[Cells(abbasi2020.epidermal)] <- paste("", Idents(abbasi2020.epidermal), sep="") # Rename the epidermal labels
abbasi2020.merged$sub.cluster[Cells(abbasi2020.immune)] <- paste("", Idents(abbasi2020.immune), sep="") # Rename the immune labels

DimPlot(abbasi2020.merged, group.by = "sample", label = FALSE)

save(abbasi2020.merged, abbasi2020.epidermal, abbasi2020.immune, abbasi2020.merged.markers, abbasi2020.epidermal.markers, abbasi2020.immune.markers, abbasi2020.merged.new.idents, abbasi2020.epidermal.new.idents, abbasi2020.immune.new.idents, file = "~/Documents/CellCellCommunicationModelling/Data/Abbasi2020/abbasi2020integrated.rdata")
