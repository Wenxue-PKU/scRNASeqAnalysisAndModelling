### Clustering of the human foreskin epidermis data analysed in Wang et al. (2020). Namely, we are interested in examining the data
### for any indication of differential adhesion that may be revealed by gene expression amongst cell types.

# Load the relevant libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

### We will need to integrate the five samples. Here, we will just use Seurat's standard workflow for integration, 
### using the parameter values specified in Wang et al. (2020).

# Read in each dataset
epidermis.list <- NULL
num_donors <- 5
for (i in 1:num_donors)
{
  dataDirectory = paste("~/Documents/CellCellCommunicationModelling/Data/Wang2020/Donor_", i, sep="")
  epidermis.data <- Read10X(data.dir=dataDirectory)
  
  # Initialise the Seurat object (with the raw, non-normalised data)
  epidermis.object <- CreateSeuratObject(counts = epidermis.data, project=paste("Donor_", i, sep=""), min.cells = 3, min.features = 200)
  epidermis.object@meta.data$sample <- paste("Sample_", i, sep="")
  epidermis.list[[i]] <- epidermis.object
}

# We need to merge the list into a single Seurat object for Harmony
epidermis.merged <- merge(epidermis.list[[1]], y = c(epidermis.list[2:5]))

# Filter out the data for over-expressed mitochondrial genes, which indicate dying/low-quality cells
for (i in 1:num_donors)
{
  epidermis.object <- epidermis.list[[i]]
  epidermis.object[["percent.mt"]] <-PercentageFeatureSet(epidermis.object, pattern="^MT-")
  
  # Subset each data
  epidermis.list[[i]] <- subset(epidermis.object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) # Subset the data
}

# Normalise and find the variable featuers in each dataset
for (i in 1:num_donors)
{
  epidermis.list[[i]] <- NormalizeData(epidermis.list[[i]], verbose = FALSE)
  epidermis.list[[i]] <- FindVariableFeatures(epidermis.list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}


######################################################################################################
### FIRST run the integration stuff on the data that has been integrated using Seurat's method
# Find the integration anchors
epidermis.anchors <- FindIntegrationAnchors(object.list = epidermis.list, dims = 1:30)

# Integrate the data
epidermis.integrated <- IntegrateData(anchorset = epidermis.anchors, dims = 1:30)

# Switch to integrated array
DefaultAssay(epidermis.integrated) <- "integrated"

# Now scale and run PCA
epidermis.integrated <- ScaleData(epidermis.integrated, verbose = FALSE)
epidermis.integrated <- RunPCA(epidermis.integrated, npcs = 30, verbose = FALSE)

# A more quantitative way of determining dimensionality, as explained explained here https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# Determine percent of variation associated with each PC
pct <- epidermis.integrated[["pca"]]@stdev / sum(epidermis.integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_epidermis_integrated <- min(co1, co2)

# Run the UMAP (all dimensionalities seem to be significant)
epidermis.integrated <- RunUMAP(epidermis.integrated, reduction = "pca", dims = 1:num_pcs_epidermis_integrated)

# Cluster and find the neighbours
epidermis.integrated <- FindNeighbors(epidermis.integrated, reduction = "pca", dims = 1:num_pcs_epidermis_integrated) 
epidermis.integrated <- FindClusters(epidermis.integrated, resolution = 0.45) # We check which resolution is best

Idents(epidermis.integrated) <- "integrated_snn_res.0.45" # 0.45 seems to be the best?

epidermis.markers <- FindAllMarkers(epidermis.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

# Get the top marker from each cluster to see how well these distinguish teh clusters
epidermis.top1.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) # Get the top 5 markers for each cluster

FeaturePlot(epidermis.integrated, features = epidermis.top1.markers$gene) # Plot the top marker gene for each cluster

epidermis.top5.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

epidermis.top10.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
epidermis.top50.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(epidermis.integrated, features = epidermis.top10.markers$gene) + NoLegend()

save(epidermis.integrated, epidermis.markers, file = "~/Documents/CellCellCommunicationModelling/Data/Wang2020/epidermiswithseuratintegration.rdata")

###########################################################################################################################
############ Now let's integrate using Harmony's methods 

# First run the standard pre-processing
epidermis.merged[["percent.mt"]] <- PercentageFeatureSet(epidermis.merged, pattern="^mt-") # Filter out low-quality cells
epidermis.merged <- subset(epidermis.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) # Subset the data
epidermis.merged <- NormalizeData(epidermis.merged, verbose = FALSE)
epidermis.merged <- FindVariableFeatures(epidermis.merged, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
epidermis.merged <- ScaleData(epidermis.merged, verbose = FALSE)
epidermis.merged <- RunPCA(epidermis.merged, npcs = 30)

# We show why we may need to integrate the data better
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = epidermis.merged, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = epidermis.merged, features = "PC_1", group.by = "sample", pt.size = .1)
p1 + p2

options(repr.plot.height = 2.5, repr.plot.width = 6)
epidermis.merged <- epidermis.merged %>% RunHarmony("sample", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(epidermis.merged, 'harmony')
harmony_embeddings[1:5, 1:5]

# We now plot these again to show how well Harmony did
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = epidermis.merged, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = epidermis.merged, features = "harmony_1", group.by = "sample", pt.size = .1)

# A more quantitative way of determining dimensionality, as explained explained here https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# Determine percent of variation associated with each PC
pct <- epidermis.merged[["harmony"]]@stdev / sum(epidermis.merged[["harmony"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_epidermis_harmony <- min(co1, co2)

# Run the UMAP (all dimensionalities seem to be significant)
epidermis.merged <- RunUMAP(epidermis.merged, reduction = "harmony", dims = 1:num_pcs_epidermis_harmony)

# Cluster and find the neighbours
epidermis.merged <- FindNeighbors(epidermis.merged, reduction = "harmony", dims = 1:num_pcs_epidermis_harmony) 
epidermis.merged <- FindClusters(epidermis.merged, resolution = seq(0.55, 0.65, 0.05)) # We check which resolution is best

# Let's look at which resolution works best
Idents(epidermis.merged) <- "RNA_snn_res.0.45"

epidermis.markers <- FindAllMarkers(epidermis.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

epidermis.top10.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
epidermis.top50.markers <- epidermis.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(epidermis.merged, features = epidermis.top10.markers$gene) + NoLegend()

# Check the feature plots of known basal markers
FeaturePlot(epidermis.merged, features = c("DSG1", "DSG3", "CDH1", "CDH3")) 

# Rename the epidermis clusters. We rename them according to also what was assigned
# in Wang et al. (2020)
# 0 - Spinous cells I 1 - Granular cells I 2 - Basal cells I (BIII) 3 - Granular cells II 4 - Spinous cells II
# 5 - Basal cells II (BII) 6 - Spinous cells III 7 - Spinous cells IV 8 - Spinous cells V
# 9 - Basal cells III (BI?) 10 - Basal cells IV (BIV) 11 - Spinous cells VI 12 - Dendritic cells
# 13 - 14 - 15 - 

save(epidermis.merged, epidermis.markers, file = "~/Documents/CellCellCommunicationModelling/Data/Wang2020/epidermiswithharmony.rdata")
