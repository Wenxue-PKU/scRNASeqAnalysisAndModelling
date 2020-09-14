### Clustering of the human foreskin epidermis data analysed in Wang et al. (2020). Namely, we are interested in examining the data
### for any indication of differential adhesion that may be revealed by gene expression amongst cell types.

# Load the relevant libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

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
  epidermis.list[[i]] <- CreateSeuratObject(counts = epidermis.data, project=paste("Donor_", i, sep=""), min.cells = 3, min.features = 200)
}

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

num_pcs_epidermis <- min(co1, co2)

# Run the UMAP (all dimensionalities seem to be significant)
epidermis.integrated <- RunUMAP(epidermis.integrated, reduction = "pca", dims = 1:num_pcs_epidermis)

# Cluster and find the neighbours
epidermis.integrated <- FindNeighbors(epidermis.integrated, reduction = "pca", dims = 1:num_pcs_epidermis) 
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

