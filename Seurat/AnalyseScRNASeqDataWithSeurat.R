### Our attempt to analyse various datasets of single-cell RNA-seq datasets that focus on wound healing and also sequence fibroblasts.

# Load the relevant libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)

#####################################################################################################################
### We first consider the dataset by Haensel et al. (2020), who analysed murine skin stem cell dynamics in healthy and
### wounded skin.
#########################################################################################################

# Load the datasets. We lost the unwounded and wounded data into separate lists first.

### Let's first analyse the unwounded data
unwounded.list <- NULL
for (i in 1:2)
{
  dataDirectory = paste("~/Documents/CellCellCommunicationModelling/Data/Haensel2020/Un-Wounded_", i, sep="")
  unwounded.data <- Read10X(data.dir=dataDirectory)
  
  # Initialise the Seurat object (with the raw, non-normalised data)
  unwounded.list[[i]] <- CreateSeuratObject(counts = unwounded.data, project=paste("UW-", i, sep=""), min.cells = 3, min.features = 200)
}

# Normalise and find the variable featuers in each dataset
for (i in 1:length(unwounded.list))
{
  unwounded.list[[i]] <- NormalizeData(unwounded.list[[i]], verbose = FALSE)
  unwounded.list[[i]] <- FindVariableFeatures(unwounded.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# Find the integration anchors
unwounded.anchors <- FindIntegrationAnchors(object.list = unwounded.list, dims = 1:30)

# Integrate the data
unwounded.integrated <- IntegrateData(anchorset = unwounded.anchors, dims = 1:30)

# Switch to integrated array
DefaultAssay(unwounded.integrated) <- "integrated"

# Now scale and run PCA
# unwounded.integrated <- RenameCells(unwounded.integrated) # Just to check 
unwounded.integrated <- ScaleData(unwounded.integrated, verbose = FALSE)
unwounded.integrated <- RunPCA(unwounded.integrated, npcs = 30, verbose = FALSE)

# Try and determine the dimensionality of the dataset
# unwounded.integrated <- JackStraw(unwounded.integrated, num.replicate = 100)
# unwounded.integrated <- ScoreJackStraw(unwounded.integrated, dims = 1:20)


# Run the UMAP (all dimensionalities seem to be significant)
unwounded.integrated <- RunTSNE(unwounded.integrated, reduction = "pca", dims = 1:10)
unwounded.integrated <- RunUMAP(unwounded.integrated, reduction = "pca", dims = 1:10)

# Cluster and find the neighbours
unwounded.integrated <- FindNeighbors(unwounded.integrated, reduction = "pca", dims = 1:10)
unwounded.integrated <- FindClusters(unwounded.integrated, resolution = 0.6)

unwounded.markers <- FindAllMarkers(unwounded.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

# Plot the UMAP/tSNE
tsneuw <- DimPlot(unwounded.integrated, reduction = "tsne")
umapuw <- DimPlot(unwounded.integrated, reduction = "umap", label = TRUE)
# tsneuw + umapuw

# We've found parameters that seem to reproduce the clusters in the paper. 

### Run the same thing for the wounded datasets
wounded.list <- NULL
for (i in 1:3)
{
  dataDirectory = paste("~/Documents/CellCellCommunicationModelling/Data/Haensel2020/Wounded_", (i), sep="")
  wounded.data <- Read10X(data.dir=dataDirectory)
  
  # Initialise the Seurat object (with the raw, non-normalised data)
  wounded.list[[i]] <- CreateSeuratObject(counts = wounded.data, project=paste("WO-", i, sep=""), min.cells = 3, min.features = 200)
}

# Normalise and find the variable featuers in each dataset
for (i in 1:length(wounded.list))
{
  wounded.list[[i]] <- NormalizeData(wounded.list[[i]], verbose = FALSE)
  wounded.list[[i]] <- FindVariableFeatures(wounded.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# Find the integration anchors
wounded.anchors <- FindIntegrationAnchors(object.list = wounded.list, dims = 1:30)

# Integrate the data
wounded.integrated <- IntegrateData(anchorset = wounded.anchors, dims = 1:30)

# Switch to integrated array
DefaultAssay(wounded.integrated) <- "integrated"

# Now scale and run PCA
wounded.integrated <- ScaleData(wounded.integrated, verbose = FALSE)
wounded.integrated <- RunPCA(wounded.integrated, npcs = 30, verbose = FALSE)
# 
# # Try and determine the dimensionality of the dataset
# wounded.integrated <- JackStraw(wounded.integrated, num.replicate = 100)
# wounded.integrated <- ScoreJackStraw(wounded.integrated, dims = 1:20)

# Run the UMAP (all dimensionalities seem to be significant)
# wounded.integrated <- RunTSNE(wounded.integrated, reduction = "pca", dims = 1:10)
wounded.integrated <- RunUMAP(wounded.integrated, reduction = "pca", dims = 1:10)

# Cluster and find the neighbours
wounded.integrated <- FindNeighbors(wounded.integrated, reduction = "pca", dims = 1:10)
wounded.integrated <- FindClusters(wounded.integrated, resolution = 0.6)

wounded.markers <- FindAllMarkers(wounded.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # Identify gene markers for each cluster

# Plot the dimensional reductions
# tsnewo <- DimPlot(wounded.integrated, reduction = "tsne")
umapwo <- DimPlot(wounded.integrated, reduction = "umap", label = TRUE)
# tsnewo + umapwo

##########################################################################################################################################
### Now analyse the dataset of Guerrero-Juarez et al. (2019), which looks at the fibroblast population during wound healing and scarring. 
### Here, it'll be a bit different, because we'll need to read in the raw data before converting that into a Seurat object
##########################################################################################################################################
dataDirectory = "~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/Shortened/"

woundedgj.data <- Read10X(data.dir=dataDirectory)

woundedgj <- CreateSeuratObject(counts = woundedgj.data, project="GJ2019", min.cells = 3, min.features = 200)

# Check for low-quality cells by counting the percentage of mitochondrial gene expression
woundedgj[["percent.mt"]] <- PercentageFeatureSet(woundedgj, pattern="^MT-")
woundedgj <- subset(woundedgj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 8) # Subset the data

# Normalise and scale the data
woundedgj <- NormalizeData(woundedgj, verbose = FALSE)
woundedgj <- FindVariableFeatures(woundedgj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

woundedgj <- ScaleData(woundedgj, verbose = FALSE)
woundedgj <- RunPCA(woundedgj, npcs = 40, verbose = FALSE)

# Cluster and find the neighbours
woundedgj<- FindNeighbors(woundedgj, reduction = "pca", dims = 1:35)
woundedgj <- FindClusters(woundedgj, resolution = 0.45)

woundedgj <- RunTSNE(woundedgj, reduction = "pca", dims = 1:35)
woundedgj <- RunUMAP(woundedgj, reduction = "pca", dims = 1:35)

umapgj <- DimPlot(woundedgj, reduction = "umap", label = TRUE)
tsnegj <- DimPlot(woundedgj, reduction = "tsne")

woundedgj.markers <- FindAllMarkers(woundedgj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

woundedgj.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

##########################################################################################################################################
### Suoqin has given us the .RData files for the Suerat objects that he used to produce the figures in Guerrero-Juarez et al. (2019) and
### Haensel et al. (2020). By his recommendation, it should be okay to subset out the fibroblasts and just subset 
##########################################################################################################################################

unwoundedHaensel <- load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/bs_cca.rdata")
woundedHaensel <- load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/sw_cca.rdata")
woundedPWD12 <- load("~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/10X_wounded_PWD12.rdata")
