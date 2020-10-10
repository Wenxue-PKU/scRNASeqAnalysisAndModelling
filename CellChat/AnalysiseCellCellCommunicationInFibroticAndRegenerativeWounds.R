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

dataDirectory = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/Fibrotic/"

fibrotic.data <- Read10X(data.dir=dataDirectory)

fibrotic <- CreateSeuratObject(counts = fibrotic.data, project="Fibrotic", min.cells = 3, min.features = 200)

# Check for low-quality cells by counting the percentage of mitochondrial gene expression
fibrotic[["percent.mt"]] <- PercentageFeatureSet(fibrotic, pattern="^mt-")

# Visualise the MT-gene count to check
VlnPlot(fibrotic, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

fibrotic <- subset(fibrotic, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # Subset the data

# Normalise and scale the data
fibrotic <- NormalizeData(fibrotic, verbose = FALSE)
fibrotic <- FindVariableFeatures(fibrotic, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

fibrotic <- ScaleData(fibrotic, verbose = FALSE)

fibrotic <- RunPCA(fibrotic, npcs = 40, verbose = FALSE) # Run the PCA

# We need to calculate hte significant PCs to determine the dimensionality of the data
fibrotic <- JackStraw(fibrotic, num.replicate = 100)
fibrotic <- ScoreJackStraw(fibrotic, dims = 1:40)

# Visualise the results to figure out how many PCs we should consier.
JackStrawPlot(fibrotic, dims = 1:40)
ElbowPlot(fibrotic, ndims = 40)

# A more quantitative way of determining dimensionality, as explained explained here https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# Determine percent of variation associated with each PC
pct <- fibrotic[["pca"]]@stdev / sum(fibrotic[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_fibrotic <- min(co1, co2)

# FInd the neighbours and cluster the data. We may play with the numbers here a bit to see what works well.
fibrotic<- FindNeighbors(fibrotic, reduction = "pca", dims = 1:num_pcs_fibrotic)

# We're going to calculate the clusters at varying resolutions (think 0.45 was the winner)
fibrotic <- FindClusters(fibrotic, resolution = 0.45)

# Run the UMAp, which doesn't rely on cluster resolution
fibrotic <- RunUMAP(fibrotic, reduction = "pca", dims = 1:num_pcs_fibrotic)

# We will inspect these and decide what's the best cluster. The way to do this is to look at each
# clustering for the resolutions and look at how well the marker genes separate the clusters
Idents(object = fibrotic) <- "RNA_snn_res.0.45" # 0.45 may be the winner right now
DimPlot(fibrotic, reduction="umap", label=TRUE) 

# Let's find the markers for this resolution
fibrotic.markers <- FindAllMarkers(fibrotic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top marker from each cluster to see how well these distinguish teh clusters
fibrotic.top1.markers <- fibrotic.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) # Get the top 5 markers for each cluster

FeaturePlot(fibrotic, features = c("Acta2")) # Plot the top marker gene for each cluster

fibrotic.top5.markers <- fibrotic.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
fibrotic.top10.markers <- fibrotic.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
fibrotic.top25.markers <- fibrotic.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)
fibrotic.top50.markers <- fibrotic.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(fibrotic, features = fibrotic.top10.markers$gene) + NoLegend()

# Rename the new fibrotic identities. We may need to subcluster some of these populations to understand what's going on.
# 0 - T cells I, 1 - Macrophage I, 2 - Myofibroblast, 3 - Fibroblasts I, 4 - Dendritic cell I, 5- NK cell/T cell II, 6 - Schwann cell
# 7 - T cells III, 8 - Fibroblast II, 9 - Macrophage II, 10 - Mast cell, 11 - Dendritic cell II, 12 - Mast cell II
# 13 - Endothelial cell 14 - Fibroblast III
fibrotic.orig.idents <- fibrotic.old.idents
fibrotic.old.idents <- levels(fibrotic) # Just in case
fibrotic.new.idents <-c("T cell I", "Macrophage I", "Myofibroblast", "Fibroblast I",
                        "Dendritic cell I", "T cell II", "Schwann cell", "T cell III",
                        "Fibroblast II", "Macrophage II", "Mast cell I", "Dendritic cell II", 
                        "Mast cell II", "Lymphatic endothelial cell", "Fibroblast III")
names(fibrotic.new.idents) <- fibrotic.old.idents

fibrotic <- RenameIdents(fibrotic, fibrotic.new.idents)

save(fibrotic, fibrotic.new.idents, fibrotic.top10.markers, fibrotic.markers, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibrotic.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibrotic.rdata")

# We're going to plot some things for the preliminary analysis
fibrotic_cluster_counts <- table(Idents(fibrotic))
num_fibrotic_cells <- sum(table(Idents(fibrotic)))
fibrotic_cluster_proportions <- 100 * fibrotic_cluster_counts / num_fibrotic_cells # Calculates the percentage of the total population that each cluster is

# Plot the umap and save it for the slides
DimPlot(fibrotic, reduction="umap") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticumap.eps', device='eps')

### Now we're in a position to run CellChat on the fibrotic dataset
fibrotic.labels <- Idents(fibrotic)
fibrotic.identity <- data.frame(group = fibrotic.labels, row.names = names(fibrotic.labels)) # Dataframe of the cell labels
fibrotic.data.input <- GetAssayData(fibrotic, assay = "RNA", slot = "data")  # Normalised data matrix

fibrotic.cc <- createCellChat(data = fibrotic.data.input, do.sparse = F)

# Add the meta-data from Seurat
fibrotic.cc <- addMeta(fibrotic.cc, meta = fibrotic.identity, meta.name = "labels") 
fibrotic.cc <- setIdent(fibrotic.cc, ident.use = "labels") # Set the labels to be the default cell identity

fibrotic.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
fibrotic.cc <- subsetData(fibrotic.cc) # We subset the expression data of signalling genes to save on computational cost
fibrotic.cc <- identifyOverExpressedGenes(fibrotic.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
fibrotic.cc <- identifyOverExpressedInteractions(fibrotic.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# fibrotic.cc <- projectData(fibrotic.cc, PPI.mouse) # Other option includes PPI.human, apparently we shouldn't do this any more

# We now infer the cell-cell communication network by calculating the communication probabilities
fibrotic.cc <- computeCommunProb(fibrotic.cc) 
fibrotic.cc <- computeCommunProbPathway(fibrotic.cc) # Calculate the probabilities at the signalling level
fibrotic.cc <- aggregateNet(fibrotic.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(3, 4, 9, 15) # Focus on the fibroblasts
pathways.show <- "MK" # All possible pathways are stored in @netP$pathways
fibroticGroupSize <- as.numeric(table(fibrotic.cc@idents)) # Get the number of cells in each group
netVisual_aggregate(fibrotic.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, vertex.size = fibroticGroupSize, pt.title = 16)
netAnalysis_contribution(fibrotic.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
fibrotic.cc <- netAnalysis_signalingRole(fibrotic.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(fibrotic.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 4
fibrotic.cc <- identifyCommunicationPatterns(fibrotic.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(fibrotic.cc, pattern = "outgoing", font.size = 3, font.size.title = 12) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(fibrotic.cc, pattern = "outgoing", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 4
fibrotic.cc <- identifyCommunicationPatterns(fibrotic.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(fibrotic.cc, pattern = "incoming", font.size = 3)

netAnalysis_dot(fibrotic.cc, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
fibrotic.cc <- computeNetSimilarity(fibrotic.cc, type = "functional", thresh = 0.25)
fibrotic.cc <- netEmbedding(fibrotic.cc, type = "functional")
fibrotic.cc <- netClustering(fibrotic.cc, type = "functional", k = 4)
netVisual_embedding(fibrotic.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
fibrotic.cc <- computeNetSimilarity(fibrotic.cc, type = "structural", thresh = 0.25)
fibrotic.cc <- netEmbedding(fibrotic.cc, type = "structural")
fibrotic.cc <- netClustering(fibrotic.cc, type = "structural")
netVisual_embedding(fibrotic.cc, type = "structural", label.size = 3.5)

### We were advised to run CellChat on just fibroblasts and immune cells, so that's what we're gonna do. 
fibrotic.subsetted <- subset(fibrotic, idents = c("T cell I", "Macrophage I", "Myofibroblast", "Fibroblast I",
                                                  "Dendritic cell I", "T cell II", "T cell III", "Fibroblast II",
                                                  "Macrophage II", "Mast cell I", "Dendritic cell II", "Mast cell II",
                                                  "Fibroblast III"))

# Now run through the CellChat analysis on the subsetted data
fibrotic.subsetted.labels <- Idents(fibrotic.subsetted)
fibrotic.subsetted.identity <- data.frame(group = fibrotic.subsetted.labels, row.names = names(fibrotic.subsetted.labels)) # Dataframe of the cell labels
fibrotic.subsetted.data.input <- GetAssayData(fibrotic.subsetted, assay = "RNA", slot = "data")  # Normalised data matrix

fibrotic.subsetted.cc <- createCellChat(data = fibrotic.subsetted.data.input, do.sparse = F)

# Add the meta-data from Seurat
fibrotic.subsetted.cc <- addMeta(fibrotic.subsetted.cc, meta = fibrotic.subsetted.identity, meta.name = "labels") 
fibrotic.subsetted.cc <- setIdent(fibrotic.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

fibroticSubsettedGroupSize <- as.numeric(table(fibrotic.subsetted.cc@idents)) # Get the number of cells in each group

fibrotic.subsetted.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
fibrotic.subsetted.cc <- subsetData(fibrotic.subsetted.cc) # We subset the expression data of signalling genes to save on computational cost
fibrotic.subsetted.cc <- identifyOverExpressedGenes(fibrotic.subsetted.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
fibrotic.subsetted.cc <- identifyOverExpressedInteractions(fibrotic.subsetted.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# fibrotic.subsetted.cc <- projectData(fibrotic.subsetted.cc, PPI.mouse) # Other option includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
fibrotic.subsetted.cc <- computeCommunProb(fibrotic.subsetted.cc) 
fibrotic.subsetted.cc <- computeCommunProbPathway(fibrotic.subsetted.cc) # Calculate the probabilities at the signalling level
fibrotic.subsetted.cc <- aggregateNet(fibrotic.subsetted.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(3, 4, 8, 13) # Focus on the fibroblasts
pathways.show <- "MK" # All possible pathways are stored in @netP$pathways
fibroticSubsettedGroupSize <- as.numeric(table(fibrotic.subsetted.cc@idents)) # Get the number of cells in each group
netVisual_aggregate(fibrotic.subsetted.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = fibroticSubsettedGroupSize, pt.title = 16)

netAnalysis_contribution(fibrotic.subsetted.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
fibrotic.subsetted.cc <- netAnalysis_signalingRole(fibrotic.subsetted.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(fibrotic.subsetted.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
fibrotic.subsetted.cc <- identifyCommunicationPatterns(fibrotic.subsetted.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(fibrotic.subsetted.cc, pattern = "outgoing", font.size = 3) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(fibrotic.subsetted.cc, pattern = "outgoing", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticsubsettedoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
fibrotic.subsetted.cc <- identifyCommunicationPatterns(fibrotic.subsetted.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(fibrotic.subsetted.cc, pattern = "incoming", font.size = 3)

netAnalysis_dot(fibrotic.subsetted.cc, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticsubsettedincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
fibrotic.subsetted.cc <- computeNetSimilarity(fibrotic.subsetted.cc, type = "functional", thresh = 0.25)
fibrotic.subsetted.cc <- netEmbedding(fibrotic.subsetted.cc, type = "functional")
fibrotic.subsetted.cc <- netClustering(fibrotic.subsetted.cc, type = "functional", k = 4)
netVisual_embedding(fibrotic.subsetted.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticsubsettedfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
fibrotic.subsetted.cc <- computeNetSimilarity(fibrotic.subsetted.cc, type = "structural", thresh = 0.25)
fibrotic.subsetted.cc <- netEmbedding(fibrotic.subsetted.cc, type = "structural")
fibrotic.subsetted.cc <- netClustering(fibrotic.subsetted.cc, type = "structural")
netVisual_embedding(fibrotic.subsetted.cc, type = "structural", label.size = 3.5)

########################################################################################################################
### Now we cluster the regenerative data.
dataDirectory = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/Regenerative/"

regenerative.data <- Read10X(data.dir=dataDirectory)

regenerative <- CreateSeuratObject(counts = regenerative.data, project="Regenerative", min.cells = 3, min.features = 200)

# Check for low-quality cells by counting the percentage of mitochondrial gene expression
regenerative[["percent.mt"]] <- PercentageFeatureSet(regenerative, pattern="^mt-")

# Visualise the MT-gene count to check
VlnPlot(regenerative, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

regenerative <- subset(regenerative, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # Subset the data

# Normalise and scale the data
regenerative <- NormalizeData(regenerative, verbose = FALSE)
regenerative <- FindVariableFeatures(regenerative, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

regenerative <- ScaleData(regenerative, verbose = FALSE)

regenerative <- RunPCA(regenerative, npcs = 40, verbose = FALSE) # Run the PCA

# We need to calculate hte significant PCs to determine the dimensionality of the data
regenerative <- JackStraw(regenerative, num.replicate = 100)
regenerative <- ScoreJackStraw(regenerative, dims = 1:40)

# Visualise the results to figure out how many PCs we should consier.
JackStrawPlot(fibrotic, dims = 1:40)
ElbowPlot(regenerative, ndims = 40)

# A more quantitative way of determining dimensionality, as explained explained here https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# Determine percent of variation associated with each PC
pct <- regenerative[["pca"]]@stdev / sum(regenerative[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_regenerative <- min(co1, co2)

# FInd the neighbours and cluster the data. We may play with the numbers here a bit to see what works well.
regenerative <- FindNeighbors(regenerative, reduction = "pca", dims = 1:num_pcs_regenerative)

# We're going to calculate the clusters at varying resolutions
regenerative <- FindClusters(regenerative, resolution = c(0.3, 0.35, 0.45))

# Run the UMAP
regenerative <- RunUMAP(regenerative, reduction = "pca", dims = 1:num_pcs_regenerative)

# Plot the UMAP at varying resolutions (0.4 is the winner at the moment)
Idents(object = regenerative) <- "RNA_snn_res.0.3"
DimPlot(regenerative, reduction="umap", label=TRUE)

# Let's find the markers for this resolution
regenerative.markers <- FindAllMarkers(regenerative, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top marker from each cluster to see how well these distinguish teh clusters
regenerative.top1.markers <- regenerative.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) # Get the top 5 markers for each cluster

FeaturePlot(regenerative, features = regenerative.top1.markers$gene) # Plot the top marker gene for each cluster

regenerative.top10.markers <- regenerative.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
regenerative.top50.markers <- regenerative.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

DoHeatmap(regenerative, features = regenerative.top10.markers$gene) + NoLegend()

regenerative.old.idents <- levels(regenerative) # Retain the old labels

# We relabel the clusters according to the cell types we identify by marker genes. From googling and literature searching, we've deduced:
# 0 - Fibroblast I, 1 - Fibroblast II, 2 - T cell I, 3 - Myofibroblast, 4 - Macrophage I
# 5 - Fibroblast III, 6 - Macrophage II, 7 - Fibroblast IV, 8 - Schwann cell,
# 9 - Lymphatic endothelial cell, 10 - Myofibroblast II, 11 - Dendritic cell
regenerative.new.idents <- c("Fibroblast I", "Fibroblast II", "T cell I", "Myofibroblast I",
                             "Macrophage I", "Fibroblast III", "Macrophage II", "Fibroblast IV",
                             "Schwann cell", "Lymphatic endothelial cell", "Myofibroblast II", "Dendritic cell")
names(regenerative.new.idents) <- regenerative.old.idents

regenerative <- RenameIdents(regenerative, regenerative.new.idents) # Rename the clusters

save(regenerative, regenerative.new.idents, regenerative.top10.markers, regenerative.markers, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerative.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerative.rdata")

# We're going to plot some things for the preliminary analysis
regenerative_cluster_counts <- table(Idents(regenerative))
num_regenerative_cells <- sum(table(Idents(regenerative)))
regenerative_cluster_proportions <- 100 * regenerative_cluster_counts / num_regenerative_cells # Calculates the percentage of the total population that each cluster is

# Plot the umap and save it for the slides
DimPlot(regenerative, reduction="umap") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeumap.eps', device='eps')

### Now we're in a position to run CellChat on the regenerative dataset
regenerative.labels <- Idents(regenerative)
regenerative.identity <- data.frame(group = regenerative.labels, row.names = names(regenerative.labels)) # Dataframe of the cell labels
regenerative.data.input <- GetAssayData(regenerative, assay = "RNA", slot = "data")  # Normalised data matrix

regenerative.cc <- createCellChat(data = regenerative.data.input, do.sparse = F)

# Add the meta-data from Seurat
regenerative.cc <- addMeta(regenerative.cc, meta = regenerative.identity, meta.name = "labels") 
regenerative.cc <- setIdent(regenerative.cc, ident.use = "labels") # Set the labels to be the default cell identity


regenerative.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
regenerative.cc <- subsetData(regenerative.cc) # We subset the expression data of signalling genes to save on computational cost
regenerative.cc <- identifyOverExpressedGenes(regenerative.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
regenerative.cc <- identifyOverExpressedInteractions(regenerative.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# regenerative.cc <- projectData(regenerative.cc, PPI.mouse) # Other option includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
regenerative.cc <- computeCommunProb(regenerative.cc) 
regenerative.cc <- computeCommunProbPathway(regenerative.cc) # Calculate the probabilities at the signalling level
regenerative.cc <- aggregateNet(regenerative.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 4, 6, 8, 11) # Focus on the fibroblasts
pathways.show <- "MK" # All possible pathways are stored in @netP$pathways
regenerativeGroupSize <- as.numeric(table(regenerative.cc@idents)) # Get the number of cells in each group
netVisual_aggregate(regenerative.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = regenerativeGroupSize, pt.title=16)
netAnalysis_contribution(regenerative.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
regenerative.cc <- netAnalysis_signalingRole(regenerative.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(regenerative.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Let us calculate the signalling pattern groups now
# Now we identify outgoing communication patterns
nPatterns <- 5
regenerative.cc <- identifyCommunicationPatterns(regenerative.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(regenerative.cc, pattern = "outgoing", font.size = 3, font.size.title = 12) 


# Visualize the communication pattern using dot plot
netAnalysis_dot(regenerative.cc, pattern = "outgoing", shape = 0, dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 4
regenerative.cc <- identifyCommunicationPatterns(regenerative.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(regenerative.cc, pattern = "incoming", font.size = 3, font.size.title = 12)


netAnalysis_dot(regenerative.cc, shape = 0, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
regenerative.cc <- computeNetSimilarity(regenerative.cc, type = "functional", thresh = 0.25)
regenerative.cc <- netEmbedding(regenerative.cc, type = "functional")
regenerative.cc <- netClustering(regenerative.cc, type = "functional", k = 4)
netVisual_embedding(regenerative.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativefunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
regenerative.cc <- computeNetSimilarity(regenerative.cc, type = "structural", thresh = 0.25)
regenerative.cc <- netEmbedding(regenerative.cc, type = "structural")
regenerative.cc <- netClustering(regenerative.cc, type = "structural")
netVisual_embedding(regenerative.cc, type = "structural", label.size = 3.5)

### Save the CellChat objects for the FULL datasets
save(fibrotic.cc, regenerative.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020cellchat.rdata")
save(fibrotic.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticallcellchat.rdata")
save(fibrotic.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticsubsettedcellchat.rdata")
save(regenerative.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativeallcellchat.rdata")
save(regenerative.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativesubsettedcellchat.rdata")
save(merged.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/mergedallcellchat.rdata")
save(merged.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/mergedsubsettedcellchat.rdata")

### We will run cell chat on a subset of the cell data---just fibroblasts and immune cells
regenerative.subsetted <- subset(regenerative, idents = c("Fibroblast I", "Fibroblast II", "T cell I", "Myofibroblast I",
                                                          "Macrophage I", "Fibroblast III", "Macrophage II", "Fibroblast IV",  "Myofibroblast II", "Dendritic cell"))

### Now we're in a position to run CellChat on the regenerative dataset
regenerative.subsetted.labels <- Idents(regenerative.subsetted)
regenerative.subsetted.identity <- data.frame(group = regenerative.subsetted.labels, row.names = names(regenerative.subsetted.labels)) # Dataframe of the cell labels
regenerative.subsetted.data.input <- GetAssayData(regenerative.subsetted, assay = "RNA", slot = "data")  # Normalised data matrix

regenerative.subsetted.cc <- createCellChat(data = regenerative.subsetted.data.input, do.sparse = F)

# Add the meta-data from Seurat
regenerative.subsetted.cc <- addMeta(regenerative.subsetted.cc, meta = regenerative.subsetted.identity, meta.name = "labels") 
regenerative.subsetted.cc <- setIdent(regenerative.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

regenerativeSubsettedGroupSize <- as.numeric(table(regenerative.subsetted.cc@idents)) # Get the number of cells in each group

regenerative.subsetted.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
regenerative.subsetted.cc <- subsetData(regenerative.subsetted.cc) # We subset the expression data of signalling genes to save on computational cost
regenerative.subsetted.cc <- identifyOverExpressedGenes(regenerative.subsetted.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
regenerative.subsetted.cc <- identifyOverExpressedInteractions(regenerative.subsetted.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# regenerative.subsetted.cc <- projectData(regenerative.subsetted.cc, PPI.mouse) # Other option includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
regenerative.subsetted.cc <- computeCommunProb(regenerative.subsetted.cc) 
regenerative.subsetted.cc <- computeCommunProbPathway(regenerative.subsetted.cc) # Calculate the probabilities at the signalling level
regenerative.subsetted.cc <- aggregateNet(regenerative.subsetted.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 4, 6, 8, 9) # Focus on the fibroblasts
pathways.show <- "MK" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(regenerative.subsetted.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = regenerativeSubsettedGroupSize, pt.title=16)
netAnalysis_contribution(regenerative.subsetted.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
regenerative.subsetted.cc <- netAnalysis_signalingRole(regenerative.subsetted.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(regenerative.subsetted.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Let us calculate the signalling pattern groups now
# Now we identify outgoing communication patterns
nPatterns <- 5
regenerative.subsetted.cc <- identifyCommunicationPatterns(regenerative.subsetted.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(regenerative.subsetted.cc, pattern = "outgoing", font.size = 3, font.size.title = 12) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(regenerative.subsetted.cc, pattern = "outgoing", shape = 0, dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativesubsettedoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
regenerative.subsetted.cc <- identifyCommunicationPatterns(regenerative.subsetted.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(regenerative.subsetted.cc, pattern = "incoming", font.size = 3)

netAnalysis_dot(regenerative.subsetted.cc, shape = 0, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativesubsettedincomingcommunicationpatterns.eps')


# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
regenerative.subsetted.cc <- computeNetSimilarity(regenerative.subsetted.cc, type = "functional", thresh = 0.25)
regenerative.subsetted.cc <- netEmbedding(regenerative.subsetted.cc, type = "functional")
regenerative.subsetted.cc <- netClustering(regenerative.subsetted.cc, type = "functional", k = 4)
netVisual_embedding(regenerative.subsetted.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativesubsettedfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
regenerative.subsetted.cc <- computeNetSimilarity(regenerative.subsetted.cc, type = "structural", thresh = 0.25)
regenerative.subsetted.cc <- netEmbedding(regenerative.subsetted.cc, type = "structural")
regenerative.subsetted.cc <- netClustering(regenerative.subsetted.cc, type = "structural", k = 3)
netVisual_embedding(regenerative.subsetted.cc, type = "structural", label.size = 3.5)

save(fibrotic.subsetted.cc, regenerative.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020cellchatsubsetted.rdata")

################################################################################################################################################
### Load the cellchat objects so that we can perform a joint manifold analysis to determine functional and structural similarity
load(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020cellchat.rdata")
load(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020cellchatsubsetted.rdata")

# Let's first analyse the full dataset
merged.cc <- mergeCellChat(list(fibrotic.cc, regenerative.cc), add.names = c("Fibrotic", "Regenerative"))

# Let's first analyse for structural similarity, and then see if we can run the same on functional similarity (I doubt it)
merged.cc <- computeNetSimilarityPairwise(merged.cc, type = "structural")
merged.cc <- netEmbedding(merged.cc, type = "structural")
merged.cc <- netClustering(merged.cc, type = "structural")

# Visualise the structural similarity
netVisual_embeddingPairwise(merged.cc, type = "structural", dot.size = c(5, 15), label.size = 6) + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityjointdata.eps')
netVisual_embeddingPairwiseZoomIn(merged.cc, type = "structural", dot.size = c(3, 9), label.size = 4) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityzoomedinjointdata.eps')

# Going to see if we can do the same for functional similarity
merged.cc <- computeNetSimilarityPairwise(merged.cc, type = "functional") # Nope can't be done.

# Visualise the pathway distance in the learnt joint manifold
rankSimilarity(merged.cc, type = "structural") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/rankedpathwaydistancesjointdata.eps')
# Identify and visualise the conserved and context-specified signalling pathways
rankNet(merged.cc, mode = "comparison") +  theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/informationflowjointdata.eps')

# Now let's analyse the subsetted datasets
merged.subsetted.cc <- mergeCellChat(list(fibrotic.subsetted.cc, regenerative.subsetted.cc), add.names = c("Fibrotic", "Regenerative"))

# Let's first analyse for structural similarity, and then see if we can run the same on functional similarity (I doubt it)
merged.subsetted.cc <- computeNetSimilarityPairwise(merged.subsetted.cc, type = "structural")
merged.subsetted.cc <- netEmbedding(merged.subsetted.cc, type = "structural")
merged.subsetted.cc <- netClustering(merged.subsetted.cc, type = "structural")

# Visualise the structural similarity
netVisual_embeddingPairwise(merged.subsetted.cc, type = "structural", dot.size = c(5, 15), label.size = 6) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityjointsubsetteddata.eps')
netVisual_embeddingPairwiseZoomIn(merged.subsetted.cc, type = "structural", dot.size = c(3, 9), label.size = 4) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityzoomedinjointsubsetteddata.eps')

# Visualise the pathway distance in the learnt joint manifold
rankSimilarity(merged.subsetted.cc, type = "structural") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/rankedpathwaydistancesjointsubsetteddata.eps')

# Identify and visualise the conserved and context-specified signalling pathways
rankNet(merged.subsetted.cc, mode = "comparison") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/informationflowjointsubsetteddata.eps')

# Plot the incoming and outgoing comunication patterns
netAnalysis_dot(merged.cc, pattern = "incoming")

save(merged.cc, merged.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020cellchatjointdata.rdata")

###############################################################################################################################################################
### I want to integrate the two datasets under the different conditions to understand the conserved markers and clusters and help with cell type identification.
gay2020.list <- NULL
conditions <- c("Fibrotic", "Regenerative")
dataDirectory = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/"

for (i in 1:length(conditions))
{
  data.directory <- paste(dataDirectory, conditions[[i]], "/", sep="")
  data <- Read10X(data.dir=data.directory)
  seurat.obj <- CreateSeuratObject(counts = data, project=conditions[[i]], min.cells = 3, min.features = 200)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern="^mt-")
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # Subset the data
  seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat.obj[['condition']] <- conditions[[i]]
  gay2020.list[[i]] <- seurat.obj
}

gay2020.anchors <- FindIntegrationAnchors(object.list = gay2020.list, dims = 1:30)

# Integrate the data
gay2020.integrated <- IntegrateData(anchorset = gay2020.anchors, dims = 1:30)

# Switch to integrated array
DefaultAssay(gay2020.integrated) <- "integrated"

gay2020.integrated <- ScaleData(gay2020.integrated, verbose = FALSE)
gay2020.integrated <- RunPCA(gay2020.integrated, npcs = 30, verbose = FALSE)

pct <- gay2020.integrated[["pca"]]@stdev / sum(gay2020.integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_integrated <- min(co1, co2)

# FInd the neighbours and cluster the data. We may play with the numbers here a bit to see what works well.
gay2020.integrated <- FindNeighbors(gay2020.integrated, reduction = "pca", dims = 1:num_pcs_integrated)

# Run the UMAP
gay2020.integrated <- RunUMAP(gay2020.integrated, reduction = "pca", dims = 1:num_pcs_integrated)
DimPlot(gay2020.integrated, reduction="umap")

# Let's cluster the data at varying resolutions to check what would be the best
gay2020.integrated <- FindClusters(gay2020.integrated, resolution = seq(0.3, 0.5, 0.05))

# Now let's determine hte optimal resolution (0.4 seems to be the winner)
Idents(object = gay2020.integrated) <- "integrated_snn_res.0.4"

# Let's find the markers for this resolution
gay2020.integrated.markers <- FindAllMarkers(gay2020.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top marker from each cluster to see how well these distinguish teh clusters

gay2020.integrated.top10.markers <- gay2020.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gay2020.integrated.top50.markers <- gay2020.integrated.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(gay2020.integrated, features = gay2020.integrated.top10.markers$gene) + NoLegend()

gay2020.integrated.old.idents <- levels(gay2020.integrated) # Retain the old labels

# We relabel the clusters according to the cell types we identify by marker genes. From googling and literature searching, we've deduced:
# 0 - Fibroblast I, 1 - Fibroblast II, 2 - Myofibroblast I, 3 - Macrophage I, 4 - Dendritic cell I
# 5 - Fibroblast III, 6 - T cell I, 7 - Dendritic cell II, 8 - Fibroblast IV, 9 - T cell II
# 10 - Fibroblast V, 11 -Macrphage II, 12 - T cell III, - 13 -Schwann cell I, 14 - Lymphatic endothelial cell
# 15 - Dendritic cell III, 16 - Mast cell, 17 - Schwann cell II, 18 - Dendritic cell IV, 19 - Myofibroblast II
gay2020.integrated.new.idents <- c("Fibroblast I", "Fibroblast II", "Myofibroblast I", "Macrophage I", "Dendritic cell I",
                                 "Fibroblast III", "T cell I", "Dendritic cell II", "Fibroblast IV", "T cell II",
                                 "Fibroblast V", "Macrophage II", "T cell III", "Schwann cell", "Lymphatic endothelial cell",
                                 "Dendritic cell III", "Mast cell", "Schwann cell II", "Dendritic cell IV", "Myofibroblast II")
names(gay2020.integrated.new.idents) <- gay2020.integrated.old.idents

gay2020.integrated <- RenameIdents(gay2020.integrated, gay2020.integrated.new.idents)

DimPlot(gay2020.integrated, reduction = "umap", label=TRUE)

save(gay2020.integrated, gay2020.integrated.new.idents, gay2020.integrated.top10.markers, gay2020.integrated.markers, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020integrated.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020integrated.rdata")

###############################################################################################################################################################
### Load the original clustered data, now that Max has given it to us.
###############################################################################################################################################################

load("~/Documents/CellCellCommunicationModelling/Data/Gay2020/wihnfromdenise.robj")

gay2020.orig.markers <- FindAllMarkers(wihnint, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gay2020.orig.top10.markers <- gay2020.orig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
gay2020.orig.top50.markers <- gay2020.orig.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

DoHeatmap(wihnint, features = gay2020.orig.top10.markers$gene) + NoLegend()

gay2020.orig.old.idents <- levels(wihnint)
# Rename the cell clusters using the cell types from Max's SkinGenes website.
gay2020.orig.new.idents <- c("Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-4", "Macrophages/DCs-1", "T cells-1",
                             "Fibroblasts-5", "Macrophages/DCs-2", "Fibroblasts-6", "T cells-2", "T cells-3", "Endothelial", 
                             "Macrophages/DCs-3", "Schwann cells", "Macrophages/DCs-4", "Macrophages/DCs-5", "Macrophages/DCs-6",
                             "Other immune cells")
names(gay2020.orig.new.idents) <- gay2020.orig.old.idents

# Rename the old clusters
wihnint <- RenameIdents(wihnint, gay2020.orig.new.idents)

# These are alternative labellings (just in case)
gay2020.orig.alt.idents <- c("Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Myofibroblasts", "Macrophages-1", "T cells-1",
                             "Fibroblasts-4", "Dendritic cells-1", "Fibroblasts-5", "T cells-2", "T cells-3", "Endothelial", 
                             "Macrophages-2", "Schwann cells", "Macrophages-3", "Dendritic cells-2", "Macrophages-4",
                             "Mast cells")

gay2020.orig.top50.markers.renamed <- gay2020.orig.top50.markers
levels(gay2020.orig.top50.markers.renamed$cluster) <- gay2020.orig.new.idents

# write.csv(gay2020.orig.top50.markers.renamed, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/Top50Markers.csv")

# Run a quick UMAP
pct <- wihnint[["pca"]]@stdev / sum(wihnint[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

num_pcs_integrated <- min(co1, co2)

# Run the UMAP
wihnint <- RunUMAP(wihnint, reduction = "pca", dims = 1:num_pcs_integrated)

# Split the original data into the fibrotic and merged
fibrotic <- subset(wihnint, subset = stim == "CTRL")
fibrotic.subsetted <- subset(fibrotic, idents = c("Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-4", "Macrophages/DCs-1", "T cells-1",
                                                  "Fibroblasts-5", "Macrophages/DCs-2", "Fibroblasts-6", "T cells-2", "T cells-3", 
                                                  "Macrophages/DCs-3", "Macrophages/DCs-4", "Macrophages/DCs-5", "Macrophages/DCs-6",
                                                  "Other immune cells"))

regenerative <- subset(wihnint, subset = stim == "STIM")
regenerative.subsetted <- subset(regenerative, idents = c("Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-4", "Macrophages/DCs-1", "T cells-1",
                                                          "Fibroblasts-5", "Macrophages/DCs-2", "Fibroblasts-6", "T cells-2", "T cells-3", 
                                                          "Macrophages/DCs-3", "Macrophages/DCs-4", "Macrophages/DCs-5", "Macrophages/DCs-6"))

####################################################################################################################################
#### Create the new CellChat objects
#######################################################################################################################################

# First the fibrotic data, full and subsetted
fibrotic.labels <- Idents(fibrotic)
fibrotic.identity <- data.frame(group = fibrotic.labels, row.names = names(fibrotic.labels)) # Dataframe of the cell labels
fibrotic.data.input <- GetAssayData(fibrotic, assay = "RNA", slot = "data")  # Normalised data matrix

fibrotic.cc <- createCellChat(data = fibrotic.data.input, do.sparse = F)

# Add the meta-data from Seurat
fibrotic.cc <- addMeta(fibrotic.cc, meta = fibrotic.identity, meta.name = "labels") 
fibrotic.cc <- setIdent(fibrotic.cc, ident.use = "labels") # Set the labels to be the default cell identity

fibroticGroupSize <- as.numeric(table(fibrotic.cc@idents)) # Get the number of cells in each group

fibrotic.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# Subsetted data
fibrotic.subsetted.labels <- Idents(fibrotic.subsetted)
fibrotic.subsetted.identity <- data.frame(group = fibrotic.subsetted.labels, row.names = names(fibrotic.subsetted.labels)) # Dataframe of the cell labels
fibrotic.subsetted.data.input <- GetAssayData(fibrotic.subsetted, assay = "RNA", slot = "data")  # Normalised data matrix

fibrotic.subsetted.cc <- createCellChat(data = fibrotic.subsetted.data.input, do.sparse = F)

# Add the meta-data from Seurat
fibrotic.subsetted.cc <- addMeta(fibrotic.subsetted.cc, meta = fibrotic.subsetted.identity, meta.name = "labels") 
fibrotic.subsetted.cc <- setIdent(fibrotic.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

fibroticSubsettedGroupSize <- as.numeric(table(fibrotic.subsetted.cc@idents)) # Get the number of cells in each group

fibrotic.subsetted.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# Now the full regenerative data
regenerative.labels <- Idents(regenerative)
regenerative.identity <- data.frame(group = regenerative.labels, row.names = names(regenerative.labels)) # Dataframe of the cell labels
regenerative.data.input <- GetAssayData(regenerative, assay = "RNA", slot = "data")  # Normalised data matrix

regenerative.cc <- createCellChat(data = regenerative.data.input, do.sparse = F)

# Add the meta-data from Seurat
regenerative.cc <- addMeta(regenerative.cc, meta = regenerative.identity, meta.name = "labels") 
regenerative.cc <- setIdent(regenerative.cc, ident.use = "labels") # Set the labels to be the default cell identity

regenerativeGroupSize <- as.numeric(table(regenerative.cc@idents)) # Get the number of cells in each group

regenerative.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

# And finally the subsetted data
regenerative.subsetted.labels <- Idents(regenerative.subsetted)
regenerative.subsetted.identity <- data.frame(group = regenerative.subsetted.labels, row.names = names(regenerative.subsetted.labels)) # Dataframe of the cell labels
regenerative.subsetted.data.input <- GetAssayData(regenerative.subsetted, assay = "RNA", slot = "data")  # Normalised data matrix

regenerative.subsetted.cc <- createCellChat(data = regenerative.subsetted.data.input, do.sparse = F)

# Add the meta-data from Seurat
regenerative.subsetted.cc <- addMeta(regenerative.subsetted.cc, meta = regenerative.identity, meta.name = "labels") 
regenerative.subsetted.cc <- setIdent(regenerative.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

regenerativeSubsettedGroupSize <- as.numeric(table(regenerative.subsetted.cc@idents)) # Get the number of cells in each group

regenerative.subsetted.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

### Let's now proceed with the CellChat analysis

