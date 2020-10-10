# Load the relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)

### We load the relevant CellChat databases.
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

# Load the original data that Max gave us
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

# Rename to fit the cell names on CellChat
fibrotic.new.idents <- c("FIB-1", "FIB-2", "FIB-3", "FIB-4", "M/DC-1", "T-1", "FIB-5", "M/DC-2", "FIB-6", "T-2", "T-3", "ENDO", "M/DC-3", "SCHWANN", "M/DC-4", "M/DC-5", "M/DC-6", "OTHER IC")
names(fibrotic.new.idents) <- levels(Idents(fibrotic))
fibrotic <- RenameIdents(fibrotic, fibrotic.new.idents)

fibrotic.subsetted.new.idents <- c("FIB-1", "FIB-2", "FIB-3", "FIB-4", "M/DC-1", "T-1", "FIB-5", "M/DC-2", "FIB-6", "T-2", "T-3", "M/DC-3", "M/DC-4", "M/DC-5", "M/DC-6", "OTHER IC")
names(fibrotic.subsetted.new.idents) <- levels(Idents(fibrotic.subsetted))
fibrotic.subsetted <- RenameIdents(fibrotic.subsetted, fibrotic.subsetted.new.idents)

# Also generate the UMAP for this data
DimPlot(fibrotic, reduction = "umap") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticallumap.eps', device='eps')

regenerative <- subset(wihnint, subset = stim == "STIM")
regenerative.subsetted <- subset(regenerative, idents = c("Fibroblasts-1", "Fibroblasts-2", "Fibroblasts-3", "Fibroblasts-4", "Macrophages/DCs-1", "T cells-1",
                                                          "Fibroblasts-5", "Macrophages/DCs-2", "Fibroblasts-6", "T cells-2", "T cells-3", 
                                                          "Macrophages/DCs-3", "Macrophages/DCs-4", "Macrophages/DCs-5", "Macrophages/DCs-6"))

# Rename to fit the cell names on CellChat
regenerative.new.idents <- c("FIB-1", "FIB-2", "FIB-3", "FIB-4", "M/DC-1", "T-1", "FIB-5", "M/DC-2", "FIB-6", "T-2", "T-3", "ENDO", "M/DC-3", "SCHWANN", "M/DC-4", "M/DC-5", "M/DC-6")
names(regenerative.new.idents) <- levels(Idents(regenerative))
regenerative <- RenameIdents(regenerative, regenerative.new.idents)

regenerative.subsetted.new.idents <- c("FIB-1", "FIB-2", "FIB-3", "FIB-4", "M/DC-1", "T-1", "FIB-5", "M/DC-2", "FIB-6", "T-2", "T-3", "M/DC-3", "M/DC-4", "M/DC-5", "M/DC-6")
names(regenerative.subsetted.new.idents) <- levels(Idents(regenerative.subsetted))
regenerative.subsetted <- RenameIdents(regenerative.subsetted, regenerative.subsetted.new.idents)

# Also generate the UMAP for this data
DimPlot(regenerative, reduction = "umap") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeallumap.eps', device='eps')

save(wihnint, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/wihnfromdenise.robj")


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
regenerative.subsetted.cc <- addMeta(regenerative.subsetted.cc, meta = regenerative.subsetted.identity, meta.name = "labels") 
regenerative.subsetted.cc <- setIdent(regenerative.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

regenerativeSubsettedGroupSize <- as.numeric(table(regenerative.subsetted.cc@idents)) # Get the number of cells in each group

regenerative.subsetted.cc@DB <- CellChatDB.use # Set the database for the fibrotic data

####################################################################################################################################
### Let's now proceed with the CellChat analysis
####################################################################################################################################

### First the fibrotic dataset

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
fibrotic.cc <- subsetData(fibrotic.cc) # We subset the expression data of signalling genes to save on computational cost
fibrotic.cc <- identifyOverExpressedGenes(fibrotic.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
fibrotic.cc <- identifyOverExpressedInteractions(fibrotic.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor

fibrotic.cc <- computeCommunProb(fibrotic.cc, raw.use = TRUE, population.size = FALSE) # These parameters were at the suggestion of Suoqin
fibrotic.cc <- computeCommunProbPathway(fibrotic.cc) # Calculate the probabilities at the signalling level
fibrotic.cc <- aggregateNet(fibrotic.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 3, 4, 7) # Focus on the fibroblasts
pathways.show <- "SEMA3" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(fibrotic.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, vertex.size = fibroticGroupSize, pt.title = 18)
netAnalysis_contribution(fibrotic.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
fibrotic.cc <- netAnalysis_signalingRole(fibrotic.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(fibrotic.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
fibrotic.cc <- identifyCommunicationPatterns(fibrotic.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(fibrotic.cc, pattern = "outgoing", font.size = 3.25, font.size.title = 14) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(fibrotic.cc, pattern = "outgoing", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticallorigoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
fibrotic.cc <- identifyCommunicationPatterns(fibrotic.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(fibrotic.cc, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(fibrotic.cc, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticallorigincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
fibrotic.cc <- computeNetSimilarity(fibrotic.cc, type = "functional", thresh = 0.25)
fibrotic.cc <- netEmbedding(fibrotic.cc, type = "functional")
fibrotic.cc <- netClustering(fibrotic.cc, type = "functional", k = 4)
netVisual_embedding(fibrotic.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticallorigfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
fibrotic.cc <- computeNetSimilarity(fibrotic.cc, type = "structural", thresh = 0.25)
fibrotic.cc <- netEmbedding(fibrotic.cc, type = "structural")
fibrotic.cc <- netClustering(fibrotic.cc, type = "structural")
netVisual_embedding(fibrotic.cc, type = "structural", label.size = 3.5)

saveRDS(fibrotic.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticallcellchat.rds")
fibrotic.cc <- readRDS(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticallcellchat.rds")

### Now the subsetted fibrotic data
fibrotic.subsetted.cc <- subsetData(fibrotic.subsetted.cc) # We subset the expression data of signalling genes to save on computational cost
fibrotic.subsetted.cc <- identifyOverExpressedGenes(fibrotic.subsetted.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
fibrotic.subsetted.cc <- identifyOverExpressedInteractions(fibrotic.subsetted.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor

# We now infer the cell-cell communication network by calculating the communication probabilities
fibrotic.subsetted.cc <- computeCommunProb(fibrotic.subsetted.cc, raw.use = TRUE, population.size = FALSE) # These parameters were at the suggestion of Suoqin
fibrotic.subsetted.cc <- computeCommunProbPathway(fibrotic.subsetted.cc) # Calculate the probabilities at the signalling level
fibrotic.subsetted.cc <- aggregateNet(fibrotic.subsetted.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 3, 4, 7) # Focus on the fibroblasts
pathways.show <- "SEMA3" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(fibrotic.subsetted.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = fibroticSubsettedGroupSize, pt.title = 16)

netAnalysis_contribution(fibrotic.subsetted.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
fibrotic.subsetted.cc <- netAnalysis_signalingRole(fibrotic.subsetted.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(fibrotic.subsetted.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
fibrotic.subsetted.cc <- identifyCommunicationPatterns(fibrotic.subsetted.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(fibrotic.subsetted.cc, pattern = "outgoing", font.size = 3.25) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(fibrotic.subsetted.cc, pattern = "outgoing", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticorigsubsettedoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
fibrotic.subsetted.cc <- identifyCommunicationPatterns(fibrotic.subsetted.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(fibrotic.subsetted.cc, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(fibrotic.subsetted.cc, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticorigsubsettedincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
fibrotic.subsetted.cc <- computeNetSimilarity(fibrotic.subsetted.cc, type = "functional", thresh = 0.25)
fibrotic.subsetted.cc <- netEmbedding(fibrotic.subsetted.cc, type = "functional")
fibrotic.subsetted.cc <- netClustering(fibrotic.subsetted.cc, type = "functional", k = 4)
netVisual_embedding(fibrotic.subsetted.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticorigsubsettedfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
fibrotic.subsetted.cc <- computeNetSimilarity(fibrotic.subsetted.cc, type = "structural", thresh = 0.25)
fibrotic.subsetted.cc <- netEmbedding(fibrotic.subsetted.cc, type = "structural")
fibrotic.subsetted.cc <- netClustering(fibrotic.subsetted.cc, type = "structural")
netVisual_embedding(fibrotic.subsetted.cc, type = "structural", label.size = 3.5)

saveRDS(fibrotic.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticsubsettedcellchat.rds")
fibrotic.subsetted.cc <- readRDS(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/fibroticsubsettedcellchat.rds")

### We now do the full regenerative dataset
regenerative.cc <- subsetData(regenerative.cc) # We subset the expression data of signalling genes to save on computational cost
regenerative.cc <- identifyOverExpressedGenes(regenerative.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
regenerative.cc <- identifyOverExpressedInteractions(regenerative.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor

# We now infer the cell-cell communication network by calculating the communication probabilities
regenerative.cc <- computeCommunProb(regenerative.cc, raw.use = TRUE, population.size = FALSE) 
regenerative.cc <- computeCommunProbPathway(regenerative.cc) # Calculate the probabilities at the signalling level
regenerative.cc <- aggregateNet(regenerative.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 3, 4, 7) # Focus on the fibroblasts
pathways.show <- "SEMA3" # All possible pathways are stored in @netP$pathways
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
netAnalysis_river(regenerative.cc, pattern = "outgoing", font.size = 3.25, font.size.title = 12) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(regenerative.cc, pattern = "outgoing", shape = 0, dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
regenerative.cc <- identifyCommunicationPatterns(regenerative.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(regenerative.cc, pattern = "incoming", font.size = 3.25, font.size.title = 12)

netAnalysis_dot(regenerative.cc, shape = 0, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
regenerative.cc <- computeNetSimilarity(regenerative.cc, type = "functional", thresh = 0.25)
regenerative.cc <- netEmbedding(regenerative.cc, type = "functional")
regenerative.cc <- netClustering(regenerative.cc, type = "functional", k = 3)
netVisual_embedding(regenerative.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
regenerative.cc <- computeNetSimilarity(regenerative.cc, type = "structural", thresh = 0.25)
regenerative.cc <- netEmbedding(regenerative.cc, type = "structural")
regenerative.cc <- netClustering(regenerative.cc, type = "structural")
netVisual_embedding(regenerative.cc, type = "structural", label.size = 3.5)

saveRDS(regenerative.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativeallcellchat.rds")
regenerative.cc <- readRDS(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativeallcellchat.rds")

### And lastly, for the individual CellChats, the subsetted regenerative data
regenerative.subsetted.cc <- subsetData(regenerative.subsetted.cc) # We subset the expression data of signalling genes to save on computational cost
regenerative.subsetted.cc <- identifyOverExpressedGenes(regenerative.subsetted.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
regenerative.subsetted.cc <- identifyOverExpressedInteractions(regenerative.subsetted.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor

# We now infer the cell-cell communication network by calculating the communication probabilities
regenerative.subsetted.cc <- computeCommunProb(regenerative.subsetted.cc, raw.use = TRUE, population.size = FALSE) 
regenerative.subsetted.cc <- computeCommunProbPathway(regenerative.subsetted.cc) # Calculate the probabilities at the signalling level
regenerative.subsetted.cc <- aggregateNet(regenerative.subsetted.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 3, 4, 7) # Focus on the fibroblasts
pathways.show <- "SEMA3" # All possible pathways are stored in @netP$pathways
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
netAnalysis_river(regenerative.subsetted.cc, pattern = "outgoing", font.size = 3.25, font.size.title = 12) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(regenerative.subsetted.cc, pattern = "outgoing", shape = 0, dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigsubsettedoutgoingcommunicationpatterns.eps')

# Now let's look at incoming patterns
nPatterns <- 5
regenerative.subsetted.cc <- identifyCommunicationPatterns(regenerative.subsetted.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(regenerative.subsetted.cc, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(regenerative.subsetted.cc, shape = 0, pattern = "incoming", dot.size = c(1.5, 4.5)) + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigsubsettedincomingcommunicationpatterns.eps')

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
regenerative.subsetted.cc <- computeNetSimilarity(regenerative.subsetted.cc, type = "functional", thresh = 0.25)
regenerative.subsetted.cc <- netEmbedding(regenerative.subsetted.cc, type = "functional")
regenerative.subsetted.cc <- netClustering(regenerative.subsetted.cc, type = "functional", k = 3)
netVisual_embedding(regenerative.subsetted.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))  + theme(text=element_text(size=18), axis.text=element_text(size=18)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativeorigsubsettedfunctionalsimilarities.eps')

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
regenerative.subsetted.cc <- computeNetSimilarity(regenerative.subsetted.cc, type = "structural", thresh = 0.25)
regenerative.subsetted.cc <- netEmbedding(regenerative.subsetted.cc, type = "structural")
regenerative.subsetted.cc <- netClustering(regenerative.subsetted.cc, type = "structural", k = 3)
netVisual_embedding(regenerative.subsetted.cc, type = "structural", label.size = 3.5)

saveRDS(regenerative.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativesubsettedcellchat.rds")
regenerative.subsetted.cc <- readRDS(file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/regenerativesubsettedcellchat.rds")

####################################################################################################################################
### We also run CellChat on merged objects of fibrotic and regenerative wound data to compare structural
### (and functional) similarity of pathways
####################################################################################################################################

### Let's first analyse the full dataset
merged.cc <- mergeCellChat(list(fibrotic.cc, regenerative.cc), add.names = c("Fibrotic", "Regenerative"))

# Let's first analyse for structural similarity, and then see if we can run the same on functional similarity (I doubt it)
merged.cc <- computeNetSimilarityPairwise(merged.cc, type = "structural")
merged.cc <- netEmbedding(merged.cc, type = "structural")
merged.cc <- netClustering(merged.cc, type = "structural")

# Visualise the structural similarity
netVisual_embeddingPairwise(merged.cc, type = "structural", dot.size = c(5, 15), label.size = 6) + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityorigjointdata.eps')
netVisual_embeddingPairwiseZoomIn(merged.cc, type = "structural", dot.size = c(3, 9), label.size = 4) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityzoomedinjointdata.eps')

# Going to see if we can do the same for functional similarity
merged.cc <- computeNetSimilarityPairwise(merged.cc, type = "functional") # Nope can't be done.

# Visualise the pathway distance in the learnt joint manifold
rankSimilarity(merged.cc, type = "structural") + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/rankedpathwaydistancesorigjointdata.eps')

# Identify and visualise the conserved and context-specified signalling pathways
rankNet(merged.cc, mode = "comparison")
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/informationfloworigjointdata.eps')

# Now let's analyse the subsetted datasets
merged.subsetted.cc <- mergeCellChat(list(fibrotic.subsetted.cc, regenerative.subsetted.cc), add.names = c("Fibrotic", "Regenerative"))

# Let's first analyse for structural similarity, and then see if we can run the same on functional similarity (I doubt it)
merged.subsetted.cc <- computeNetSimilarityPairwise(merged.subsetted.cc, type = "structural")
merged.subsetted.cc <- netEmbedding(merged.subsetted.cc, type = "structural")
merged.subsetted.cc <- netClustering(merged.subsetted.cc, type = "structural")

# Visualise the structural similarity
netVisual_embeddingPairwise(merged.subsetted.cc, type = "structural", dot.size = c(5, 15), label.size = 6) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityorigjointsubsetteddata.eps')
netVisual_embeddingPairwiseZoomIn(merged.subsetted.cc, type = "structural", dot.size = c(3, 9), label.size = 4) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/structuralsimilarityzoomedinorigjointsubsetteddata.eps')

# Visualise the pathway distance in the learnt joint manifold
rankSimilarity(merged.subsetted.cc, type = "structural") + theme(text=element_text(size=18), axis.text=element_text(size=18))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/rankedpathwaydistancesorigjointsubsetteddata.eps')

# Identify and visualise the conserved and context-specified signalling pathways
rankNet(merged.subsetted.cc, mode = "comparison") + theme(text=element_text(size=18), axis.text=element_text(size=18))
ggsave('~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/informationfloworigjointsubsetteddatafixed.eps')

saveRDS(c(merged.cc, merged.subsetted.cc), file = "~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020origcellchatjointdata.rds")
mergedCellChats <- readRDS("~/Documents/CellCellCommunicationModelling/Data/Gay2020/gay2020origcellchatjointdata.rds")
merged.orig.cc <- mergedCellChats[[1]]
merged.orig.subsetted.cc <- mergedCellChats[[2]]

# We were asked to look at expression of Crabp1, which marks upper fibroblasts
FeaturePlot(fibrotic, features = c("Crabp1")) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/fibroticcrabp1featureplot.eps', device='eps')

FeaturePlot(regenerative, features = c("Crabp1")) + theme(text=element_text(size=24), axis.text=element_text(size=24))
ggsave(filename = '~/Documents/Presentations/PostdocTalks/CellChatAnalysisOfGay2020_Sep2020/regenerativecrabp1featureplot.eps', device='eps')
