### Application of Suoqin Jin's CellChat package to various datas on unwounded and wounded murine skin tissue data.
### We're following the tutorial found here:https://github.com/sqjin/CellChat/blob/master/vignettes/walkthrough_wound.html

# Load the relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # Other options include ECM-Receptor and Cell-Cell Contact

### Let's first look at the unwounded data from Haensel et al. (2020)
load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/bs_cca.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
unwounded.data.input <- GetAssayData(ep, assay = "RNA", slot = "data")  # Normalised data matrix
# unwounded.data.input <- ep@data
unwounded.labels <- Idents(ep)  # Get the cell group labels
# unwounded.identity <- data.frame(group = ep@ident, row.names = names(ep@ident))
unwounded.identity <- data.frame(group = unwounded.labels, row.names = names(unwounded.labels)) # Dataframe of the cell labels

# Create the CellChat object now
unwounded.cc <- createCellChat(data = unwounded.data.input, do.sparse = F)

# Add the meta-data from Seurat
unwounded.cc <- addMeta(unwounded.cc, meta = unwounded.identity, meta.name = "labels") 
unwounded.cc <- setIdent(unwounded.cc, ident.use = "labels") # Set the labels to be the default cell identity

unwoundedGroupSize <- as.numeric(table(unwounded.cc@idents)) # Get the number of cells in each group

# Set the current object to be the 'secreted' cell chat data
unwounded.cc.secreted <- unwounded.cc

unwounded.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
# unwounded.cc <- subsetData(unwounded.cc) # We subset the expression data of signalling genes to save on computational cost
unwounded.cc <- identifyOverExpressedGenes(unwounded.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
unwounded.cc <- identifyOverExpressedInteractions(unwounded.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# unwounded.cc <- projectData(unwounded.cc, PPI.mouse) # Other option includes PPI.human (we've been told maybe we don't need this for this data.)

# We now infer the cell-cell communication network by calculating the communication probabilities
unwounded.cc <- computeCommunProb(unwounded.cc, raw.use = TRUE, population.size = FALSE) 
unwounded.cc <- computeCommunProbPathway(unwounded.cc) # Calculate the probabilities at the signalling level
unwounded.cc <- aggregateNet(unwounded.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- c(1, 2, 3, 4) # Focus on the fibroblast cells
pathways.show <- "FGF" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(unwounded.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = unwoundedGroupSize)
netAnalysis_contribution(unwounded.cc.secreted, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
unwounded.cc.secreted <- netAnalysis_signalingRole(unwounded.cc.secreted, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(unwounded.cc.secreted, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns = 5 
unwounded.cc <- identifyCommunicationPatterns(unwounded.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(unwounded.cc.secreted, pattern = "outgoing") 

# Visualize the communication pattern using dot plot
netAnalysis_dot(unwounded.cc, pattern = "outgoing")

# Now let's look at incoming patterns
nPatterns <- 5
unwounded.cc <- identifyCommunicationPatterns(unwounded.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(unwounded.cc, pattern = "incoming")

netAnalysis_dot(unwounded.cc, pattern = "incoming")

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
unwounded.cc <- computeNetSimilarity(unwounded.cc, type = "functional", thresh = 0.25)
unwounded.cc <- netEmbedding(unwounded.cc, type = "functional")
unwounded.cc <- netClustering(unwounded.cc, type = "functional", k = 5)
netVisual_embedding(unwounded.cc, type = "functional", pathway.remove.show = F, label.size = 3.5)

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
unwounded.cc <- computeNetSimilarity(unwounded.cc, type = "structural", thresh = 0.25)
unwounded.cc <- netEmbedding(unwounded.cc, type = "structural")
unwounded.cc <- netClustering(unwounded.cc, type = "structural")
netVisual_embedding(unwounded.cc, type = "structural", label.size = 3.5)

### We now look at the wounded data to compare the differences in cell-cell communications
load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/sw_cca.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
wounded.data.input <- GetAssayData(ep, assay = "RNA", slot = "data")  # Normalised data matrix
# wounded.data.input <- ep@data
wounded.labels <- Idents(ep)  # Get the cell group labels
levels(wounded.labels) <- c("Epidermal basal", "Proliferative epidermal basal", "Epidermal spinous", "HF/HFSC", "Fibroblast I", "Fibroblast II", "Myofibroblast", "Macrophage I", "Macrophage II", "Macrophage III", "Dendritic/Langerhans cell", "T cell", "Endothelial", "Skeletal muscle")
wounded.identity <- data.frame(group = wounded.labels, row.names = names(wounded.labels)) # Dataframe of the cell labels
# wounded.identity <- data.frame(group = ep@ident, row.names = names(ep@ident))
unique(wounded.identity$group)

ep.old.idents <- levels(ep) # Just in case
ep.new.idents <-c("Epidermal basal", "Proliferative epidermal basal", "Epidermal spinous", "HF/HFSC", "Fibroblast I", "Fibroblast II", "Myofibroblast", "Macrophage I", "Macrophage II", "Macrophage III", "Dendritic/Langerhans cell", "T cell", "Endothelial", "Skeletal muscle")
names(ep.new.idents) <- ep.old.idents

ep <- RenameIdents(ep, ep.new.idents) # Rename the labels to fix the typo

# Create the CellChat object now
wounded.cc <- createCellChat(data = wounded.data.input, do.sparse = F)

# Add the meta-data from Seurat
wounded.cc <- addMeta(wounded.cc, meta = wounded.identity, meta.name = "labels") 
wounded.cc <- setIdent(wounded.cc, ident.use = "labels") # Set the labels to be the default cell identity

woundedGroupSize <- as.numeric(table(wounded.cc@idents)) # Get the number of cells in each group

wounded.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
wounded.cc <- subsetData(wounded.cc) # We subset the expression data of signalling genes to save on computational cost
wounded.cc <- identifyOverExpressedGenes(wounded.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
wounded.cc <- identifyOverExpressedInteractions(wounded.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
wounded.cc <- projectData(wounded.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

# We now infer the cell-cell communication network by calculating the communication probabilities
wounded.cc <- computeCommunProb(wounded.cc, raw.use = TRUE, population.size = FALSE) 
wounded.cc <- computeCommunProbPathway(wounded.cc) # Calculate the probabilities at the signalling level
wounded.cc <- aggregateNet(wounded.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

vertex.receiver <- c(5,6,7) # Focus on the fibroblasts
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(wounded.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, vertex.size = woundedGroupSize, pt.title=32)
netAnalysis_contribution(wounded.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
wounded.cc <- netAnalysis_signalingRole(wounded.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(wounded.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
wounded.cc <- identifyCommunicationPatterns(wounded.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(wounded.cc, pattern = "outgoing") 

# Visualize the communication pattern using dot plot
netAnalysis_dot(wounded.cc, pattern = "outgoing")

# Now let's look at incoming patterns
nPatterns <- 5
wounded.cc <- identifyCommunicationPatterns(wounded.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(wounded.cc, pattern = "incoming")

netAnalysis_dot(wounded.cc, pattern = "incoming")

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
wounded.cc <- computeNetSimilarity(wounded.cc, type = "functional", thresh = 0.25)
wounded.cc <- netEmbedding(wounded.cc, type = "functional")
wounded.cc <- netClustering(wounded.cc, type = "functional", k = 5)
netVisual_embedding(wounded.cc, type = "functional", pathway.remove.show = F, label.size = 3.5)

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
wounded.cc <- computeNetSimilarity(wounded.cc, type = "structural", thresh = 0.25)
wounded.cc <- netEmbedding(wounded.cc, type = "structural")
wounded.cc <- netClustering(wounded.cc, type = "structural")
netVisual_embedding(wounded.cc, type = "structural", label.size = 3.5)

save(unwounded.cc, wounded.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020adjustedcellchat.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020cellchat.rdata")

# Merge the CellChat objects
merged.cc <- mergeCellChat(list(unwounded.cc.secreted, wounded.cc), add.names = c("Unwounded", "Wounded"))

save(merged.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020mergedadjustedcellchat.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020mergedadjustedcellchat.rdata")

# Let's analyse for structural similarity
merged.cc <- computeNetSimilarityPairwise(merged.cc, type = "structural")
merged.cc <- netEmbedding(merged.cc, type = "structural")
merged.cc <- netClustering(merged.cc, type = "structural")

# Visualise the structural similarity
netVisual_embeddingPairwise(merged.cc, type = "structural", dot.size = c(5, 15), label.size = 6) + theme(text=element_text(size=24), axis.text=element_text(size=24)) 
netVisual_embeddingPairwiseZoomIn(merged.cc, type = "structural", dot.size = c(3, 9), label.size = 4) + coord_fixed() + theme(text=element_text(size=24), axis.text=element_text(size=24)) 

# Visualise the pathway distance in the learnt joint manifold
rankSimilarity(merged.cc, type = "structural") + theme(text=element_text(size=24), axis.text=element_text(size=24))

# Identify and visualise the conserved and context-specified signalling pathways
rankNet(merged.cc, mode = "comparison")

### Run CellChat on the subsetted data

# Subset the wounded data for just fibroblasts and immune cells
ep.subsetted <- subset(ep, idents = c("Fibroblast I", "Fibroblast II", "Myofibroblast", "Macrophage I",
                                      "Macrophage II", "Macrophage III", "Dendritic/Langerhans cell", "T cell"))

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
wounded.subsetted.data.input <- GetAssayData(ep.subsetted, assay = "RNA", slot = "data")  # Normalised data matrix
# wounded.data.input <- ep@data
wounded.subsetted.labels <- Idents(ep.subsetted)  # Get the cell group labels
wounded.subsetted.identity <- data.frame(group = wounded.subsetted.labels, row.names = names(wounded.subsetted.labels)) # Dataframe of the cell labels
# wounded.identity <- data.frame(group = ep@ident, row.names = names(ep@ident))
unique(wounded.subsetted.identity$group)

# Create the CellChat object now
wounded.subsetted.cc <- createCellChat(data = wounded.subsetted.data.input, do.sparse = F)

# Add the meta-data from Seurat
wounded.subsetted.cc <- addMeta(wounded.subsetted.cc, meta = wounded.subsetted.identity, meta.name = "labels") 
wounded.subsetted.cc <- setIdent(wounded.subsetted.cc, ident.use = "labels") # Set the labels to be the default cell identity

woundedSubsettedGroupSize <- as.numeric(table(wounded.subsetted.cc@idents)) # Get the number of cells in each group

wounded.subsetted.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
wounded.subsetted.cc <- subsetData(wounded.subsetted.cc) # We subset the expression data of signalling genes to save on computational cost
wounded.subsetted.cc <- identifyOverExpressedGenes(wounded.subsetted.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
wounded.subsetted.cc <- identifyOverExpressedInteractions(wounded.subsetted.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
wounded.subsetted.cc <- projectData(wounded.subsetted.cc, PPI.mouse) # Other option includes PPI.human (optional step that Suoqin told us to maybe dismiss)

# We now infer the cell-cell communication network by calculating the communication probabilities
wounded.subsetted.cc <- computeCommunProb(wounded.subsetted.cc, raw.use = TRUE, population.size = FALSE) 
wounded.subsetted.cc <- computeCommunProbPathway(wounded.subsetted.cc) # Calculate the probabilities at the signalling level
wounded.subsetted.cc <- aggregateNet(wounded.subsetted.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

vertex.receiver <- c(1, 2, 3) # Focus on the fibroblasts
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(wounded.subsetted.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = woundedSubsettedGroupSize, pt.title=32)
netAnalysis_contribution(wounded.subsetted.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
wounded.subsetted.cc <- netAnalysis_signalingRole(wounded.subsetted.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(wounded.subsetted.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
wounded.subsetted.cc <- identifyCommunicationPatterns(wounded.subsetted.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(wounded.subsetted.cc, pattern = "outgoing") 

# Visualize the communication pattern using dot plot
netAnalysis_dot(wounded.subsetted.cc, pattern = "outgoing")

# Now let's look at incoming patterns
nPatterns <- 5
wounded.subsetted.cc <- identifyCommunicationPatterns(wounded.subsetted.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(wounded.subsetted.cc, pattern = "incoming")

netAnalysis_dot(wounded.subsetted.cc, pattern = "incoming")

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
wounded.subsetted.cc <- computeNetSimilarity(wounded.subsetted.cc, type = "functional", thresh = 0.25)
wounded.subsetted.cc <- netEmbedding(wounded.subsetted.cc, type = "functional")
wounded.subsetted.cc <- netClustering(wounded.subsetted.cc, type = "functional", k = 4)
netVisual_embedding(wounded.subsetted.cc, type = "functional", pathway.remove.show = F, label.size = 3.5)

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
wounded.subsetted.cc <- computeNetSimilarity(wounded.subsetted.cc, type = "structural", thresh = 0.25)
wounded.subsetted.cc <- netEmbedding(wounded.subsetted.cc, type = "structural")
wounded.subsetted.cc <- netClustering(wounded.subsetted.cc, type = "structural")
netVisual_embedding(wounded.subsetted.cc, type = "structural", label.size = 3.5)

save(wounded.subsetted.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020woundedsubcellchat.rdata")

### We now look at the large-wounded data (PWD12) from Guerrero-Juarez et al. (2019)
load(file = "~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/10X_wounded_PWD12.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
pwd12.data.input <- GetAssayData(ep, assay = "RNA", slot = "data")  # Normalised data matrix
pwd12.labels <- Idents(ep)  # Get the cell group labels
levels(pwd12.labels) <- c("Fibroblast I", "Fibroblast II", "Myeloid cell", "Fibroblast III", "Endothelial", "Fibroblast IV", "T lymphocyte", "B lymphocyte", "Fibroblast V", "Schwann cell", "Erythrocyte", "Dendritic cell", "Lymphatic endothelial cell")
pwd12.identity <- data.frame(group = pwd12.labels, row.names = names(pwd12.labels)) # Dataframe of the cell labels
unique(pwd12.identity$group)

# Create the CellChat object now
pwd12.cc <- createCellChat(data = pwd12.data.input, do.sparse = F)

# Add the meta-data from Seurat
pwd12.cc <- addMeta(pwd12.cc, meta = pwd12.identity, meta.name = "labels") 
pwd12.cc <- setIdent(pwd12.cc, ident.use = "labels") # Set the labels to be the default cell identity

pwd12GroupSize <- as.numeric(table(pwd12.cc@idents)) # Get the number of cells in each group

pwd12.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
pwd12.cc <- subsetData(pwd12.cc) # We subset the expression data of signalling genes to save on computational cost
pwd12.cc <- identifyOverExpressedGenes(pwd12.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
pwd12.cc <- identifyOverExpressedInteractions(pwd12.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
pwd12.cc <- projectData(pwd12.cc, PPI.mouse) # Other options includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
pwd12.cc <- computeCommunProb(pwd12.cc) 
pwd12.cc <- computeCommunProbPathway(pwd12.cc) # Calculate the probabilities at the signalling level
pwd12.cc <- aggregateNet(pwd12.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

vertex.receiver <- c(1, 2, 4, 6, 9) # Focus on the epidermal and HF cells on the left
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(pwd12.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout="circle", vertex.size = pwd12GroupSize)
netAnalysis_contribution(pwd12.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
pwd12.cc <- netAnalysis_signalingRole(pwd12.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(pwd12.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 4
pwd12.cc <- identifyCommunicationPatterns(pwd12.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(pwd12.cc, pattern = "outgoing") 

# Visualize the communication pattern using dot plot
netAnalysis_dot(pwd12.cc, pattern = "outgoing")

# Now let's look at incoming patterns
nPatterns <- 5
pwd12.cc <- identifyCommunicationPatterns(pwd12.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(pwd12.cc, pattern = "incoming")

netAnalysis_dot(pwd12.cc, pattern = "incoming")

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
pwd12.cc <- computeNetSimilarity(pwd12.cc, type = "functional", thresh = 0.25)
pwd12.cc <- netEmbedding(pwd12.cc, type = "functional")
pwd12.cc <- netClustering(pwd12.cc, type = "functional", k = 4)
netVisual_embedding(pwd12.cc, type = "functional", pathway.remove.show = F, label.size = 3.5)

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
pwd12.cc <- computeNetSimilarity(pwd12.cc, type = "structural", thresh = 0.25)
pwd12.cc <- netEmbedding(pwd12.cc, type = "structural")
pwd12.cc <- netClustering(pwd12.cc, type = "structural")
netVisual_embedding(pwd12.cc, type = "structural", label.size = 3.5)

save(pwd12.cc, file = "~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/guerrerojuarez2019cellchat.rdata")
load("~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/guerrerojuarez2019cellchat.rdata")

