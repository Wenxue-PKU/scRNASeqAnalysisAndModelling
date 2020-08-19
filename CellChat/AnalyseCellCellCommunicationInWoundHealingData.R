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
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

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

unwounded.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
unwounded.cc <- subsetData(unwounded.cc) # We subset the expression data of signalling genes to save on computational cost
future::plan("multiprocess", workers = 4) # Run in parallel as well
unwounded.cc <- identifyOverExpressedGenes(unwounded.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
unwounded.cc <- identifyOverExpressedInteractions(unwounded.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
unwounded.cc <- projectData(unwounded.cc, PPI.mouse) # Other option includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
unwounded.cc <- computeCommunProb(unwounded.cc) 
unwounded.cc <- computeCommunProbPathway(unwounded.cc) # Calculate the probabilities at the signalling level
unwounded.cc <- aggregateNet(unwounded.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# Plot the communications 
vertex.receiver <- seq(1,9) # Focus on the epidermal and HF cells on the left
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(unwounded.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "circle", vertex.size = unwoundedGroupSize)
netAnalysis_contribution(unwounded.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
unwounded.cc <- netAnalysis_signalingRole(unwounded.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(unwounded.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns = 5 
unwounded.cc <- identifyCommunicationPatterns(unwounded.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(unwounded.cc, pattern = "outgoing") 

# Visualize the communication pattern using dot plot
netAnalysis_dot(unwounded.cc, pattern = "outgoing")

### We now look at the wounded data to compare the differences in cell-cell communications
load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/sw_cca.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
wounded.data.input <- GetAssayData(ep, assay = "RNA", slot = "data")  # Normalised data matrix
# wounded.data.input <- ep@data
wounded.labels <- Idents(ep)  # Get the cell group labels
wounded.identity <- data.frame(group = wounded.labels, row.names = names(wounded.labels)) # Dataframe of the cell labels
# wounded.identity <- data.frame(group = ep@ident, row.names = names(ep@ident))
unique(wounded.identity$group)

# Create the CellChat object now
wounded.cc <- createCellChat(data = wounded.data.input, do.sparse = F)

# Add the meta-data from Seurat
wounded.cc <- addMeta(wounded.cc, meta = wounded.identity, meta.name = "labels") 
wounded.cc <- setIdent(wounded.cc, ident.use = "labels") # Set the labels to be the default cell identity

woundedGroupSize <- as.numeric(table(wounded.cc@idents)) # Get the number of cells in each group

wounded.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
wounded.cc <- subsetData(wounded.cc) # We subset the expression data of signalling genes to save on computational cost
future::plan("multiprocess", workers = 4) # Run in parallel as well
wounded.cc <- identifyOverExpressedGenes(wounded.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
wounded.cc <- identifyOverExpressedInteractions(wounded.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
wounded.cc <- projectData(wounded.cc, PPI.mouse) # Other option includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
wounded.cc <- computeCommunProb(wounded.cc) 
wounded.cc <- computeCommunProbPathway(wounded.cc) # Calculate the probabilities at the signalling level
wounded.cc <- aggregateNet(wounded.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

vertex.receiver <- seq(1,4) # Focus on the epidermal and HF cells on the left
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
p2 <- netVisual_aggregate(wounded.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "circle", vertex.size = woundedGroupSize)
netAnalysis_contribution(wounded.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

save(unwounded.cc, wounded.cc, file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/haensel2020cellchat.rdata")

### We now look at the large-wounded data (PWD12) from Guerrero-Juarez et al. (2019)
load(file = "~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/10X_wounded_PWD12.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

# We need to get the normalised gene expression matrix and the cell group labels to create the CellChat object
pwd12.data.input <- GetAssayData(ep, assay = "RNA", slot = "data")  # Normalised data matrix
pwd12.labels <- Idents(ep)  # Get the cell group labels
levels(pwd12.labels) <- c("Fibroblast I", "Fibroblast II", "Myeloid cell", "Fibroblast III", "Endothelial", "Fibroblast IV", "T lymphocyte", "B lymphocyte", "Fibroblst V", "Schwann cell", "Erythrocyte", "Dendritic cell", "Lymphatic endothelial cell")
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
future::plan("multiprocess", workers = 4) # Run in parallel as well
pwd12.cc <- identifyOverExpressedGenes(pwd12.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
pwd12.cc <- identifyOverExpressedInteractions(pwd12.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
pwd12.cc <- projectData(pwd12.cc, PPI.mouse) # Other options includes PPI.human

# We now infer the cell-cell communication network by calculating the communication probabilities
pwd12.cc <- computeCommunProb(pwd12.cc) 
pwd12.cc <- computeCommunProbPathway(pwd12.cc) # Calculate the probabilities at the signalling level
pwd12.cc <- aggregateNet(pwd12.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

vertex.receiver <- c(1, 2, 4, 6, 9) # Focus on the epidermal and HF cells on the left
pathways.show <- "TGFb" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(pwd12.cc, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = pwd12GroupSize)
netAnalysis_contribution(pwd12.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

save(pwd12.cc, file = "~/Documents/CellCellCommunicationModelling/Data/GuerreroJuarez2019/guerrerojuarez2019cellchat.rdata")
