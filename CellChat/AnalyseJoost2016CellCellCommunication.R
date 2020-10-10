### Application of Suoqin Jin's CellChat package to various datas on intrafollicular epidermis, which also includes
### hair follicle cells from Joost et al. (2016) (https://doi.org/10.1016/j.cels.2016.08.010). What's interesting about this data
### is that we actually clustered and identified cell types using Scanpy for once, so we'll be importing the data slightly differently!

# Load the relevant packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(reticulate)

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

# Import the data
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/Scanpy/joost2016.h5ad")
# access normalized data matrix
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))

data.input <- data.input - min(data.input) # Need to make sure the data isn't negative (it is in this case)
# access meta data
meta.data <- py_to_r(ad_object$obs)
identity <- data.frame(group = meta.data$leiden, row.names = row.names(meta.data))

# Create the cellchat object
joost.cc <-createCellChat(data = data.input, do.sparse = F)
joost.cc <- addMeta(joost.cc, meta = identity, meta.name = "labels")
joost.cc <- setIdent(joost.cc, ident.use = "labels") # set "labels" as default cell identity
levels(joost.cc@idents) # show factor levels of the cell labels

### Let's first look at secreted signalling
joost.cc@DB <- CellChatDB.use # Set the database for the data

joost.cc <- subsetData(joost.cc) # We subset the expression data of signalling genes to save on computational cost
joost.cc <- identifyOverExpressedGenes(joost.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Scanpy has an effect on this)
joost.cc <- identifyOverExpressedInteractions(joost.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# joost.cc <- projectData(joost.cc, PPI.mouse) # Other option includes PPI.human, apparently we shouldn't do this any more

# We now infer the cell-cell communication network by calculating the communication probabilities
joost.cc <- computeCommunProb(joost.cc, raw.use = TRUE, population.size = FALSE) 
joost.cc <- computeCommunProbPathway(joost.cc) # Calculate the probabilities at the signalling level
joost.cc <- aggregateNet(joost.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

joostGroupSize <- as.numeric(table(joost.cc@idents)) # Get the number of cells in each group

vertex.receiver <- c(4, 5) # Focus on the epidermal basal cells
pathways.show <- "EGF" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(joost.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = joostGroupSize, pt.title = 18)
netAnalysis_contribution(joost.cc, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
joost.cc <- netAnalysis_signalingRole(joost.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(joost.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(joost.cc, pattern = "outgoing", font.size = 3.25, font.size.title = 14) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(joost.cc, pattern = "outgoing", dot.size = c(1.5, 4.5))

# Now let's look at incoming patterns
nPatterns <- 5
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(joost.cc, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(joost.cc, pattern = "incoming", dot.size = c(1.5, 4.5))

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
joost.cc <- computeNetSimilarity(joost.cc, type = "functional", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "functional")
joost.cc <- netClustering(joost.cc, type = "functional", k = 4)
netVisual_embedding(joost.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
joost.cc <- computeNetSimilarity(joost.cc, type = "structural", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "structural")
joost.cc <- netClustering(joost.cc, type = "structural")
netVisual_embedding(joost.cc, type = "structural", label.size = 3.5)

# We will save this cellchat object as the 'secreted' version, and will now look at cell-cell contact signalling
joost.cc.secreted <- joost.cc

### Let us look at cell-cell contact too
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # Other options include ECM-Receptor and Cell-Cell Contact
joost.cc@DB <- CellChatDB.use # Set the database for the data

# joost.cc <- subsetData(joost.cc) # We subset the expression data of signalling genes to save on computational cost
joost.cc <- identifyOverExpressedGenes(joost.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Scanpy has an effect on this)
joost.cc <- identifyOverExpressedInteractions(joost.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# joost.cc <- projectData(joost.cc, PPI.mouse) # Other option includes PPI.human, apparently we shouldn't do this any more

# We now infer the cell-cell communication network by calculating the communication probabilities
joost.cc <- computeCommunProb(joost.cc, raw.use = TRUE, population.size = FALSE) 
joost.cc <- computeCommunProbPathway(joost.cc) # Calculate the probabilities at the signalling level
joost.cc <- aggregateNet(joost.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

joostGroupSize <- as.numeric(table(joost.cc@idents)) # Get the number of cells in each group

vertex.receiver <- c(4, 5) # Focus on the epidermal basal cells
pathways.show <- "CDH" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(joost.cc.cellcellcontact, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = joostGroupSize, pt.title = 18)
netAnalysis_contribution(joost.cc.cellcellcontact, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
joost.cc <- netAnalysis_signalingRole(joost.cc., slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(joost.cc.cellcellcontact, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 5
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(joost.cc, pattern = "outgoing", font.size = 3.25, font.size.title = 14) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(joost.cc.cellcellcontact, pattern = "outgoing", dot.size = c(1.5, 4.5))

# Now let's look at incoming patterns
nPatterns <- 5
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(joost.cc.cellcellcontact, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(joost.cc.cellcellcontact, pattern = "incoming", dot.size = c(1.5, 4.5))

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
joost.cc <- computeNetSimilarity(joost.cc, type = "functional", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "functional")
joost.cc <- netClustering(joost.cc, type = "functional", k = 4)
netVisual_embedding(joost.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
joost.cc <- computeNetSimilarity(joost.cc, type = "structural", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "structural")
joost.cc <- netClustering(joost.cc, type = "structural")
netVisual_embedding(joost.cc, type = "structural", label.size = 3.5)

joost.cc.cellcellcontact <- joost.cc

### Finally, let us consider cell-ECM receptor signalling!
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # Other options include ECM-Receptor and Cell-Cell Contact
joost.cc@DB <- CellChatDB.use # Set the database for the data

joost.cc <- subsetData(joost.cc) # We subset the expression data of signalling genes to save on computational cost
joost.cc <- identifyOverExpressedGenes(joost.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Scanpy has an effect on this)
joost.cc <- identifyOverExpressedInteractions(joost.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# joost.cc <- projectData(joost.cc, PPI.mouse) # Other option includes PPI.human, apparently we shouldn't do this any more

# We now infer the cell-cell communication network by calculating the communication probabilities
joost.cc <- computeCommunProb(joost.cc, raw.use = TRUE, population.size = FALSE) 
joost.cc <- computeCommunProbPathway(joost.cc) # Calculate the probabilities at the signalling level
joost.cc <- aggregateNet(joost.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

joostGroupSize <- as.numeric(table(joost.cc@idents)) # Get the number of cells in each group

vertex.receiver <- c(4, 5) # Focus on the epidermal basal cells
pathways.show <- "EGF" # All possible pathways are stored in @netP$pathways
netVisual_aggregate(joost.cc, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "circle", vertex.size = joostGroupSize, pt.title = 18)
netAnalysis_contribution(joost.cc.secreted, signaling = pathways.show) # Show the primary contributors to the selected pathway

# Let's now look at the signalling roles of each cell group
joost.cc <- netAnalysis_signalingRole(joost.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(joost.cc, signaling = pathways.show, width = 12, height = 2.5, font.size = 10) # Visualise the roles

# Now we identify outgoing communication patterns
nPatterns <- 3
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(joost.cc.secreted, pattern = "outgoing", font.size = 3.25, font.size.title = 14) 

# Visualize the communication pattern using dot plot
netAnalysis_dot(joost.cc.secreted, pattern = "outgoing", dot.size = c(1.5, 4.5))

# Now let's look at incoming patterns
nPatterns <- 3
joost.cc <- identifyCommunicationPatterns(joost.cc, pattern = "incoming", k = nPatterns)

netAnalysis_river(joost.cc, pattern = "incoming", font.size = 3.25)

netAnalysis_dot(joost.cc, pattern = "incoming", dot.size = c(1.5, 4.5))

# We now identify signalling groups based on functional similarity, which is when the sender/receiver groups are similar
joost.cc <- computeNetSimilarity(joost.cc, type = "functional", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "functional")
joost.cc <- netClustering(joost.cc, type = "functional", k = 4)
netVisual_embedding(joost.cc, type = "functional", pathway.remove.show = F, label.size = 4.5, dot.size = c(4, 12))

# We now identify signalling groups based on structural simlarity, which is to consider the signalling network structure
joost.cc <- computeNetSimilarity(joost.cc, type = "structural", thresh = 0.25)
joost.cc <- netEmbedding(joost.cc, type = "structural")
joost.cc <- netClustering(joost.cc, type = "structural")
netVisual_embedding(joost.cc, type = "structural", label.size = 3.5)

joost.cc.ecm <- joost.cc

# Save the CellChat objects for later use.
save(joost.cc.secreted, joost.cc.cellcellcontact, joost.cc.ecm, file="~/Documents/CellCellCommunicationModelling/Data/Joost2016/joost2016cellchat.rds")
