### Application of Suoqin Jin's CellChat package to various datas on unwounded, large wound centre and periphery, and small wound
### murine skin tissue data from Abbasi et al. (2020) (https://doi.org/10.1016/j.stem.2020.07.008).
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

# Load the data
load("~/Documents/CellCellCommunicationModelling/Data/Abbasi2020/abbasi2020integrated.rdata")

# Subset the data by various experimental conditions

