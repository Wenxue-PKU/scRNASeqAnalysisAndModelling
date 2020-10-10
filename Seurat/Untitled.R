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

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact).
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