### Application of NicheNet, which infers ligand-target-gene interactions, to the skin and wound healing data from Haensel et al. (2020) and Guerrero-Juarez (2019)
### Note that at the moment, NicheNet only has a database for human ligand-target-gene interactions, and it kind of makes an assumption that the mouse genes are in
### one-to-one correspondence. A lot of the steps are based on the following NicheNet tutorial https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md

# Load the relevant packages
library(nichenetr)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(Seurat)

# Load the relevant NicheNet databases
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")) # Prescribed database of ligand-receptor intearctions. I wonder if I can use CellChat's DB and combine it with NicheNet's ligand-target-gene matrix?
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# We need to also convert the human genes to mouse genes, as that's the model system we'll be analysing.
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

### Let's first look at the unwounded data from Haensel et al. (2020)
load(file = "~/Documents/CellCellCommunicationModelling/Data/Haensel2020/bs_cca.rdata") # Load the data, which includes a Seurat V2 object
ep <- UpdateSeuratObject(object = ep) # Update this object from Seurat V2 to a V3 object

cell_types = levels(Idents(ep)) # Isolate the cell types, so we can define the sender/receiver cells

# Define the sender and receiver cells and get their relevant expressed genes. Here, let's try epidermal-to-fibroblast signalling
sender_cell_types = c("Epidermal basal I", "Epidermal basal II", "Epidermal basal III", "Proliferative epidermal basal")  # Probably don't need spinous.\
# sender_cell_types = c("Fibroblast I", "Fibroblast II")
list_expressed_genes_sender = sender_cell_types %>% unique() %>% lapply(get_expressed_genes, ep, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# receiver_cells = c("Fibroblast I", "Fibroblast II", "Myofibroblast") # Not sure whether we need to incorporate myofibroblasts, but we'll see
receiver_cells = cell_types[10:15]  # Probably don't need spinous.
expressed_genes_receiver = get_expressed_genes(receiver_cells, ep, pct = 0.10) 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#### We now define the gene set of interest. In general, this may depend on the experimental condition in question. 
### For instance! unwounded or wounding, and seeing how targeting of fibroblasts during wounding changes could be interesting,
## but for here, I think the gene set of interest iwll be the same as the receiver cells. Maybe this means we change how the background set goes.
ep_receiver = subset(ep, idents = receiver_cells)

DE_table_receiver = FindMarkers(object = ep_receiver, ident.1 = c("Fibroblast I", "Fibroblast II"), min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Now define the potential ligands and receptors, which are those corresponding to the by the expressed genes from the sender and receiver populations.
ligands = lr_network %>% pull(from) %>% unique() # Ligands defined by the 'from' variable
receptors = lr_network %>% pull(to) %>% unique() # Receptors defined by the 'to' variable

expressed_ligands = intersect(ligands, expressed_genes_sender) # Weird formatting issues. NicheNet wrote them as all capitals while in this dataset they're not.
expressed_receptors = intersect(receptors, expressed_genes_receiver) # Weird formatting issues

# The potential ligands are defined as those in the database that correspond to those with the expressed genes in the sender and receiver groups.
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique() 

# We're now ready to analyse the ligand activity with NicheNet!n Again, note that we had to format the geneset to be all caps
# rownames(ligand_target_matrix) = toupper(rownames(ligand_target_matrix))
# colnames(ligand_target_matrix) = toupper(colnames(ligand_target_matrix))
ligand_activities = predict_ligand_activities(geneset = toupper(geneset_oi), background_expressed_genes = toupper(background_expressed_genes), ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

# We use the top 20 ligands to construct the ligand-receptor network
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# Need to reformat these for the dotplot
best_upstream_ligands_formatted <- best_upstream_ligands

# Just in case the genes are in all caps, this reformats everything
for (i in seq(1, length(best_upstream_ligands)))
{
  gene <- best_upstream_ligands[i]
  gene_split <- unlist(strsplit(gene, ""))
  gene_lowercase <- paste(tolower(gene_split[2:length(gene_split)]), collapse="")
  best_upstream_ligands_formatted[i] <- paste(gene_split[1], gene_lowercase, sep="")
}

DotPlot(ep, features = best_upstream_ligands_formatted %>% rev(), cols = "RdYlBu") + RotatedAxis()

# Let's now infer the target genes involved
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 25) %>% bind_rows() %>% drop_na() # I changed n here because 200 was too high

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network # Plot hte ligand-target network
