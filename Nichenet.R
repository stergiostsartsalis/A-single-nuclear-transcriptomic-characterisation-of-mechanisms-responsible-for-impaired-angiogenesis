setwd("~/VASC/GerritsSmith/NicheNet_wPVMAstro/")

library(nichenetr)
library(Seurat)
library(tidyverse)





seuratObj<-qs::qread( '../vasc_PVM_Astro_NicheNet.qs')


############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ BASIC VIGNETTE############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 

# Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))



#Because the expression data is of mouse origin, we will convert the NicheNet network gene symbols from human to mouse based on one-to-one orthology:
# lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
# rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
# 
# ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
# 
# weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()


# Define a “sender/niche” cell population and a “receiver/target” cell population 
#present in your expression data and determine which genes are expressed in both 
#populations


## receiver
receiver = "EC"

seuratObj<-SetIdent(seuratObj,value = 'cluster_celltype')

expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.05)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


## sender
sender_celltypes = c('EC',"FB",'PC','PVM','Astro')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


### define gene sets of interest

genesets<-readRDS("../../DEGs/DEGs_categorical_gerritssmithoxford/DEGs_endo_enrichment.RDS")
geneset_oi<-c(genesets$DOWN,genesets$UP)

hubs_enrichment<-readRDS("~/VASC/GerritsSmith/MEGENA_endothelial_ALL/hubs_enrichment.RDS")
geneset_oi<-hubs_enrichment$c1_21


# genesets<-readRDS("/Volumes/stsartsa/home/VASC/GerritsSmith/MEGENA_endothelial_ALL/MEGENA.output.RDS")
# DEGs_DOWN_enrichment_modules_endo_all <- read.delim("/Volumes/stsartsa/home/VASC/GerritsSmith/MEGENA_endothelial_ALL/DEGs_DOWN_enrichment_modules_endo_all.txt")
# DEGs_UP_enrichment_modules_endo_all <- read.delim("/Volumes/stsartsa/home/VASC/GerritsSmith/MEGENA_endothelial_ALL/DEGs_UP_enrichment_modules_endo_all.txt")
# module.table <- read.delim("/Volumes/stsartsa/home/VASC/GerritsSmith/MEGENA_endothelial_ALL/module.table.txt")
# geneset_oi<-genesets$module.output$modules$c1_172

DEGs_fibro<-readRDS('../DEGs/original/DEGs_fibro.RDS')
DEGs_peric<-readRDS('../DEGs/original/DEGs_peric.RDS')

genesets<-readRDS("~/VASC/GerritsSmith/DEGs/DEGs_categorical_gerritssmithoxford/DEGs_endo_enrichment.RDS")
geneset_oi<-c('SPRED2','ITPR1','RASAL2','TEK','FGF2','DUSP16','SPTBN1')


#### choose either the strictest DEGs or a more permissive DEG list
#all_degs<-unique(c(DEGs_fibro$UP,DEGs_fibro$DOWN,DEGs_peric$UP,DEGs_peric$DOWN,geneset_oi))

all_degs<-readRDS("~/VASC/GerritsSmith/NicheNet_wPVMAstro/poential_deg_ligans.RDS")


## define set of potential ligands

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

####

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#


#Perform NicheNet ligand activity analysis: rank the potential ligands based 
#on the presence of their target genes in the gene set of interest (compared to the 
#background set of genes)

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

write.table(ligand_activities,'ligand_activites_EC_DEGs.txt',sep='\t',row.names = F)


#The number of top-ranked ligands that are further used to predict active target
#genes and construct an active ligand-receptor network is here 20.

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

## we add the DEG ligands in the "best"list
best_upstream_ligands<-unique(c(intersect(ligand_activities$test_ligand,all_degs),best_upstream_ligands))

#These ligands are expressed by one or more of the input sender cells.
#To see which cell population expresses which of these top-ranked ligands, you can run the following:


d<-DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()


ggplot2::ggsave('dotplot_ligands_DEGs_all_lower_expr_threshold005.tiff',height = 4,width = 14)










#Infer receptors and top-predicted target genes of ligands that are top-ranked
#in the ligand activity analysis

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.5)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()


#Note that not all ligands from the top 20 are present in this ligand-target heatmap. 
#The left-out ligands are ligands that don’t have target genes with high enough 
#regulatory potential scores. Therefore, they did not survive the used cutoffs. 
#To include them, you can be less stringent in the used cutoffs.
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))


ggsave('predicted_target_genes_endo_DEGs_all_lowerthresholds.tiff',p_ligand_target_network,height = 7,width = 7)
saveRDS(p_ligand_target_network$data,'predicted_target_genes_endo_DEGs_all_lowerthresholds_regulatory_potential.RDS')
#Receptors of top-ranked ligands

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()


p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

ggsave('predicted_receptors_DEGs_all_lower_thresholds.tiff',p_ligand_receptor_network,height = 7,width = 7)

saveRDS(p_ligand_receptor_network$data,'predicted_receptors_DEGs_all_lower_thresholds_prior_interact_potential.RDS')

#Receptors of top-ranked ligands, but after considering only bona fide 
#ligand-receptor interactions documented in literature and publicly available databases

lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()


p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")

ggsave('predicted_receptors_strict_DEGs_all_lowerthresholds.tiff',p_ligand_receptor_network_strict,height = 7,width = 7)


# Add log fold change information of ligands from sender cells

#In some cases, it might be possible to also check upregulation of ligands in sender cells. 
#This can add a useful extra layer of information next to the ligand activities defined by NicheNet, 
#because you can assume that some of the ligands inducing DE in receiver cells, will be DE themselves 
#in the sender cells.

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seuratObj = seuratObj, condition_colname = "diagnosis", condition_oi ='AD', condition_reference = 'CTR', expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

qs::qsave(DE_table_all, 'DE_table_all_NicheNet_ECfromEC_PVM_FB_PC.qs')
# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc


# change colors a bit to make them more stand out
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
p_ligand_lfc



# Summary visualizations of the NicheNet analysis

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot











############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
########### Assess how well top-ranked ligands can predict a gene set of interest#########
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 

#Read in expression data of interacting cells
seuratObj =readRDS('../CellChat_wPVM/vasc_PVM_for_CellChat.RDS')


#Load the ligand-target model we want to use
 
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))


# Load the gene set of interest and background of genes

genesets<-readRDS("~/VASC/GerritsSmith/DEGs/DEGs_categorical_gerritssmithoxford/DEGs_endo_enrichment.RDS")
geneset_oi<-genesets$DOWN


### choose cell types
## receiver
receiver = "EC"

seuratObj<-SetIdent(seuratObj,value = 'cluster_celltype')

expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


## sender
sender_celltypes = c('EC',"FB",'PC','PVM')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()







#Perform NicheNet’s ligand activity analysis on the gene set of interest

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% pull(test_ligand)


######For the top 20 ligands, we will now build a multi-ligand model that uses all top-ranked ligands to predict 
#whether a gene belongs to the p-EMT program of not. This classification model will be trained via 
#cross-validation and returns a probability for every gene.
# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k = 3 # 3-fold
n = 20 # 2 rounds

gene_predictions_top20_list = seq(n) %>% lapply(assess_rf_class_probabilities, folds = k, geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix)

#Evaluate now how well the target gene probabilies accord to the gene set assignments

# get performance: auroc-aupr-pearson
target_prediction_performances_cv = gene_predictions_top20_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% bind_rows() %>% mutate(round=seq(1:nrow(.)))

#Evaluate now whether genes belonging to the gene set are more likely to be top-predicted. 
#We will look at the top 5% of predicted targets here.

# get performance: how many p-EMT genes and non-p-EMT-genes among top 5% predicted targets
target_prediction_performances_discrete_cv = gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>% bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(gene_predictions_top20_list), each = 2))



#What is the fraction of p-EMT genes that belongs to the top 5% predicted targets?


target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()

#What is the fraction of non-p-EMT genes that belongs to the top 5% predicted targets?
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()

# We see that the p-EMT genes are enriched in the top-predicted target genes. 
#To test this, we will now apply a Fisher’s exact test for every cross-validation round and report the average p-value.
target_prediction_performances_discrete_fisher = gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
target_prediction_performances_discrete_fisher %>% unlist() %>% mean()


#Finally, we will look at which p-EMT genes are well-predicted in every cross-validation round.

top_predicted_genes = seq(length(gene_predictions_top20_list)) %>% lapply(get_top_predicted_genes,gene_predictions_top20_list) %>% purrr::reduce(full_join, by = c("gene","true_target"))
View(top_predicted_genes %>% filter(true_target))









############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############## Differential NicheNet analysis between conditions of interest############## 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 

library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #
#Read in expression data of interacting cells
seuratObj =readRDS('../CellChat_wPVM/vasc_PVM_for_CellChat.RDS')


#For the Differential NicheNet, we need to compare at least 2 niches or conditions to each other. 
#In this case, the 2 niches are the pEMT-high-niche and the pEMT-low-niche. We will adapt the 
#names of the cell types based on their niche of origin.

seuratObj@meta.data$celltype_aggregate = paste(seuratObj@meta.data$cluster_celltype, seuratObj@meta.data$diagnosis,sep = "_") # user adaptation required on own dataset



seuratObj@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)

celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seuratObj = SetIdent(seuratObj, value = seuratObj[[celltype_id]])


#Read in the NicheNet ligand-receptor network and ligand-target matrix

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)



#Define the niches/microenvironments of interest

niches = list(
  "AD" = list(
    "sender" = c("EC_AD", "PVM_AD", "PC_AD", "FB_AD"),
    "receiver" = c("EC_AD")),
  "CTR" = list(
    "sender" = c("EC_CTR", "PVM_CTR", "PC_CTR", "FB_CTR"),
    "receiver" = c("EC_CTR"))
)

#Calculate differential expression between the niches

assay_oi = "RNA" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seuratObj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types

DE_receiver = calculate_niche_de(seuratObj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets

DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))


expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")


#Combine sender-receiver DE based on L-R pairs:
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)



####### spatial information



include_spatial_info_sender = F # if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
specificity_score_spatial = "lfc"
# this is how this should be defined if you don't have spatial info
# mock spatial info
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
} 
if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj  = seuratObj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  
}
## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"
if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj  = seuratObj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}















#Calculate ligand activities and infer active ligand-target links

lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = calculate_niche_de_targets(seuratObj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()


geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))






lfc_cutoff = 0.75 

specificity_score_targets = "min_lfc"

DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()




top_n_target = 250

niche_geneset_list = list(
  "AD" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "CTR" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
)

ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)



#Calculate (scaled) expression of ligands, receptors and targets across cell types of interest
#(log expression values and expression fractions)

features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot = suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)


exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))


# Expression fraction and receptor

exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 


#Prioritization of ligand-receptor and ligand-target links
prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)


prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)

# Visualization of the Differential NicheNet output

top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche





############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
#################### Single-cell NicheNet’s ligand activity analysis ##################### 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 

library(nichenetr)
library(tidyverse)

#
#Read in expression data of interacting cells
seuratObj =readRDS('../CellChat_wPVM/vasc_PVM_for_CellChat.RDS')


#Load the ligand-target model we want to use

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))


# Load the gene set of interest and background of genes

genesets<-readRDS("/Volumes/stsartsa/home/VASC/GerritsSmith/DEGs/DEGs_categorical_gerritssmithoxford/DEGs_endo_enrichment.RDS")
geneset_oi<-genesets$DOWN



### choose cell types
## receiver
receiver = "EC"

seuratObj<-SetIdent(seuratObj,value = 'cluster_celltype')

expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


## sender
sender_celltypes = c('EC',"FB",'PC','PVM')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#Load the ligand-target model we want to use

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

#Perform NicheNet’s single-cell ligand activity analysis

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network$from %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)
receptors = lr_network$to %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
head(potential_ligands)


# Perform NicheNet’s single-cell ligand activity analysis

expression_scaled = as.matrix(seuratObj) %>% .[background_expressed_genes,WhichCells(seuratObj,idents = receptors==ADP)] %>% scale_quantile()


#Perform NicheNet’s single-cell ligand activity analysis
ligand_activities = predict_single_cell_ligand_activities(cell_ids = malignant_hn5_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

#Ligand prioritization by regression analysis
cell_scores_tbl = tibble(cell = malignant_hn5_ids, score = expression_scaled[malignant_hn5_ids,"TGFBI"])


normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)



output_correlation_analysis = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)


#####
inner_join(cell_scores_tbl,normalized_ligand_activities) %>% ggplot(aes(score,TNC)) + geom_point() + geom_smooth(method = "lm")









############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############## Inferring ligand-to-target signaling pathsm############## 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 




weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))




ligands_all = "APOE" # this can be a list of multiple ligands if required
targets_all = c('ANGPT2','FGF2','HIF1A','VEGFC')

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

DiagrammeR::render_graph(graph_min_max, layout = "tree")





#####data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)









#######################
#######################
#######################
#######################
#######################
#######################
#######################
######## circular plot


ligands_to_show<-c('FGF2','ROBO1','ANGPT2','PSAP','ANGPT1','EDN3','VEGFA','APOE','GPNMB','JAM3','PIK3CB','HSPG2','TF','PECAM1','AGT','NCAM1','CDH5','ITGAM','TNC')
targets_to_show<-c('ANGPT2','CFLAR','HIF1A','KLF9','NRXN3','RAC1','SEMA3B','SLC38A2','SPRED2','SPTBN1','TRIP10','ITPR1','MAGI1','MECOM','MEF2C','NCOR2','PICALM','FES','GNA14','LDHA','TPM1','PIK3CB','TEK','CDLK1','TCEA2','FGF2','ABCA1','ATP9A','CSGALNACT1','DLGAP4','MACROD2','PHC2','RASAL2','SESN3','ACTN1','MAP4')


######## circular plot

df.net <- readRDS('~/VASC/GerritsSmith/NicheNet_wPVMAstro/EC_DEGs/predicted_target_genes_endo_DEGs_all_lowerthresholds_regulatory_potential.RDS')
colnames(df.net)[1:2]<-c('ligand','target')
df.net<-df.net[df.net$ligand%in%ligands_to_show,]
df.net<-df.net[df.net$target%in%targets_to_show,]
df.net<-df.net[df.net$score!=0,]


df.net$ligand<-paste(df.net$ligand,'ligand',sep = '_')
df.net$target<-paste(df.net$target,'target',sep = '_')


df.net<-df.net[order(df.net$target),]
df.net<-df.net[order(df.net$ligand),]

lig<-unique(as.character(df.net$ligand))
tar<-unique(as.character(df.net$target))
gene_order<-c(lig[order(lig)],tar[order(tar,decreasing = T)])
#gene_order<-gene_order[c(1:24,38:39,25:37,40:74)]
### normal names

real_names<-structure(unique(c(as.character(df.net$ligand),as.character(df.net$target))))
names(real_names)<-as.character(real_names)
real_names<-real_names[gene_order]
for (x in 1:length(real_names))
{real_names[x]<-strsplit(real_names[x], split = '_')[[1]][1]}
  


## group by celltype/DEG

LR<-structure(unique(c(as.character(df.net$ligand),as.character(df.net$target))))
names(LR)<-as.character(LR)
LR<-LR[gene_order]
LR[grep('ligand',LR)]<-'ligand'
LR[grep('target',LR)]<-'target'

col_LR = c("ligand" = "black", "target" = "grey")

### dotplot seurat for circular heatmap
#matr<-d$data ## dotplot
matr<-qs::qread('matrix_dotplot_ligands_DEGs_all_lower_expr_threshold005.qs')
matr<-matr[,c('id','features.plot','pct.exp')]

mat<-as.data.frame(tidyr::pivot_wider(matr,names_from = 'features.plot',values_from = 'pct.exp'))
rownames(mat)<-mat$id
mat$id<-NULL

mat<-mat[,real_names[grep('ligand',names(real_names))]]
m2<-matrix(0,nrow=nrow(mat),ncol = length(grep('target',names(real_names))))

mat<-cbind(mat,m2)

colnames(mat)<-names(real_names)

mat<-t(mat)

mat<-as.data.frame(mat)
mat<-log(mat+1)

###########
matr<-d$data
matr<-matr[,c('id','features.plot','avg.exp.scaled')]

mat2<-as.data.frame(tidyr::pivot_wider(matr,names_from = 'features.plot',values_from = 'avg.exp.scaled'))
rownames(mat2)<-mat2$id
mat2$id<-NULL

mat2<-mat2[,real_names[grep('ligand',names(real_names))]]
m2<-matrix(0,nrow=nrow(mat2),ncol = length(grep('target',names(real_names))))

mat2<-cbind(mat2,m2)

colnames(mat2)<-names(real_names)
mat2<-t(mat2)



############# target genes (DEG) annotation
DEG <- read.delim("../../../../MAST_DGE/de/diagnosis_glmer_mod1/GerritsSmithOxford_endo_ensg/GerritsSmithOxford_endo_ensg_all_diagnosis_cngeneson_pc_mito_brain_region_sex_random_effect_manifest_glmer.tsv")
DEG<-DEG[DEG$padj<=0.1,]
DEG<-DEG[DEG$gene %in% real_names[LR=='target'],]

## remove duplicate MAGI1
DEG<-DEG[DEG$ensembl_gene_id!='ENSG00000282956',]
rownames(DEG)<-DEG$gene
DEG<-DEG[real_names[LR=='target'],]

DEG_logfc<-structure(c(rep(0,length(LR[LR=='ligand'])),DEG$logFC))
names(DEG_logfc)<-names(real_names)
DEG_logfc<-as.data.frame(DEG_logfc)

DEG_padj<-structure(c(rep(0,length(LR[LR=='ligand'])),-log10(DEG$padj)))
names(DEG_padj)<-names(real_names)
DEG_padj<-as.data.frame(DEG_padj)




#####
library(circlize)
library(ComplexHeatmap)
library(gridBase)
source('../../../../Basic_R_scRNAseq/circos_heatmap_initialize.R')


#col_fun = colorRamp2(range(df.net$score), c("#FFEEEE", "#FF0000"), transparency = 0)
col_fun = colorRamp2(range(df.net$score), hcl_palette = 'Lajolla', transparency = 0)

lgd_links<-Legend(title='Regulatory\npotential',col_fun = col_fun,direction = 'horizontal')

# lgd_links = Legend(labels = c('NDC','AD','NDC+AD'), legend_gp = gpar(fill=c('blue','red','darkolivegreen1')),
#                    title_position = "topleft", title = 'Diagnosis specificity\n(link color)')
# lgd_links_deg<-Legend(labels = c('DOWN','UP'), legend_gp = gpar(fill=c('cadetblue2','darkorange2')),
#                       title_position = "topleft",title='DEG (AD-NDC)')
# lgd_links_celltype<-Legend(labels = c('EC','FB','PC'), legend_gp = gpar(fill=c('darkgoldenrod3','burlywood2','blanchedalmond')),
#                            title_position = "topleft",title='Cell type')

col_fun1 = colorRamp2(c(-0.76, 0, 1.36), c('blue','white','red'))
lgd_links1<-Legend(title='DGE logFC',col_fun = col_fun1,direction = 'horizontal')

col_fun2 = colorRamp2(c(0, 4, 5.458),  hcl_palette = 'Oslo')
lgd_links2<-Legend(title=expression(bold("DGE\n-log"[10]~"(padj)")), col_fun = col_fun2,direction = 'horizontal')

col_fun3 = colorRamp2(c(-1.11, 0, 1.78), c('green','white','orange'))
lgd_links3<-Legend(title='Scaled\naverage\nexpression',col_fun = col_fun3,direction = 'horizontal')




lgd_links_LR<-Legend(labels = c('ligand','receptor'), legend_gp = gpar(fill=c('black','grey')),
                     title_position = "topleft",title='L/R')

tiff('legends.tiff',height = 10,width = 10,units = 'in',res=300)

lgd<-packLegend(lgd_links,lgd_links1,lgd_links2,lgd_links3,lgd_links_LR)
packLegend(lgd_links,lgd_links1,lgd_links2,lgd_links3,lgd_links_LR)

#a(lgd, just = c("left"))
draw(lgd)
dev.off()

#pdf('test.pdf', height = 9, width = 10.5, dpi)
tiff('circular_plot2.tiff', height = 35,width = 35,units = 'in',res=300)



####### add legend to the plot

plot.new()
# circle_size = unit(1, "snpc") # snpc unit gives you a square region
# 
# pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
#                       just = c("left", "center")))
# par(omi = gridOMI(), new = TRUE)

 #circos.par(canvas.xlim=c(-1,0),canvas.ylim=c(0,1))
circos.par(gap.degree=c(rep(2,length(LR[LR=='ligand'])-1),20,rep(2,length(LR[LR=='target'])-1),20))


sector_width<-vector(length = length(LR))
names(sector_width)<-names(LR)

## this allows for the sector width to be proportional to the number of LR pairs
# for (i in names(LR)){
#   sector_width[i]<-length(grep(i,c(df.net$ligand,df.net$target)))
# }

## this allows for the sector width to be proportional to the total score of the LR pairs
 for (i in names(LR)){
   sector_width[i]<-sum(df.net$score[df.net$ligand==i | df.net$target==i])
 }

xlim<-as.data.frame(cbind(rep(0,nrow(mat)),sector_width))
rownames(xlim)<-rownames(mat)
#circos.initialize(sectors = names(LR),xlim = xlim)
circos.heatmap.initialize(mat,split=factor(names(LR),levels=names(LR)),cluster=F,cell_width = sector_width)

# chordDiagram(df.net, col = col_fun, annotationTrack = c('grid'),annotationTrackHeight = .03 ,preAllocateTracks = list(list(track.height = .6),list(track.height=.03)),scale = F, group = LR, grid.col = grid.col)
# chordDiagram(df.net, col = col_fun,  group = LR, grid.col = grid.col)

circos.track(track.index=1,ylim=c(0,1),bg.border=NA)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], labels =real_names[CELL_META$sector.index] ,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=5)}, bg.border = NA)

circos.heatmap(DEG_logfc, col = col_fun1,track.height = 0.06,bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,cluster = F,rownames.cex = 1)
circos.heatmap(DEG_padj, col = col_fun2,track.height = 0.06,bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,cluster = F,rownames.cex = 1)
circos.heatmap(mat2, col = col_fun3,track.height = 0.2,bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,cluster = F,rownames.cex = 1)


circos.heatmap(LR, col = col_LR, track.height = 0.03,cluster = F)



############## links
df<-df.net

df.net$p1_1<-NA
df.net$p1_2<-NA
df.net$p2_1<-NA
df.net$p2_2<-NA

for(ll in unique(df.net$ligand)) 
  {
  
  width.ll<-get.cell.meta.data('xlim',sector.index = ll)
  
  width.ll<-width.ll[2]-width.ll[1]
  width.unit<-width.ll/(sum(df.net$score[df.net$ligand==ll])*nrow(df.net[df.net$ligand==ll,]))
  
  n1=0
  n2=0
  for (ii in grep(T,df.net$ligand==ll))
    {
  n1<-n2
    n2=n2+df.net$score[ii]*width.unit*nrow(df.net[df.net$ligand==ll,])
    
    df.net$p1_1[ii]<-n1
    df.net$p1_2[ii]<-n2     
  }
}


for(tt in unique(df.net$target)) 
{
  
  width.tt<-get.cell.meta.data('xlim',sector.index = tt)
  
  width.tt<-width.tt[2]-width.tt[1]
  width.unit<-width.tt/(sum(df.net$score[df.net$target==tt])*nrow(df.net[df.net$target==tt,]))
  
  n1=width.tt
  n2=width.tt
  for (ii in grep(T,df.net$target==tt))
  {
    n2<-n1
    n1=n2-df.net$score[ii]*width.unit*nrow(df.net[df.net$target==tt,])
    
    df.net$p2_1[ii]<-n1
    df.net$p2_2[ii]<-n2     
  }
}


source('../../../../Basic_R_scRNAseq/circlize_links.R')

for (iii in 1:nrow(df.net)){
  row_from<-df.net$ligand[iii]
  row_to<-df.net$target[iii]
  col<-col_fun(df.net$score[iii])
  point1<-as.vector(c(df.net$p1_1[iii],df.net$p1_2[iii]))
  point2<-as.vector(c(df.net$p2_1[iii],df.net$p2_2[iii]))
                                                                                                                                                                                                 
  circos.link( sector.index1=row_from,
                      sector.index2=row_to,
                      col = col,
                      point1=point1,
                      point2=point2)
}





circos.track(track.index = 2, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 51) { # the last sector
    cn = 'logFC'
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)

circos.track(track.index = 3, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 51) { # the last sector
    cn = expression("-log"[10]~"padj")
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)

circos.track(track.index = 4, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 18) { # the last sector
    cn = rev(colnames(mat2))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                1:n - 0.5, cn,
                cex = 5, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)


circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()
circos.clear()






#upViewport()

#draw(lgd, x = circle_size, just = c("left"))


dev.off()


### for source

write.table(df.net,'dfnet_for_source.txt',sep='\t',row.names = F)
write.table(DEG_logfc,'deglogfc_for_source.txt',sep='\t',col.names = NA)
write.table(DEG_padj,'degpadj_for_source.txt',sep='\t', col.names =  NA)
write.table(mat2,'ligand_cell_expr_for_source.txt',sep='\t',col.names = NA )



