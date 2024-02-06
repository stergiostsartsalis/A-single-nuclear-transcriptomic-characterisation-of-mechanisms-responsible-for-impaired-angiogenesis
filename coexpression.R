
library(Seurat)
library(ggplot2)
library(MEGENA)
library(DGCA)
library(effsize)
library(bc3net)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(AUCell)
library(dplyr)
library(muscat)


setwd("~/VASC/GerritsSmith/MEGENA_fibroblast_ALL")


vasc_ALL<-qs::qread('../GerritsSmith_vascALL_ensg_sce.qs')

sce.rowdata <- read.delim("~/VASC/GerritsSmith/final/SCE/final_sce/sce-rowdata.tsv")

rownames(vasc_ALL)<-as.character(sce.rowdata$gene)

fibro<-vasc_ALL[,colData(vasc_ALL)$cluster_celltype=='fibro']

fibro_pb<-make_pseudobulk(fibro, pseudobulk_ID = 'manifest', pb_columns = c('manifest','individual','brain_region','diagnosis','group','sex','dataset'))

saveRDS(fibro_pb,'fibro_pb.RDS')


fibro_ALL_mat<-as.data.frame(fibro_pb$sumDat)

fibro_ALL_mat = fibro_ALL_mat[apply(fibro_ALL_mat,1,function(x) sum(x==0))<ncol(fibro_ALL_mat)*0.5,]


##########

library(limma)

fibro_ALL_mat<-voom(fibro_ALL_mat)$E
##########

### check for batch effect



metadata<-fibro_pb$annot_pb
rownames(metadata)<-as.character(metadata$group_sample)
metadata<-metadata[colnames(fibro_ALL_mat),]


library(PCAtools)

pca<-pca(fibro_ALL_mat,metadata = metadata)
biplot(pca,lab = metadata$dataset)

library(sva)

model_mat<-model.matrix(~as.factor(diagnosis)+as.factor(sex), data = metadata)

fibro_ALL_mat<-ComBat(fibro_ALL_mat, batch = metadata$dataset, mod = model_mat)



pca<-pca(fibro_ALL_mat,metadata = metadata)
biplot(pca,lab = metadata$diagnosis)


###################################

 #fibro_ALL_mat<-as.data.frame(fibro_ALL@assays$RNA@data)
 #mean_expr<-rowMeans(fibro_ALL_mat)
# 
#fibro_ALL_mat_non_zero<-fibro_ALL_mat[mean_expr!=0,]
# 


fibro_ALL_mat_non_zero<-fibro_pb$sumDat


## first makeDesign matrix to apply the filtering in each study 
 # library(DGCA)
 # design<-unique(cbind(fibro_ALL$sample,fibro_ALL$Study))
 # rownames(design)<-as.character(design[,1])
 # design<-as.data.frame(design)
 # design$V1<-NULL
 # design<-makeDesign(design$V2)
# 
# 
# fibro_ALL_filtered<-filterGenes(fibro_ALL_mat, filterTypes = c('central'),filterCentralPercentile = .1,filterCentralType = 'median',filterDispersionType ='variance',  filterDispersionPercentile  =  .1, sequential = T)





write.table(fibro_ALL_mat, 'fibro_ALL_filt.txt',sep = '\t')























method = "pearson" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2



library(doParallel)

registerDoParallel(cores = 8)

#ijw <- calculate.correlation(fibro_ALL_filtered , doPerm = cor.perm,output.corTable = FALSE,output.permFDR = FALSE, doPar = T, num.cores = 8)
#saveRDS(ijw, file = 'ijw.RDS')
ijw<-readRDS('ijw.RDS')

PFN<-calculate.PFN(ijw[,1:3], doPar=TRUE, num.cores=48)
saveRDS(PFN, file = 'PFNpar.RDS')
g <- graph.data.frame(PFN,directed = FALSE)

MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = T,num.cores = 48,n.perm = hub.perm,
                           save.output = T)
saveRDS(MEGENA.output, file = 'MEGENA.output.RDS')



if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


# summarize results

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                       min.size = 10,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)
if (!is.null(annot.table))
{
  # update annotation to map to gene symbols
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}
print(head(summary.output$modules,2))
print(summary.output$module.table)

module.table <- summary.output$module.table
colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".
hierarchy.obj <- plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
                                       arrow.size = 0.03,node.label.color = "blue")

write.table(module.table,'module.table.txt',sep='\t', row.names=F)
##hierarchy plot with a vector according which we colour the modules on the plot
# Limma_quantification_variable_fibro_ALL_MEGENA_modules_blockSample_group_ADCON <- read.delim("~/VASC/GerritsSmith/MEGENA_fibroblast_ALL/Limma_quantification_variable_fibro_ALL_MEGENA_modules_blockSample_group_ADCON.txt", row.names=1)
# datacol<-as.data.frame(cbind(Limma_quantification_variable_fibro_ALL_MEGENA_modules_blockSample_group_ADCON$logFC,Limma_quantification_variable_fibro_ALL_MEGENA_modules_blockSample_group_ADCON$adj.P.Val))

enrichm_deg_UP<-read.delim("~/VASC/GerritsSmith/MEGENA_fibroblast_ALL//DEGs_UP_enrichment_modules_fibro_all.txt", row.names=1)




























#########
#########
#########
#########
######### enrichment in DEGs MAST


DEGs <- readRDS("~/VASC/GerritsSmith/DEGs/DEGs_fibro.RDS")


DEGs_UP<-DEGs$UP

DEGs_DOWN<-DEGs$DOWN
###############################################

DEGs_UP_enrichment_table<-matrix(nrow=length(modules), ncol = 6)

colnames(DEGs_UP_enrichment_table)<-c('Module_name', 'TermID', 'Genes', 'All', 'pval', 'padj')

input_genes<-rownames(fibro$sumDat)

t=0
for (mn in names(modules))
{ t=t+1
mod_genes<-modules[mn]
mod_genes<-unlist(mod_genes)

p<-enrichment(mod_genes, reference =input_genes, genesets = as.data.frame(DEGs_UP), adj = 'bonferroni')
DEGs_UP_enrichment_table[t,] <-c(mn,p$TermID,p$genes,p$all, p$pval, p$padj )         
}
DEGs_UP_enrichment_table[,'padj']<-p.adjust(DEGs_UP_enrichment_table[,'pval'], method = 'fdr')
write.table(DEGs_UP_enrichment_table, 'DEGs_UP_enrichment_modules_fibro_all.txt', sep = "\t", row.names = F)


##########

###################
###################

DEGs_DOWN_enrichment_table<-matrix(nrow=length(modules), ncol = 6)

colnames(DEGs_DOWN_enrichment_table)<-c('Module_name', 'TermID', 'Genes', 'All', 'pval', 'padj')

input_genes<-rownames(fibro_ALL_mat_non_zero)

t=0
for (mn in names(modules))
{ t=t+1
mod_genes<-modules[mn]
mod_genes<-unlist(mod_genes)

p<-enrichment(mod_genes, reference =input_genes, genesets = as.data.frame(DEGs_DOWN), adj = 'bonferroni')
DEGs_DOWN_enrichment_table[t,] <-c(mn,p$TermID,p$genes,p$all, p$pval, p$padj )         
}
DEGs_DOWN_enrichment_table[,'padj']<-p.adjust(DEGs_DOWN_enrichment_table[,'pval'], method = 'fdr')
write.table(DEGs_DOWN_enrichment_table, 'DEGs_DOWN_enrichment_modules_fibro_all.txt', sep = "\t", row.names = F)








##################### two-layer network of GWAS genes



source('~/Basic_R_scRNAseq/pathway_analysis_enrichr.r')
source('~/Basic_R_scRNAseq/dotplot_enrichr.r')
library(dplyr)
library(ggplot2)


PFN<-readRDS('PFN.RDS')

Important_gene_sets<-readRDS('~/AD/Gene_sets/Important_GeneSets.RDS')
GWAS<-Important_gene_sets[['GWAS_extended']]

GWAS<-intersect(GWAS, unique(c(as.character(PFN$row),as.character(PFN$col))))

two_layer_network<-list()

for (g in GWAS)
{
  
  immed_netw<-c(as.character(PFN$row[PFN$col==g]),as.character(PFN$col[PFN$row==g]))  
  two_layer_network[[g]]<-unique(c(as.character(PFN$row[PFN$col %in% immed_netw]),as.character(PFN$col[PFN$row %in% immed_netw]),immed_netw))
  
  rm(immed_netw)
  
  
  
}
saveRDS(two_layer_network,'two_layer_network.RDS')

enrichment_gwas_networks<-list()

for (e in names(two_layer_network))
  
{
  enrichment_gwas_networks[[e]]<-pathway_analysis_enrichr(two_layer_network[[e]])
  
  
  
  
}
saveRDS(enrichment_gwas_networks,'enrichment_gwas_networks.RDS')




### extract for cytoscape

GWAS_gene<-'CCN2'
immed_netw<-c(as.character(PFN$row[PFN$col==GWAS_gene]),as.character(PFN$col[PFN$row==GWAS_gene]))  
PFN_GWAS<-unique(rbind(PFN[PFN$col==GWAS_gene,],PFN[PFN$row==GWAS_gene,],PFN[PFN$col %in% immed_netw,],PFN[PFN$row %in% immed_netw,]))

write.table(PFN_GWAS, file = paste(GWAS_gene,'edgelist.txt',sep='_'), row.names = F, sep = '\t', quote = F)




########### enrichment of two layer networks in DEGs


#########
#########
#########
######### enrichment in DEGs MAST


DEGs<- unlist(readRDS("~/VASC/GerritsSmith/DEGs/DEGs_fibro.RDS"))

#DEGs_UP<-DEGs$UP

#DEGs_DOWN<-DEGs$DOWN
###############################################

DEGs_enrichment_table<-matrix(nrow=length(two_layer_network), ncol = 6)

colnames(DEGs_enrichment_table)<-c('GWAS_network', 'TermID', 'Genes', 'All', 'pval', 'padj')

input_genes<-rownames(fibro_ALL_mat_non_zero)

t=0
for (mn in names(two_layer_network))
{ t=t+1
mod_genes<-two_layer_network[mn]
mod_genes<-unlist(mod_genes)

p<-enrichment(mod_genes, reference =input_genes, genesets = as.data.frame(DEGs), adj = 'bonferroni')
DEGs_enrichment_table[t,] <-c(mn,p$TermID,p$genes,p$all, p$pval, p$padj )         
}
DEGs_enrichment_table[,'padj']<-p.adjust(DEGs_enrichment_table[,'pval'], method = 'fdr')
write.table(DEGs_enrichment_table, 'DEGs_enrichment_two_layer_network_fibro_all.txt', sep = "\t", row.names = F)











    





####### module enrichment table

enrichment_table<-list()

for (i in names(enrichment_modules))
  
{
  enrichment_table[[i]]<-rbind(enrichment_modules[[i]][["GO_Molecular_Function_2018"]],enrichment_modules[[i]][["GO_Cellular_Component_2018"]],enrichment_modules[[i]][["GO_Biological_Process_2018"]],enrichment_modules[[i]][["KEGG_2019_Human"]],enrichment_modules[[i]][["WikiPathways_2019_Human"]],enrichment_modules[[i]][["Reactome_2016"]])
  
  
  
}

enrichment_table<-dplyr::bind_rows(enrichment_table, .id = "module")

write.table(enrichment_table, 'enrichment_table.txt',sep='\t',row.names = F)








####### module enrichment table

enrichment_table<-list()

for (i in names(enrichment_gwas_networks))
  
{
  enrichment_table[[i]]<-rbind(enrichment_gwas_networks[[i]][["GO_Molecular_Function_2018"]],enrichment_gwas_networks[[i]][["GO_Cellular_Component_2018"]],enrichment_gwas_networks[[i]][["GO_Biological_Process_2018"]],enrichment_gwas_networks[[i]][["KEGG_2019_Human"]],enrichment_gwas_networks[[i]][["WikiPathways_2019_Human"]],enrichment_gwas_networks[[i]][["Reactome_2016"]])
  
  
  
}

enrichment_table<-dplyr::bind_rows(enrichment_table, .id = "GWAS gene")

write.table(enrichment_table, 'two_layer_network_enrichment_table.txt',sep='\t',row.names = F)


#two_layer_network<-stack(two_layer_network)

two_layer_network<-plyr::ldply(two_layer_network, rbind)

write.table(two_layer_network, 'two_layer_networks.txt',sep='\t', row.names=F)








enrichment_table <- read.delim("~/VASC/GerritsSmith/MEGENA_fibroblast_ALL//enrichment_table.txt")

df<-enrichment_table[c(91,949,1242,1765,108,1278,1279,647,1286,1820,1500,1552,1049,3,114,412,1843,1871),]

df<-enrichment_table[enrichment_table$geneset%in%df$geneset,]

rownames(df)<-NULL






df$module<-as.numeric(sapply(as.character(df$module),function(x)strsplit(x, split = '1_', fixed = T)[[1]][2]))


dt.p <- df %>%
  #dplyr::filter(database %in% c("GO", "Reactome")) %>%
  dplyr::select(c(description, module, odds_ratio)) %>%
  reshape2::dcast(description ~ module, value.var = "odds_ratio",fun.aggregate = sum) 


rownames(dt.p) <- dt.p$description
dt.p$description <- NULL
dt.p<-dt.p[as.character(unique(df$description)),]

colnames(dt.p)<-sapply(colnames(dt.p),function(x)paste('M',x,sep=''))


### pheatmap
dt.p[is.na(dt.p)]<-0

## prepare the matrix for the circular heatmap
dt.p<-as.data.frame(t(dt.p))


dt.p<-dt.p[intersect(rownames(DEG_enr),rownames(dt.p)),]

DEG_enr<-DEG_enr[intersect(rownames(DEG_enr),rownames(dt.p)),]
DEG_enr<--log10(DEG_enr)
DEG_enr[DEG_enr<1]<-0


## circular heatmap
mat1<-rbind(dt.p,dt.p,dt.p,dt.p)
split<-c(rep(1,nrow(dt.p)),rep(2,nrow(dt.p)),rep(3,nrow(dt.p)),rep(4,nrow(dt.p)))

tiff('circular_heatmap_modules.tiff',width = 85,height = 50,units = 'cm',res = 300)

circos.clear()



col_fun1 = colorRamp2(seq(0,20,by=0.5), c('#FFFFFFFF',rev(viridis::viridis(length(seq(0.5,20,by=0.5))))))

circos.par(start.degree=180,gap.after=c(90,0,0,40),canvas.xlim=c(-1,0),canvas.ylim=c(0,1))
circos.heatmap(mat1, split = split, col = col_fun1, track.height = 0.5,bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,rownames.side = 'outside',cluster = F,rownames.cex = 4)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = rev(colnames(mat1))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 3, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)

mat2<-rbind(DEG_enr,DEG_enr,DEG_enr,DEG_enr)

col_fun2 = colorRamp2(c(0, 8), c('#FFFFFF','#FF0000'))




circos.heatmap(mat2, col = col_fun2, track.height = 0.05,bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,cluster = F)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = rev(colnames(mat2))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 3, adj = c(0, 0.5), facing = "inside")
  }
}, bg.border = NA)


#circos.heatmap(dt.p, split=split,col = col_fun1, track.height = 0.4, 
#bg.border = NA, bg.lwd = 2, bg.lty = 2,show.sector.labels = F,rownames.side = 'outside',rownames.cex = 0.3)
circos.clear()
dev.off()

# circos.updatePlotRegion(sector.index = 1, bg.lwd=0)
# circos.updatePlotRegion(sector.index = 3, bg.lwd=0)
tiff('circular_heatmap_modules_legends.tiff',width = 25,height = 20,units = 'cm',res = 300)

lgd_1<-Legend(title='Enrichment\nodds ratio',col_fun = col_fun1,direction = 'horizontal')
lgd_2 = Legend(title = expression(bold("-log"[10]~"padj")), col_fun = col_fun2,direction = 'horizontal')

lgd<-packLegend(lgd_1,lgd_2)

draw(lgd)


dev.off()

### to separately create a plot without text and 90 degree split width, put much lower gap values (eg 10 in the first and last position and do not write text)



## table for source data
if(all(rownames(mat1)==rownames(mat2))==T){
  write.table(cbind(mat1,mat2),'module_enrichment_table_for_source_data.txt',sep = '\t',col.names = NA)
  
}







