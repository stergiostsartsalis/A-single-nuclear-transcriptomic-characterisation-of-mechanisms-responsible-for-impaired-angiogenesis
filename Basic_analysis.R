library(SingleCellExperiment)
library(Seurat)
library(ggplot2)


setwd("~/VASC/GerritsSmith")

GerritsSmith<-qs::qread('GerritsSmith_seurat.qs')

DimPlot(GerritsSmith, reduction = "UMAP_Liger", group.by = "clusters",label = F)+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())


ggsave('dimplot_GerritsSmith_ALL.tiff',height=4,width=4)




#### dimplot with larger celltypes

a<-as.data.frame(as.character(GerritsSmith$clusters))
colnames(a)<-'clusters'
a$clusters_for_markers<-a$clusters
a$clusters_for_markers[a$clusters %in% c('1','2','9')]<-'Micro'
a$clusters_for_markers[a$clusters %in% c('3','4','5')]<-'Astro'
a$clusters_for_markers[a$clusters %in% c('6')]<-'Oligo'
a$clusters_for_markers[a$clusters %in% c('10','7','8')]<-'Vasc'

GerritsSmith$clusters_for_markers<-as.character(a$clusters_for_markers)
GerritsSmith<-SetIdent(GerritsSmith, value='clusters_for_markers')

DimPlot(GerritsSmith, reduction = "UMAP_Liger", group.by = "clusters_for_markers",label = F)+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())

ggsave('dimplot_GerritsSmith_ALL_Celltypes.tiff',height=4,width=4)

source_data<-as.data.frame(cbind(as.character(GerritsSmith@meta.data$barcode),as.character(GerritsSmith$clusters_for_markers),as.character(GerritsSmith$sex),as.character(GerritsSmith$brain_region),as.character(GerritsSmith$diagnosis),GerritsSmith@reductions$UMAP_Liger@cell.embeddings))
colnames(source_data)[1:5]<-c('barcode','cell type','sex','brain_region','diagnosis')
source_data$`cell type`[source_data$`cell type`==11]<-'Neurons1'
source_data$`cell type`[source_data$`cell type`==12]<-'Neurons2'
source_data$`cell type`[source_data$`cell type`==13]<-'Lymphocytes1'
source_data$`cell type`[source_data$`cell type`==14]<-'Lymphocytes2'

source_data$diagnosis[source_data$diagnosis=='CTR']<-'NDC'
source_data$brain_region[source_data$brain_region=='EC']<-'EntC'

write.table(source_data,file = 'source_1a_S1.txt',sep='\t',row.names = F)

##########################
GerritsSmith<-qs::qread('GerritsSmith_vascALL_seurat.qs')

sce.rowdata <- read.delim("~/VASC/GerritsSmith/final/SCE/final_sce/sce-rowdata.tsv")

rownames(GerritsSmith)<-as.character(sce.rowdata$gene)

GerritsSmith<-as.Seurat(GerritsSmith,data=NULL)
GerritsSmith<-NormalizeData(GerritsSmith)
GerritsSmith<-ScaleData(GerritsSmith, features =unique(c("FLT1", "CLDN5", "NOSTRIN", "EBF1", "IFI27", "MYL9", "VWF","CD163", "PTPRC", "PDGFRB", "CD44", "TMEM119",'PECAM1','ACTA2','KCNJ8','CD248','ANPEP','RGS5','ZIC1','BMX','VEGFC','NR2F2','CSPG4','ABCC5','PDLIM3','PDGFRB','CSPG4','ANPEP','RGS5','ABCC9','CD248','S1PR3','PDLIM3','ACTA2','TAGLN','MYLK','MYH11','SNCG','PERGL','PLN','PDGFRA','LUM','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','SPP1')))



DimPlot(GerritsSmith, reduction = "UMAP_LIGER", group.by = "clusters",label = F)+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())


ggsave('GerritsSmith_vascALL.tiff',height=4,width=4)



DimPlot(GerritsSmith, reduction = "UMAP_PCA", group.by = "manifest",label=T)+NoLegend()
DimPlot(GerritsSmith, reduction = "UMAP_Liger", group.by = "clusters",label=T)+NoLegend()


#general feature

FeaturePlot(object = GerritsSmith, features = c("CD74", "GFAP",  "PLP1", "GAD1", "GAD2",'RBFOX3','MIAT','MEG3', "FLT1",'PDGFRA','PDGFRB','ACTA2'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')

FeaturePlot(object = GerritsSmith, features = c("CD74", "GFAP", "MBP", "PLP1", "GAD1", "GAD2", "VCAN", "LPAR6"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')
#microglia
FeaturePlot(object = GerritsSmith, features = c("HLA-DRA", "CX3CR1", "C1QB", "CSF1R", "CD74", "LPAR6", "C3", "DOCK8", "SYK", "P2RY12"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')
#Astro
FeaturePlot(object = GerritsSmith, features = c("AQP4", "SLC1A2", "GFAP", "FGFR3", "SLC14A1", "TNC", "CD44"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')
#Inhib Neuron
FeaturePlot(object = GerritsSmith, features = c("GAD1", "GAD2", "SST", "PVALB", "VIP", "SV2C"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')
#Excite Neuron
FeaturePlot(object = GerritsSmith, features = c("MIAT", "MEG3", "FEZF2", "RBFOX3"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')
#Oligo
FeaturePlot(object = GerritsSmith, features = c("MOBP", "MBP", "PLP1", "SOX10"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 
#OPC
FeaturePlot(object = GerritsSmith, features = c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "SOX6", "SMOC1"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 
#Endo
FeaturePlot(object = GerritsSmith, features = c("FLT1", "CLDN5", "NOSTRIN", "EBF1", "IFI27", "MYL9", "VWF"), pt.size = 0.05,reduction = 'UMAP_Liger', slot='data') 
#peri vas macro & vascular
FeaturePlot(object = GerritsSmith, features = c("CD163", "PTPRC", "PDGFRB", "CD44", "TMEM119",'PECAM1'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 
FeaturePlot(object = GerritsSmith, features = c('ACTA2','KCNJ8','CD248','ANPEP','RGS5','ZIC1','BMX','VEGFC','NR2F2','CSPG4','ABCC5','PDLIM3'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 
#PVM
#peri vas macro & vascular
FeaturePlot(object = GerritsSmith, features = c("CD163", "DPYD",'NAMPT','F13A1','IQGAP2'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 


## pericytes from mouse paper (Amy)
FeaturePlot(object = GerritsSmith, features = c('PDGFRB','CSPG4','ANPEP','RGS5','ABCC9','CD248','S1PR3'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 

## SMC from mouse paper (Amy)
FeaturePlot(object = GerritsSmith, features = c('PDLIM3','ACTA2','TAGLN','MYLK','MYH11','SNCG','PERGL','PLN'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 

## FB from mouse paper (Amy)
FeaturePlot(object = GerritsSmith, features = c('PDGFRA','LUM','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','SPP1'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 

#B cells

FeaturePlot(object = GerritsSmith, features = c("MS4A1"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')

# T cells

FeaturePlot(object = GerritsSmith, features = c("IL7R", "CCR7", "S100A4", "CD8A"), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')

## top markers from scflow for cluster 14

FeaturePlot(object = GerritsSmith, features = c("TXK", "MCTP2", "SKAP1", "SAMD3",'CD247','CD4'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data')



dev.off()



FeaturePlot(GerritsSmith, features=c('FLT1','VWF','NOSTRIN','COL1A1','COL12A1','COL6A1','PDGFRB','ACTA2','RGS5'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 

ggsave('FeaturePlot_test.tiff',height = 9,width=9)


### gene expression for source

f<-FeaturePlot(GerritsSmith, features=c("HLA-DRA", "CX3CR1", "C1QB", "CD74", "P2RY12","AQP4","SLC1A2","GFAP","SLC14A1","PLP1","MOBP","GAD1","RBFOX3",'FLT1','VWF','NOSTRIN','COL1A1','COL12A1','COL6A1','PDGFRB','ACTA2','RGS5'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data') 
source<-cbind(as.data.frame(f[[1]][['data']]),as.data.frame(f[[2]][['data']]),as.data.frame(f[[3]][['data']]),as.data.frame(f[[4]][['data']]),as.data.frame(f[[5]][['data']]),as.data.frame(f[[6]][['data']]),as.data.frame(f[[7]][['data']]),as.data.frame(f[[8]][['data']]),as.data.frame(f[[9]][['data']]),as.data.frame(f[[10]][['data']]),as.data.frame(f[[11]][['data']]),as.data.frame(f[[12]][['data']]),as.data.frame(f[[13]][['data']]),as.data.frame(f[[14]][['data']]),as.data.frame(f[[15]][['data']]),as.data.frame(f[[16]][['data']]),as.data.frame(f[[17]][['data']]),as.data.frame(f[[18]][['data']]),as.data.frame(f[[19]][['data']]),as.data.frame(f[[20]][['data']]),as.data.frame(f[[21]][['data']]),as.data.frame(f[[22]][['data']]))
source<-source[,c(1:4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88)]
write.table(source,'feature_plot_source.txt',sep='\t',col.names = NA)


##### Feature plots unenriched

setwd("~/VASC/GerritsSmith")
library(Seurat)
library(ggplot2)

sce<-qs::qread('~/AD/Basic_data_Glia/unenriched/unenriched_AD/final_sce.qs')
sce<-sce[,colData(sce)$brain_region!='mTemp']
gc()

sce<-as.Seurat(sce, data=NULL)
gc()

sce<-NormalizeData(sce)

gc()

cluster_celltype<-as.data.frame(sce$subcluster_celltype)
colnames(cluster_celltype)<-'subcluster_celltype'
cluster_celltype$cluster_celltype<-as.character(cluster_celltype$subcluster_celltype)
cluster_celltype$cluster_celltype[grep('EN',cluster_celltype$subcluster_celltype,ignore.case = F)]<-'EN'
cluster_celltype$cluster_celltype[grep('IN',cluster_celltype$subcluster_celltype,ignore.case = F)]<-'IN'
cluster_celltype$cluster_celltype[cluster_celltype$subcluster_celltype=='Micro']<-'MGL'
cluster_celltype$cluster_celltype[cluster_celltype$subcluster_celltype=='Astro']<-'AST'
cluster_celltype$cluster_celltype[cluster_celltype$subcluster_celltype=='Oligo']<-'OLG'
cluster_celltype$cluster_celltype[cluster_celltype$subcluster_celltype=='Endo']<-'VASC'
cluster_celltype$cluster_celltype[cluster_celltype$subcluster_celltype=='Micro']<-'MGL'




sce$cluster_celltype<-as.character(cluster_celltype$cluster_celltype)


DimPlot(sce,reduction = 'UMAP_Liger',group.by = 'cluster_celltype',label = T)+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())

ggsave('dimplot_unenriched_ALL.tiff',height=4,width=4)

embeddings<-as.matrix(sce@reductions$UMAP_Liger@cell.embeddings)
colnames(embeddings)<-c('umapliger_1','umapliger_2')

sce[['umapliger']]<-Seurat::CreateDimReducObject(embeddings = embeddings,key = 'umapliger_')

qs::qsave(sce,'~/AD/Basic_data_Glia/unenriched/unenriched_AD/final_seurat_ECSSC.qs')


FeaturePlot(object = sce, features = c("CD74", "GFAP",  "PLP1", "PCDH15",'RBFOX3', "GAD2",'MIAT','MEG3', "FLT1",'RGS5','ACTA2','COL1A1'), pt.size = 0.05, reduction = 'umapliger', slot='data')

ggsave('FeaturePlot_unenriched.tiff',height = 9,width=12)



### source data unenriched
object<-qs::qread('~/AD/Basic_data_Glia/unenriched/unenriched_AD/final_seurat_ECSSC.qs')
source_data<-as.data.frame(cbind(as.character(object@meta.data$barcode),as.character(object$cluster_celltype),as.character(object$sex),as.character(object$brain_region),as.character(object$diagnosis),object@reductions$UMAP_Liger@cell.embeddings))
colnames(source_data)[1:5]<-c('barcode','cell type','sex','brain_region','diagnosis')


source_data$diagnosis[source_data$diagnosis=='Control']<-'NDC'
source_data$brain_region[source_data$brain_region=='EC']<-'EntC'

write.table(source_data,file = 'source_S8_S9.txt',sep='\t',row.names = F)


### gene expression for source

f<-FeaturePlot(object = object, features = c("CD74", "GFAP",  "PLP1", "PCDH15",'RBFOX3', "GAD2",'MIAT','MEG3', "FLT1",'RGS5','ACTA2','COL1A1'), pt.size = 0.05, reduction = 'umapliger', slot='data')
source<-cbind(as.data.frame(f[[1]][['data']]),as.data.frame(f[[2]][['data']]),as.data.frame(f[[3]][['data']]),as.data.frame(f[[4]][['data']]),as.data.frame(f[[5]][['data']]),as.data.frame(f[[6]][['data']]),as.data.frame(f[[7]][['data']]),as.data.frame(f[[8]][['data']]),as.data.frame(f[[9]][['data']]),as.data.frame(f[[10]][['data']]),as.data.frame(f[[11]][['data']]),as.data.frame(f[[12]][['data']]))
source<-source[,c(1:4,8,12,16,20,24,28,32,36,40,44,48)]
write.table(source,'feature_plot_source_S10.txt',sep='\t',col.names = NA)



### Heatmap

vasc_ALL<-qs::qread('./GerritsSmith_vascALL_seurat.qs')
reclustered<-qs::qread('./reclustered_seurat_symbol.qs')
pericytes<-colnames(subset(reclustered,subset=clusters %in% c('4','7','10','22')))
smc<-colnames(subset(reclustered,subset=clusters %in% c('19')))

endofb<-colnames(subset(vasc_ALL, subset=cluster_celltype %in% c('endo','fibro')))
### pericytes are defined in the reclustered object
pericytes<-setdiff(pericytes,endofb)
smc<-setdiff(smc,endofb)

a<-as.data.frame(vasc_ALL$cluster_celltype)
colnames(a)<-'old'
a$new<-NA
a[pericytes,'new']<-'PC'
a[smc,'new']<-'SMC'

a[a$old=='endo','new']<-'EC'
a[a$old=='fibro','new']<-'FB'
vasc_ALL[['new_celltype']]<-as.character(a$new)

vasc_ALL<-subset(vasc_ALL, subset=new_celltype %in% c('EC','FB','PC','SMC'))

vasc_ALL[['cluster_celltype']]<-vasc_ALL[['new_celltype']]

saveRDS(vasc_ALL,'vasc_ALL_new_celltype_for_heatmap.RDS')

gc()



#################




GerritsSmith <- readRDS("./vasc_ALL_new_celltype_for_heatmap.RDS")
GerritsSmith<-NormalizeData(GerritsSmith)
GerritsSmith<-ScaleData(GerritsSmith, features =unique(c("FLT1", "CLDN5", "NOSTRIN", "EBF1", "IFI27", "MYL9", "VWF","CD163", "PTPRC", "PDGFRB", "CD44", "TMEM119",'PECAM1','ACTA2','KCNJ8','CD248','ANPEP','RGS5','ZIC1','BMX','VEGFC','NR2F2','CSPG4','ABCC5','PDLIM3','PDGFRB','CSPG4','ANPEP','NES','RGS5','GRM8','ABCC9','CD248','S1PR3','PDLIM3','ACTA2','TAGLN','MYLK','MYH11','SNCG','PERGL','PLN','PDGFRA','LUM','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','SPP1')))

library(viridis)
DoHeatmap(subset(GerritsSmith, downsample=15000), features =c("FLT1", "CLDN5", "NOSTRIN", "IFI27", "VWF",'PECAM1','NES','PDGFRA','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','ACTA2','ANPEP','RGS5','ZIC1','NR2F2','CSPG4','PDGFRB','CSPG4','ANPEP','RGS5','ABCC9','S1PR3','ACTA2','TAGLN','MYLK','MYH11'),group.by = 'cluster_celltype',label=F, lines.width = 40, slot = 'scale.data')+
  theme(axis.text.y=element_text(size=20,face='bold'))+
  scale_fill_viridis()+
  theme(legend.text = element_text(size=20,face='bold'))

ggsave('Markers_vasc_Heatmap_average.tiff',height=10,width = 8)

source_data<-AverageExpression(GerritsSmith,slot = 'scale.data',features =c("FLT1", "CLDN5", "NOSTRIN", "IFI27", "VWF",'PECAM1','NES','PDGFRA','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','ACTA2','ANPEP','RGS5','ZIC1','NR2F2','CSPG4','PDGFRB','CSPG4','ANPEP','RGS5','ABCC9','S1PR3','ACTA2','TAGLN','MYLK','MYH11'),group.by = 'cluster_celltype')
write.table(source_data$RNA,file = 'source_data_1b.txt',sep = '\t',col.names = NA)
######### average expression heatmap

GerritsSmith<-SetIdent(GerritsSmith, value = 'cluster_celltype')
average<-AverageExpression(GerritsSmith, slot = 'data',return.seurat = T,features =unique(c("FLT1", "CLDN5", "NOSTRIN", "IFI27", "VWF",'PECAM1','PDGFRA','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','ACTA2','ANPEP','RGS5','GRM8','ZIC1','NR2F2','CSPG4','PDGFRB','CSPG4','ANPEP','NES','RGS5','ABCC9','S1PR3','ACTA2','TAGLN','MYLK','MYH11') ))
cluster<-as.data.frame(colnames(average))
colnames(cluster)<-'cluster'
rownames(cluster)<-cluster$cluster
cluster$cluster<-as.character(rownames(cluster))
average[['cluster']]<-as.character(cluster$cluster)
average$cluster<-factor(average$cluster,levels = c('EC','FB','PC','SMC'))

library(viridis)
DoHeatmap(average,c("FLT1", "CLDN5", "NOSTRIN", "IFI27", "VWF",'PECAM1','NES','PDGFRA','DCN','COL3A1','COL5A1','COL8A2','COL12A1','MMP2','COL6A1','COL1A1','ANPEP','RGS5','GRM8','ZIC1','NR2F2','CSPG4','PDGFRB','ANPEP','RGS5','ABCC9','S1PR3','ACTA2','TAGLN','MYLK','MYH11'),label=F,  slot = 'scale.data',group.by = 'cluster',group.bar = T, draw.lines = F)+
  theme(axis.text.y=element_text(size=20,face='bold'))+
  scale_fill_viridis()+
  theme(legend.text = element_text(size=15,face='bold'))+
  theme(legend.title = element_text(size=15, face='bold'))

ggsave('Markers_vasc_Heatmap_average.tiff',height=10,width = 8)



############## Feature plot reclustered
setwd('~/VASC/GerritsSmith/reclustered/')
vasc_ALL<-qs::qread('../GerritsSmith_vascALL_seurat.qs')
reclustered<-qs::qread('./reclustered_seurat_symbol.qs')
pericytes<-colnames(subset(reclustered,subset=clusters %in% c('4','7','10','22')))
smc<-colnames(subset(reclustered,subset=clusters %in% c('19')))
endo<-colnames(subset(vasc_ALL,subset=cluster_celltype %in% 'endo'))
fibro<-colnames(subset(vasc_ALL,subset=cluster_celltype %in% 'fibro'))

a<-as.data.frame(reclustered$cluster_celltype)
colnames(a)<-'old'
a$new<-NA
a[pericytes,'new']<-'PC'
a[smc,'new']<-'SMC'

a[endo,'new']<-'EC'
a[fibro,'new']<-'FB'
reclustered[['new_celltype']]<-as.character(a$new)

reclustered<-subset(reclustered, subset=new_celltype %in% c('EC','FB','PC','SMC'))

reclustered[['cluster_celltype']]<-reclustered[['new_celltype']]


FeaturePlot(reclustered, features=c('FLT1','VWF','NOSTRIN','COL1A1','COL12A1','COL6A1','PDGFRB','ACTA2','RGS5'), pt.size = 0.05, reduction = 'UMAP_Liger', slot='data',raster = F) 

ggsave('FeaturePlot_vasc_markers_reclustered.tiff',height = 9,width=9)




######### stacked violin plot
source('~/Basic_R_scRNAseq/StackedVlnPlot.R')

GerritsSmith<-SetIdent(GerritsSmith,value='cluster_celltype')

#Gerrits$cluster_celltype<-factor(GerritsSmith$cluster_celltype,levels=)

pdf('stackedViolin_vasc.pdf',height = 18,width = 6)
print(StackedVlnPlot(GerritsSmith, features=c('FLT1','VWF','NOSTRIN','COL1A1','COL12A1','PDGFRB','ACTA2'),group.by='cluster_celltype',pt.size = 1))

dev.off()









########## markers comparison with mouse atlas paper
GerritsSmith<-SetIdent(GerritsSmith, value='cluster_celltype')
markersEC<-FindMarkers(GerritsSmith, ident.1 = 'endo',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')


markersFB<-FindMarkers(GerritsSmith, ident.1 = 'fibro',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')



markers<-list()
markers[['EC']]<-markersEC
markers[['FB']]<-markersFB

reclustered<-SetIdent(reclustered,value='cluster_celltype')

markers[['PC']]<-FindMarkers(reclustered, ident.1 = 'PC',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')

markers[['SMC']]<-FindMarkers(reclustered, ident.1 = 'SMC',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')
saveRDS(markers, 'markers_EC_FB.RDS')


### all gerritssmith markers
GerritsSmith<-qs::qread('GerritsSmith_seurat.qs')
a<-as.data.frame(as.character(GerritsSmith$clusters))
colnames(a)<-'clusters'
a$clusters_for_markers<-NA
a$clusters_for_markers[a$clusters %in% c('1','2','9')]<-'Micro'
a$clusters_for_markers[a$clusters %in% c('3','4','5')]<-'Astro'
a$clusters_for_markers[a$clusters %in% c('6')]<-'Oligo'
a$clusters_for_markers[a$clusters %in% c('10','7','8')]<-'Vasc'

GerritsSmith$clusters_for_markers<-as.character(a$clusters_for_markers)
GerritsSmith<-SetIdent(GerritsSmith, value='clusters_for_markers')

markers<-FindAllMarkers(GerritsSmith, only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')
saveRDS(markers,'GerritsSmith_allMarkers_wholeDataset_byCelltype_and_small_unknown_clusters.RDS')





EC_mouse_markers<-read.delim('~/Gene_sets/Endothelial_specific_markers_mouse_atlas_paper.txt',header=T)
FB_mouse_markers<-read.delim('~/Gene_sets/FIbroblast_specific_markers_mouse_atlas_paper.txt',header=T)
PC_mouse_markers<-read.delim('~/Gene_sets/PC_specific_markers_mouse_atlas_paper.txt',header=T)
SMC_mouse_markers<-read.delim('~/Gene_sets/SMC_specific_markers_mouse_atlas_paper.txt',header=T)



mouse_vasc_markers<-list()
mouse_vasc_markers[['EC']]<-as.character(EC_mouse_markers$Gene.Symbol[1:100])
mouse_vasc_markers[['FB']]<-as.character(FB_mouse_markers$Gene.Symbol[1:100])
mouse_vasc_markers[['PC']]<-as.character(PC_mouse_markers$Gene.Symbol[1:100])
mouse_vasc_markers[['SMC']]<-as.character(SMC_mouse_markers$Gene.Symbol[1:100])

saveRDS(mouse_vasc_markers,'~/Gene_sets/mouse_vasc_markers.RDS')

length(intersect(rownames(markers$EC),PCSMC_mouse_markers$Gene.Symbol))

length(intersect(PCSMC_mouse_markers$Gene.Symbol,FB_mouse_markers$Gene.Symbol))

#### enrichment with bc3net
library(bc3net)
all_genes<-rownames(GerritsSmith)

enrichment_EC_markers<-enrichment(rownames(markers$EC), reference =all_genes, genesets = mouse_vasc_markers, adj = 'bonferroni')

enrichment_FB_markers<-enrichment(rownames(markers$FB), reference =all_genes, genesets = mouse_vasc_markers, adj = 'bonferroni')

enrichment_PC_markers<-enrichment(rownames(markers$PC), reference =all_genes, genesets = mouse_vasc_markers, adj = 'bonferroni')

enrichment_SMC_markers<-enrichment(rownames(markers$SMC), reference =all_genes, genesets = mouse_vasc_markers, adj = 'bonferroni')

write.table(rbind(enrichment_EC_markers,enrichment_FB_markers,enrichment_PC_markers,enrichment_SMC_markers),'enrichment_markers_human_mouse_vasc.txt',sep='\t',row.names = F)



### dotplot enrichment markers

enrichment_markers_human_mouse_vasc <- read.delim("~/VASC/GerritsSmith/enrichment_markers_human_mouse_vasc.txt")

enrichment_markers_human_mouse_vasc$Cluster_markers<-'NA'

enrichment_markers_human_mouse_vasc$Cluster_markers[1:4]<-'EC'
enrichment_markers_human_mouse_vasc$Cluster_markers[5:8]<-'FB'
enrichment_markers_human_mouse_vasc$Cluster_markers[9:12]<-'PC'
enrichment_markers_human_mouse_vasc$Cluster_markers[13:16]<-'SMC'


####
library(ggplot2)
ggplot(enrichment_markers_human_mouse_vasc, aes(x=TermID,y=Cluster_markers))+geom_point(aes(colour=padj,size=genes))+
  scale_size(name = "size", range = c(0, 8)) +
  scale_fill_gradient(
    low = "navy", high = "gold", name = "padj",
    guide = guide_colorbar(reverse = F),
    aesthetics = c("colour"))+
  xlab('Cluster markers (from  )')+
  ylab('Cluster markers\n(present study)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 25),axis.text.y = element_text(size = 25))+
  theme(axis.title = element_text(size = 20))+ ## larger labels for paper: 30
  scale_x_discrete(labels=c('EC','FB','PC','SMC'))+
  theme(legend.text = element_text(size=15), legend.title = element_text(size=15))

ggsave('Dotplot_enrichment_mouse_markers.tiff',dpi = 300,compression='lzw',width = 5,height = 6)







############# HUMAN markers Garcia

Human_vasc_markers <- read_excel("~/Gene_sets/Human_vasc_markers.xlsx", 
                                 sheet = "Post Mortem Vascular Subcluster", 
                                 skip = 1)


human_vasc_markers<-list()
human_vasc_markers[['EC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('aEndo','capEndo','vEndo')]))
human_vasc_markers[['FB']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Fib1','Fib2','Fib3')]))

human_vasc_markers[['PC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Per1','Per2','Per3','Per4')]))
human_vasc_markers[['SMC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('aSMC','vSMC')]))


human_vasc_markers[['EC_arterial']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('aEndo')]))
human_vasc_markers[['EC_capillary']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('capEndo')]))
human_vasc_markers[['EC_venous']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('vEndo')]))

human_vasc_markers[['PC1']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Per1')]))
human_vasc_markers[['PC2']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Per2')]))
human_vasc_markers[['PC3']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Per3')]))
human_vasc_markers[['PC4']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('Per4')]))

human_vasc_markers[['SMC_arterial']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('aSMC')]))
human_vasc_markers[['SMC_venous']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c('vSMC')]))



saveRDS(human_vasc_markers,'human_vasc_markers.RDS')



### only top 100 markers

human_vasc_markers$EC<-human_vasc_markers$EC[1:100]
human_vasc_markers$FB<-human_vasc_markers$FB[1:100]
human_vasc_markers$PC<-human_vasc_markers$PC[1:100]
human_vasc_markers$SMC<-human_vasc_markers$SMC[1:100]

human_vasc_markers$EC_arterial<-human_vasc_markers$EC_arterial[1:100]
human_vasc_markers$EC_capillary<-human_vasc_markers$EC_capillary[1:100]
human_vasc_markers$EC_venous<-human_vasc_markers$EC_venous[1:100]
human_vasc_markers$PC1<-human_vasc_markers$PC1[1:100]
human_vasc_markers$PC2<-human_vasc_markers$PC2[1:100]
human_vasc_markers$PC3<-human_vasc_markers$PC3[1:100]
human_vasc_markers$PC4<-human_vasc_markers$PC4[1:100]
human_vasc_markers$SMC_arterial<-human_vasc_markers$SMC_arterial[1:100]
human_vasc_markers$SMC_venous<-human_vasc_markers$SMC_venous[1:100]





#### enrichment with bc3net
library(bc3net)
all_genes<-rownames(GerritsSmith)

enrichment_EC_markers<-enrichment(rownames(markers$EC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_FB_markers<-enrichment(rownames(markers$FB), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_PC_markers<-enrichment(rownames(markers$PC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_SMC_markers<-enrichment(rownames(markers$SMC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

write.table(rbind(enrichment_EC_markers,enrichment_FB_markers,enrichment_PC_markers,enrichment_SMC_markers),'enrichment_markers_human_vasc_with_subclusters.txt',sep='\t',row.names = F)








### dotplot enrichment markers

enrichment_markers_human_mouse_vasc <- read.delim("~/VASC/GerritsSmith/enrichment_markers_human_vasc.txt")

enrichment_markers_human_mouse_vasc$Cluster_markers<-'NA'

enrichment_markers_human_mouse_vasc$Cluster_markers[1:13]<-'EC'
enrichment_markers_human_mouse_vasc$Cluster_markers[14:26]<-'FB'
enrichment_markers_human_mouse_vasc$Cluster_markers[27:39]<-'PC'
enrichment_markers_human_mouse_vasc$Cluster_markers[40:52]<-'SMC'


enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$Cluster_markers %in% c('SMC')&enrichment_markers_human_mouse_vasc$TermID %in% c('SMC_venous','SMC_arterial'),]



####
library(ggplot2)
ggplot(enrichment_markers_human_mouse_vasc, aes(x=TermID,y=Cluster_markers))+geom_point(aes(colour=padj,size=genes))+
  scale_size(name = "size", range = c(0, 8),limits = c(20,100)) +
  scale_fill_gradient(
    low = "navy", high = "gold", name = "padj",
    guide = guide_colorbar(reverse = F),
    aesthetics = c("colour"),
    limits=c(0,1))+
  xlab('Cluster markers (from  )')+
  ylab('Cluster markers (present study)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 25),axis.text.y = element_text(size = 25))+
  theme(axis.title = element_text(size = 20))+
  #scale_x_discrete(labels=c('EC','FB','PC','SMC'))+
  #scale_x_discrete(labels=c('arterial','capillary','venous'))+
  theme(legend.text = element_text(size=15), legend.title = element_text(size=15))

ggsave('Dotplot_enrichment_human_markers_subclusters_SMC.tiff',dpi = 300,compression='lzw',width = 5,height = 6)













############# HUMAN markers Yang

Human_vasc_markers <- readxl::read_excel("~/Gene_sets/Yang_vascular_celltypes.xlsx", 
                                         sheet = 1, 
                                         skip = 0)


human_vasc_markers<-list()
human_vasc_markers[['EC']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c("BEC, Arterial","BEC, Capillary", "BEC, Venous")]))
human_vasc_markers[['FB_periv']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c('Perivascular Fibroblast')]))
human_vasc_markers[['FB_mening']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c('Meningeal Fibroblast')]))

human_vasc_markers[['PC']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c('Pericyte')]))
human_vasc_markers[['SMC']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c('SMC')]))



human_vasc_markers[['EC_arterial']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c("BEC, Arterial")]))
human_vasc_markers[['EC_capillary']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c("BEC, Capillary")]))
human_vasc_markers[['EC_venous']]<-unique(as.character(Human_vasc_markers$Gene[Human_vasc_markers$`Cell Type` %in% c("BEC, Venous")]))

saveRDS(human_vasc_markers,'human_vasc_markers_Yang.RDS')


Yang_subclusters <- read_excel("Yang_subclusters.xlsx", 
                               +     sheet = "endo")

human_vasc_markers_Yang_subclusters<-list()
human_vasc_markers_Yang_subclusters[['capillary']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Capillary']
human_vasc_markers_Yang_subclusters[['arterial']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Arterial']
human_vasc_markers_Yang_subclusters[['venous']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Venous']
human_vasc_markers_Yang_subclusters[['tip-like']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Tip-like/ Proteostatic']

Yang_subclusters <- read_excel("Yang_subclusters.xlsx",     sheet = "mural")
human_vasc_markers_Yang_subclusters[['Peric1']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Pericyte 1']
human_vasc_markers_Yang_subclusters[['Peric2']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='Pericyte 2']
human_vasc_markers_Yang_subclusters[['SMC']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='SMC']
human_vasc_markers_Yang_subclusters[['aaSMC']]<-Yang_subclusters$Gene[Yang_subclusters$`Cell subtype`=='aaSMC']


saveRDS(human_vasc_markers_Yang_subclusters,'human_vasc_markers_Yang_subclusters.RDS')



### only top 80 markers

human_vasc_markers$EC<-human_vasc_markers$EC[1:80]
human_vasc_markers$FB_periv<-human_vasc_markers$FB_periv[1:80]
human_vasc_markers$FB_mening<-human_vasc_markers$FB_mening[1:80]
human_vasc_markers$PC<-human_vasc_markers$PC[1:80]
human_vasc_markers$SMC<-human_vasc_markers$SMC[1:80]

# 
# markers$EC<-markers$EC[1:100,]
# markers$FB<-markers$FB[1:100,]
# markers$PC<-markers$PC[1:100,]
# markers$SMC<-markers$SMC[1:100,]


human_vasc_markers$capillary<-human_vasc_markers$capillary[1:50]
human_vasc_markers$arterial<-human_vasc_markers$arterial[1:50]
human_vasc_markers$venous<-human_vasc_markers$venous[1:50]
human_vasc_markers$`tip-like`<-human_vasc_markers$`tip-like`[1:50]
human_vasc_markers$Peric1<-human_vasc_markers$Peric1[1:50]
human_vasc_markers$Peric2<-human_vasc_markers$Peric2[1:50]
human_vasc_markers$SMC<-human_vasc_markers$SMC[1:50]
human_vasc_markers$aaSMC<-human_vasc_markers$aaSMC[1:50]




#### enrichment with bc3net
library(bc3net)
sce<-qs::qread('GerritsSmith_vascALL_seurat.qs')
all_genes<-rownames(sce)
enrichment_EC_markers<-enrichment(rownames(markers$EC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_FB_markers<-enrichment(rownames(markers$FB), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_PC_markers<-enrichment(rownames(markers$PC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

enrichment_SMC_markers<-enrichment(rownames(markers$SMC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')

write.table(rbind(enrichment_EC_markers,enrichment_FB_markers,enrichment_PC_markers,enrichment_SMC_markers),'enrichment_markers_human_vasc_subclusters_Yang.txt',sep='\t',row.names = F)








### dotplot enrichment markers

#enrichment_FB_markers$Cluster_markers<-'FB'

 enrichment_markers_human_mouse_vasc$Cluster_markers[1:8]<-'EC'
 enrichment_markers_human_mouse_vasc$Cluster_markers[9:16]<-'FB'
 enrichment_markers_human_mouse_vasc$Cluster_markers[17:24]<-'PC'
 enrichment_markers_human_mouse_vasc$Cluster_markers[25:32]<-'SMC'

 #enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$Cluster_markers %in% c('SMC')&enrichment_markers_human_mouse_vasc$TermID %in% c('SMC','aaSMC'),]
 
 enrichment_markers_human_mouse_vasc$TermID[enrichment_markers_human_mouse_vasc$TermID=='FB_periv']<-'FB'
 enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$TermID %in% c('FB_mening','FB_periv'),]
 enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$Cluster_markers=='FB',]
 
 
 
####
library(ggplot2)
 
 ggplot(enrichment_markers_human_mouse_vasc, aes(x=TermID,y=Cluster_markers))+geom_point(aes(colour=padj,size=genes))+
   scale_size(name = "size", range = c(0, 8),limits = c(10,80)) +
   scale_fill_gradient(
     low = "navy", high = "gold", name = "padj",
     guide = guide_colorbar(reverse = F),
     aesthetics = c("colour"),
     limits=c(0,.1))+
   xlab('Cluster markers (from  )')+
   ylab('Cluster markers (present study)')+
   theme_bw()+
   theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 25),axis.text.y = element_text(size = 25))+
   theme(axis.title = element_text(size = 20))+
   #scale_x_discrete(labels=c('EC','FB_periv','FB_mening','PC','SMC'))+
   theme(legend.text = element_text(size=15), legend.title = element_text(size=15))
 
 ggsave('Dotplot_enrichment_human_markers_Yang.tiff',dpi = 300,compression='lzw',width = 5,height = 6)
 
 ggsave('Dotplot_enrichment_human_markers_Yang_all_subclusters.tiff',dpi = 300,compression='lzw',width = 5,height = 6)
 ggsave('Dotplot_enrichment_Yang_FB_markers.tiff',dpi = 300,compression='lzw',width = 5,height = 6)
 

# enrichment_FB_markers$padj<-(-log10(enrichment_FB_markers$padj))
# 
# ggplot(enrichment_FB_markers, aes(x=TermID,y=Cluster_markers))+geom_point(aes(colour=padj,size=genes))+
#   scale_size(name = "size", range = c(0, 8), limits = c(20,80)) +
#   scale_fill_gradient(
#     high = "navy", low = "gold", name = expression("-log"[10]*"(padj)"),
#     guide = guide_colorbar(reverse = F),
#     aesthetics = c("colour"),limits=c(30,105))+
#   xlab('FB Subcluster markers (Yang et al)')+
#   ylab('Cluster markers \n(present study)')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 25),axis.text.y = element_text(size = 25))+
#   theme(axis.title = element_text(size = 20))+
#   scale_x_discrete(labels=c('Mening','Periv'))+
#   theme(legend.text = element_text(size=15), legend.title = element_text(size=15))
# 
# ggsave('Dotplot_enrichment_Yang_FB_markers.tiff',dpi = 300,compression='lzw',width = 6,height = 5)
# 


 
 
 
 
 
 
 
 
 
 ############# HUMAN markers Sun
 
 Human_vasc_markers <- readxl::read_excel("~/Gene_sets/Sun_marker_genes.xlsx", 
                                          sheet = 1, 
                                          skip = 0)
 
 
 human_vasc_markers<-list()
 human_vasc_markers[['EC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c("Endo")]))
 human_vasc_markers[['FB']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c("Fib")]))
 human_vasc_markers[['PC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c("Per")]))
 human_vasc_markers[['SMC']]<-unique(as.character(Human_vasc_markers$gene[Human_vasc_markers$cluster %in% c("SMC")]))
 
 
 
 
 saveRDS(human_vasc_markers,'human_vasc_markers_Sun.RDS')
 
 
 Sun_subclusters <- readxl::read_excel("~/Gene_sets/Sun_marker_genes.xlsx", 
                                    sheet = "subtype")
 
 human_vasc_markers_Sun_subclusters<-list()
 human_vasc_markers_Sun_subclusters[['capillary']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='capEndo']
 human_vasc_markers_Sun_subclusters[['arterial']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='aEndo']
 human_vasc_markers_Sun_subclusters[['venous']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='vEndo']
 human_vasc_markers_Sun_subclusters[['Peric1']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='Per1']
 human_vasc_markers_Sun_subclusters[['Peric2']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='Per2']
 human_vasc_markers_Sun_subclusters[['FB1']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='Fib1']
 human_vasc_markers_Sun_subclusters[['FB2']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='Fib2']
 human_vasc_markers_Sun_subclusters[['FB3']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='Fib3']
 
  human_vasc_markers_Sun_subclusters[['aSMC']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='aSMC']
 human_vasc_markers_Sun_subclusters[['vSMC']]<-Sun_subclusters$gene[Sun_subclusters$cluster=='vSMC']
 
 
 
 saveRDS(human_vasc_markers_Sun_subclusters,'human_vasc_markers_Sun_subclusters.RDS')
 
 
 
 ### only top 80 markers
 
 human_vasc_markers$EC<-human_vasc_markers$EC[1:80]
 human_vasc_markers$FB<-human_vasc_markers$FB[1:80]
 human_vasc_markers$PC<-human_vasc_markers$PC[1:80]
 human_vasc_markers$SMC<-human_vasc_markers$SMC[1:80]
 
 # 
 # markers$EC<-markers$EC[1:100,]
 # markers$FB<-markers$FB[1:100,]
 # markers$PC<-markers$PC[1:100,]
 # markers$SMC<-markers$SMC[1:100,]
 
 
 human_vasc_markers$capillary<-human_vasc_markers$capillary[1:50]
 human_vasc_markers$arterial<-human_vasc_markers$arterial[1:50]
 human_vasc_markers$venous<-human_vasc_markers$venous[1:50]
 human_vasc_markers$Peric1<-human_vasc_markers$Peric1[1:50]
 human_vasc_markers$Peric2<-human_vasc_markers$Peric2[1:50]
 human_vasc_markers$FB1<-human_vasc_markers$FB1[1:50]
 human_vasc_markers$FB2<-human_vasc_markers$FB2[1:50]
 human_vasc_markers$FB3<-human_vasc_markers$FB3[1:50]
 human_vasc_markers$aSMC<-human_vasc_markers$aSMC[1:50]
 human_vasc_markers$vSMC<-human_vasc_markers$vSMC[1:50]
 
 
 
 
 
 ### enrichment with bc3net
 library(bc3net)
 sce<-qs::qread('GerritsSmith_vascALL_seurat.qs')
 all_genes<-rownames(sce)
 enrichment_EC_markers<-enrichment(rownames(markers$EC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')
 
 enrichment_FB_markers<-enrichment(rownames(markers$FB), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')
 
 enrichment_PC_markers<-enrichment(rownames(markers$PC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')
 
 enrichment_SMC_markers<-enrichment(rownames(markers$SMC), reference =all_genes, genesets = human_vasc_markers, adj = 'bonferroni')
 
 write.table(rbind(enrichment_EC_markers,enrichment_FB_markers,enrichment_PC_markers,enrichment_SMC_markers),'enrichment_markers_human_vasc_subclusters_Sun.txt',sep='\t',row.names = F)
 
 
 
 
 
 
 
 ### dotplot enrichment markers
 
 #enrichment_FB_markers$Cluster_markers<-'FB'
 
 enrichment_markers_human_mouse_vasc$Cluster_markers[1:10]<-'EC'
 enrichment_markers_human_mouse_vasc$Cluster_markers[11:20]<-'FB'
 enrichment_markers_human_mouse_vasc$Cluster_markers[21:30]<-'PC'
 enrichment_markers_human_mouse_vasc$Cluster_markers[31:40]<-'SMC'
 
 enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$TermID %in% c('aSMC','vSMC'),]
 enrichment_markers_human_mouse_vasc<-enrichment_markers_human_mouse_vasc[enrichment_markers_human_mouse_vasc$Cluster_markers=='SMC',]
 
 
 ####
 library(ggplot2)
 
 ggplot(enrichment_markers_human_mouse_vasc, aes(x=TermID,y=Cluster_markers))+geom_point(aes(colour=padj,size=genes))+
   scale_size(name = "size", range = c(0, 8),limits = c(10,80)) +
   scale_fill_gradient(
     low = "navy", high = "gold", name = "padj",
     guide = guide_colorbar(reverse = F),
     aesthetics = c("colour"),
     limits=c(0,.1))+
   xlab('Cluster markers (from  )')+
   ylab('Cluster markers (present study)')+
   theme_bw()+
   theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 25),axis.text.y = element_text(size = 25))+
   theme(axis.title = element_text(size = 20))+
   #scale_x_discrete(labels=c('EC','FB_periv','FB_mening','PC','SMC'))+
   theme(legend.text = element_text(size=15), legend.title = element_text(size=15))
 
 ggsave('Dotplot_enrichment_human_markers_Sun.tiff',dpi = 300,compression='lzw',width = 5,height = 6)
 
 ggsave('Dotplot_enrichment_human_markers_Sun_all_subclusters_all.tiff',dpi = 300,compression='lzw',width = 5,height = 6)
 
 
 
 











################# subclustering

setwd('~/VASC/GerritsSmith/reclustered/')

reclustered<-qs::qread('reclustered_seurat_symbol.qs')






DimPlot(reclustered, reduction='UMAP_Liger',label=T,group.by='clusters')+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  theme(plot.title =element_blank())


ggsave('dimplot_reclustered_cluster_numbers_all_cluster_numbers.tiff',height=4,width=4)



######### only PCSMC cluster numbers


a<-as.data.frame(as.character(reclustered$clusters))
colnames(a)<-c('cluster')
a$cluster<-as.character(a$cluster)
a[!a$cluster %in% c('4','7','10','22','19'),'cluster']<-NA
a$cluster[a$cluster=='4']<-'PC1'
a$cluster[a$cluster=='7']<-'PC2'
a$cluster[a$cluster=='10']<-'PC3'
a$cluster[a$cluster=='22']<-'PC4'
a$cluster[a$cluster=='19']<-'SMC'



reclustered[['clusters']]<-as.character(a$cluster)

d<-DimPlot(reclustered, reduction='UMAP_Liger',label=F,group.by='clusters')+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  theme(plot.title =element_blank())


ggsave(d,'dimplot_reclustered_cluster_numbers_only_PCSMC.tiff',height=4,width=4)

## for source

write.table(as.data.frame(d[[1]][['data']]),'reclustered_dimplot_with_cluster_numbers_only_PCSMC_for_source.txt',sep='\t', row.names = T,col.names = NA)




####### with cluster names
setwd('~/VASC/GerritsSmith/reclustered/')

reclustered<-qs::qread('reclustered_seurat_symbol.qs')

a<-as.data.frame(cbind(as.character(reclustered$clusters),reclustered$cluster_celltype))
colnames(a)<-c('cluster','cell')
a$cell<-as.character(a$cell)
a$cluster<-as.character(a$cluster)
a[a$cluster %in% c('19'),'cell']<-'SMC'
a[a$cluster %in% c('4','7','10','22'),'cell']<-'PC'
a[a$cluster %in% c('13','8','1','2','9','5','20'),'cell']<-'EC'
a[a$cluster %in% c('18','21','23','16','15','3','12','6','14','17','11'),'cell']<-'FB'

#
#a[!a$cluster %in% c('4','7','10','22','19'),'cell']<-NA
#

reclustered[['cluster_celltype']]<-as.character(a$cell)

DimPlot(reclustered, reduction='UMAP_Liger',label=F,group.by='cluster_celltype')+NoLegend()+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  theme(plot.title =element_blank())

ggsave('dimplot_reclustered_only_PCSMC_REVISION_sameFont_as_GlobalUMAP_withAllColours.tiff',height=4,width=4)




source_data<-as.data.frame(cbind(as.character(reclustered@meta.data$barcode),as.character(reclustered$cluster_celltype),reclustered@reductions$UMAP_Liger@cell.embeddings))
colnames(source_data)[1:2]<-c('barcode','cell type')


write.table(source_data,file = 'source_1c.txt',sep='\t',row.names = F)


##############

reclustered<-SetIdent(reclustered,value='clusters')
markers19<-FindMarkers(reclustered, ident.1='19',ident.2=c('22','4','7','10'),only.pos=T)

markers7<-FindMarkers(reclustered, ident.1='7',ident.2=c('22','4','19','10'),only.pos=T)

markers4<-FindMarkers(reclustered, ident.1='4',ident.2=c('22','7','19','10'),only.pos=T)


markers22<-FindMarkers(reclustered, ident.1='22',ident.2=c('4','7','19','10'),only.pos=T)

markers10<-FindMarkers(reclustered, ident.1='10',ident.2=c('4','7','19','22'),only.pos=T)

mark<-list()

mark[['cl19']]<-markers19
mark[['cl7']]<-markers7
mark[['cl4']]<-markers4
mark[['cl22']]<-markers22
mark[['cl10']]<-markers10


saveRDS(mark,'markers_PCSMC.RDS')


reclustered<-subset(reclustered,subset=clusters%in%c('PC1','PC2','PC3','PC4','SMC'))
VlnPlot(reclustered,features=c('ACTA2','RGS5'),group.by='clusters',pt.size=0)


ggsave('VlnPlot_pcsmc.tiff',height=4,width=6)

VlnPlot(reclustered,features=c('MYH11','GRM8'),group.by='clusters',pt.size=0)

ggsave('VlnPlot_pcsmc2.tiff',height=4,width=6)


##### marker comparisons for source
reclustered<-SetIdent(reclustered,value='clusters')
markers19<-FindMarkers(reclustered, ident.1='19',ident.2=c('22','4','7','10'),only.pos=T, test.use = 'MAST')

markers7<-FindMarkers(reclustered, ident.1='7',ident.2=c('19'),only.pos=T, test.use = 'MAST')


markers4<-FindMarkers(reclustered, ident.1='4',ident.2=c('19'),only.pos=T, test.use = 'MAST')



markers22<-FindMarkers(reclustered, ident.1='22',ident.2=c('19'),only.pos=T, test.use = 'MAST')


markers10<-FindMarkers(reclustered, ident.1='10',ident.2=c('19'),only.pos=T, test.use = 'MAST')


mark<-list()

mark[['cl19']]<-markers19
mark[['cl7']]<-markers7
mark[['cl4']]<-markers4
mark[['cl22']]<-markers22
mark[['cl10']]<-markers10


saveRDS(mark,'markers_PCSMC_corresponding_to_vlnplot_for_source.RDS')




### gene expression for source
source<-VlnPlot(reclustered,features=c('ACTA2','RGS5','MYH11','GRM8'),group.by='clusters',pt.size=0)

sc<-cbind(source[[1]][['data']],source[[2]][['data']],source[[3]][['data']],source[[4]][['data']])

sc<-sc[,-c(2,4,6)]
write.table(sc,'vlnPCSMC_for_source.txt',sep='\t',row.names = T,col.names = NA)
########




a<-as.data.frame(cbind(as.character(reclustered$clusters),reclustered$cluster_celltype))
colnames(a)<-c('cluster','cell')
a$cell<-as.character(a$cell)
a$cluster<-as.character(a$cluster)
a[a$cluster %in% c('19'),'cell']<-'SMC'
a[a$cluster %in% c('4','7','10','22'),'cell']<-'PC'

reclustered[['cluster_celltype']]<-as.character(a$cell)

qs::qsave(reclustered, 'reclustered_seurat_symbol.qs')



########### violin plots with astro/micro cell markers


GerritsSmith <- readRDS("./vasc_ALL_new_celltype_for_heatmap.RDS")


VlnPlot(GerritsSmith,features=c("HLA-DRA", "CX3CR1", "C1QB", "CSF1R", "CD74", "LPAR6", "C3", "DOCK8", "SYK", "P2RY12"),group.by='cluster_celltype',pt.size=0.3)
VlnPlot(GerritsSmith,features=c("AQP4", "SLC1A2", "GFAP", "FGFR3", "SLC14A1", "TNC", "CD44"),group.by='cluster_celltype',pt.size=0.3)





########## cluster markers

markers_EC_FB <- readRDS("~/VASC/GerritsSmith/markers_EC_FB.RDS")

reclustered<-SetIdent(reclustered,value='cluster_celltype')

markers_EC_FB[['PC']]<-FindMarkers(reclustered, ident.1 = 'PC',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')

markers_EC_FB[['SMC']]<-FindMarkers(reclustered, ident.1 = 'SMC',only.pos = T,slot = 'data',assay = 'RNA',test.use = 'MAST')

saveRDS(markers_EC_FB,'markers_EC_FB.RDS')








####### dotplot GWAS
library(Seurat)
vasc <- readRDS("~/VASC/GerritsSmith/vasc_ALL_new_celltype_for_heatmap.RDS")
#vasc<-subset(vasc, subset=cluster_celltype!='SMC')


Important_gene_sets<-readRDS('~/AD/Gene_sets/Important_GeneSets.RDS')
GWAS<-Important_gene_sets[['GWAS_extended']]
GWAS<-GWAS[order(GWAS,decreasing = T)]

DotPlot(vasc,features = GWAS,group.by = 'cluster_celltype',scale.by = 'size',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .05)+
  coord_flip()+theme(axis.text.x = element_text(angle = 90))

DotPlot(vasc,features = GWAS,group.by = 'cluster_celltype',scale.by = 'size',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .05)+
  coord_flip()+theme(axis.text.x = element_text(angle = 90))

ggsave('Dotplot_GWAS_vasc_ECFBPCSMC.tiff',width = 4,height = 10)

source_data<-DotPlot(vasc,features = GWAS,group.by = 'cluster_celltype',scale.by = 'size',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .05)+
  coord_flip()+theme(axis.text.x = element_text(angle = 90))
source_data<-source_data$data
write.table(source_data,'source_data_1H.txt',sep = '\t',col.names = NA)

##### AUCell 
library(AUCell)
human_vasc_markers_Yang <- readRDS("~/VASC/GerritsSmith/human_vasc_markers_Yang.RDS")

cells_rankings <- readRDS("~/VASC/GerritsSmith/MEGENA_fibroblast_ALL/cells_rankings_fibro.RDS")

cells_AUC<-AUCell_calcAUC(human_vasc_markers_Yang[2:3], cells_rankings)

saveRDS(cells_AUC,'FB_markers_Yang_AUC.RDS')



############### UMAP plots by sex etc
library(ggplot2)
object<-qs::qread('GerritsSmith_seurat.qs')
DimPlot(object, reduction='UMAP_Liger',split.by = 'sex',group.by = 'sex',raster=F)+NoLegend()

ggplot2::ggsave('dimplot_GerritsSmith_bySex2.tiff',height=4,width=4)


a<-as.data.frame(as.character(object$diagnosis))
a$`as.character(object$diagnosis)`<-as.character(a$`as.character(object$diagnosis)`)
a$`as.character(object$diagnosis)`[a$`as.character(object$diagnosis)`=='CTR']<-'NDC'
object$diagnosis<-as.character(a$`as.character(object$diagnosis)`)
###
DimPlot(object, reduction='UMAP_Liger',split.by = 'diagnosis',group.by = 'diagnosis',raster=F)+NoLegend()
ggplot2::ggsave('dimplot_GerritsSmith_byDiagnosis2.tiff',height=4,width=4)

###

a<-as.data.frame(as.character(object$brain_region))
a$`as.character(object$brain_region)`<-as.character(a$`as.character(object$brain_region)`)
a$`as.character(object$brain_region)`[a$`as.character(object$brain_region)`=='EC']<-'EntC'
object$brain_region<-as.character(a$`as.character(object$brain_region)`)

DimPlot(object, reduction='UMAP_Liger',split.by = 'brain_region',group.by = 'brain_region',raster = F)+NoLegend()+ggtitle('brain region')
ggplot2::ggsave('dimplot_GerritsSmith_byRegion2.tiff',height=4,width=8)

###
DimPlot(object, reduction='UMAP_Liger',split.by = 'dataset',group.by = 'dataset')+NoLegend()
ggplot2::ggsave('dimplot_GerritsSmith_byDataset.tiff',height=4,width=4)





############### UMAP plots by sex etc in the unenriched dataset
library(ggplot2)
sce<-qs::qread('~/AD/Basic_data_Glia/unenriched/unenriched_AD/final_seurat_ECSSC.qs')


DimPlot(sce, reduction='UMAP_Liger',split.by = 'sex',group.by = 'sex',raster = F)+NoLegend()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title = element_text(size = 20))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())

ggplot2::ggsave('dimplot_unenriched_bySex2.tiff',height=4,width=8)


a<-as.data.frame(as.character(sce$diagnosis))
a$`as.character(sce$diagnosis)`<-as.character(a$`as.character(sce$diagnosis)`)
a$`as.character(sce$diagnosis)`[a$`as.character(sce$diagnosis)`=='Control']<-'NDC'
sce$diagnosis<-as.character(a$`as.character(sce$diagnosis)`)
###
DimPlot(sce, reduction='UMAP_Liger',split.by = 'diagnosis',group.by = 'diagnosis',raster = F)+NoLegend()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title = element_text(size=20))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())

ggplot2::ggsave('dimplot_unenriched_byDiagnosis2.tiff',height=4,width=8)
###

a<-as.data.frame(as.character(sce$brain_region))
a$`as.character(sce$brain_region)`<-as.character(a$`as.character(sce$brain_region)`)
a$`as.character(sce$brain_region)`[a$`as.character(sce$brain_region)`=='EC']<-'EntC'
sce$brain_region<-as.character(a$`as.character(sce$brain_region)`)
DimPlot(sce, reduction='UMAP_Liger',split.by = 'brain_region',group.by = 'brain_region',raster = F)+NoLegend()+ggtitle('brain region')+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(axis.title = element_text(size = 20))+
  xlab('umap_liger_1')+
  ylab('umap_liger_2')+
  theme(plot.title =element_blank())

ggplot2::ggsave('dimplot_unenriched_byRegion2.tiff',height=4,width=8)







###### vlnplot production of important genes
library(ggplot2)
object =readRDS('~/VASC/GerritsSmith//CellChat_wPVM/vasc_PVM_for_CellChat.RDS')
object =qs::qread('~/VASC/GerritsSmith/GerritsSmith_seurat.qs')


VlnPlot(object,features=c('VEGFA','ANGPT2','FGF2'),group.by='cluster_celltype',pt.size=0)


ggplot2::ggsave('~/VASC/GerritsSmith/VlnPlot_VEFGA_ANGPT2_FGF2.tiff',height=4,width=6)



DotPlot(object,features = c('VEGFA','ANGPT2','FGF2'),group.by = 'cluster_celltype',scale.by = 'size',cols = RColorBrewer::brewer.pal(8,'RdBu')[c(8,1)], dot.min = .05)+
  coord_flip()+theme(axis.text.x = element_text(angle = 90))

ggplot2::ggsave('~/VASC/GerritsSmith/DotPlot_VEFGA_ANGPT2_FGF2_GerritsSmithALL.tiff',height=4,width=6)




########## EWCE 


######### EWCE 

library(EWCE)
vasc_ALL<-qs::qread('./GerritsSmith_vascALL_seurat.qs')

ctd <- readRDS('~/Bipolar/Proteomics/ctd_fromEWCE.RDS')








hits<-readRDS("~/VASC/GerritsSmith/markers_EC_FB.RDS")



gene_set<-rownames(hits$FB[1:400,])
background<-rownames(vasc_ALL)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "mouse",
                                                genelistSpecies = "human",
                                                hits =gene_set, 
                                                reps = 10000,
                                                annotLevel = 1,
                                                bg=background)
# ,
#                                                 method = 'homologene')

knitr::kable(full_results$results)



write.table(full_results$results,'EWCE_results_table_FB_markers.txt',sep='\t',row.names = F)





plot_list <- EWCE::ewce_plot(total_res = full_results$results,
                             mtc_method = "fdr",
                             ctd = ctd)

p<-plot_list$plain
p$data$CellType<-as.character(rownames(p$data))

ggsave('EWCE_FB_markers.tiff',plot = p,width = 7,height = 5)


print(plot_list$plain)




########### stackplot cell numbers

#GerritsSmith<-qs::qread('GerritsSmith_vascALL_seurat.qs')

setwd("~/VASC/GerritsSmith")
library(graphics)
library(SingleCellExperiment)
library(Seurat)

#### 
GerritsSmith<-qs::qread('./GerritsSmith_ensg_sce.qs')
sce.rowdata <- read.delim("~/VASC/GerritsSmith/final/SCE/final_sce/sce-rowdata.tsv")
rownames(GerritsSmith)<-as.character(sce.rowdata$gene)

reclustered<-qs::qread('./reclustered/reclustered_seurat_symbol.qs')
pericytes<-colnames(subset(reclustered,subset=clusters %in% c('4','7','10','22')))
smc<-colnames(subset(reclustered,subset=clusters %in% c('19')))

a<-as.data.frame(as.character(GerritsSmith$clusters))
rownames(a)<-colnames(GerritsSmith)
colnames(a)<-'clusters'
a$clusters_for_markers<-a$clusters
a$clusters_for_markers[a$clusters %in% c('1','2','9')]<-'Micro'
a$clusters_for_markers[a$clusters %in% c('3','4','5')]<-'Astro'
a$clusters_for_markers[a$clusters %in% c('6')]<-'Oligo'
a$clusters_for_markers[a$clusters %in% c('8')]<-'EC'
a$clusters_for_markers[a$clusters %in% c('10')]<-'PC/SMC'
a$clusters_for_markers[a$clusters %in% c('7')]<-'FB'
a[intersect(pericytes,rownames(a)),'clusters_for_markers']<-'PC'
a[intersect(smc,rownames(a)),'clusters_for_markers']<-'SMC'


GerritsSmith$cluster_celltype<-as.character(a$clusters_for_markers)
sce<-GerritsSmith[,colData(GerritsSmith)$cluster_celltype%in%c('PC','SMC','EC','FB','Astro','Micro')]






mat<-as.data.frame(cbind(as.character(sce$brain_region),as.character(sce$cluster_celltype),as.character(sce$sex),as.character(sce$diagnosis)))
colnames(mat)<-c('brain_region','cell_type','sex','diagnosis')
mat$brain_region[mat$brain_region=='EC']<-'EntC'
mat$cell_type<-factor(mat$cell_type,levels = c('EC','FB','PC','SMC','Astro','Micro'))
mat<-mat[!mat$cell_type%in%c('Astro','Micro'),]

library(ggplot2)

ggplot(data=mat, aes(x=brain_region,fill=cell_type))+geom_bar(stat = 'count', position = 'fill')+
  theme_classic()+theme(axis.text.y=element_text(size=25),axis.text.x = element_text(angle = 45,hjust = 1,size = 30), axis.title =element_text(size=30))+
 # scale_y_discrete(labels=c('EC','FB','PC/SMC'))+
  xlab('Brain region')+ylab('% Cell count')+
  theme(legend.text = element_text(size = 25),legend.title = element_text(size = 25))+
  labs(fill='Cell type')
ggsave(filename = 'stackplot_cell_count_byRegion_larger_font.tiff',width = 6,height = 6)


ggplot(data=mat, aes(x=sex,fill=cell_type))+geom_bar(stat = 'count', position = 'fill')+
  theme_classic()+theme(axis.text=element_text(size=20), axis.title =element_text(size=25))+
  # scale_y_discrete(labels=c('EC','FB','PC/SMC'))+
  xlab('Cluster')+ylab('% Cell count')+
  theme(legend.text = element_text(size = 15),legend.title = element_text(size = 20))
ggsave(filename = 'stackplot_cell_count_bySex.tiff',width = 6,height = 6)


ggplot(data=mat, aes(x=diagnosis,fill=cell_type))+geom_bar(stat = 'count', position = 'fill')+
  theme_classic()+theme(axis.text=element_text(size=20), axis.title =element_text(size=25))+
  # scale_y_discrete(labels=c('EC','FB','PC/SMC'))+
  xlab('Cluster')+ylab('% Cell count')+
  theme(legend.text = element_text(size = 15),legend.title = element_text(size = 20))
ggsave(filename = 'stackplot_cell_count_byDiagnosis.tiff',width = 6,height = 6)


### for source

t<-table(mat$brain_region,mat$cell_type)
t<-t/rowSums(t)*100

write.table(t,'stackplot_for_source.txt',sep='\t',row.names = T)
