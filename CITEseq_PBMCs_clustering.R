#CITE-seq

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(Polychrome)
library(patchwork)
library(BiocManager)
library(sctransform)
library(tidyverse)
library(hdf5r)
library(glmGamPoi)
library(SeuratDisk)
library(scran)
library(future) 
library(data.table)
library(tidyr)
library(polychrome)

sample_info <- data.frame( 
  'disease' = c("SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","SSc","HD","HD","HD","HD","HD","HD"),
  'path' = paste0(fileroot, c('SSc1_before.h5',"SSc2.h5","SSc3.h5","SSc4.h5","SSc5.h5",'SSc6.h5',"SSc7.h5","SSc8.h5","SSc9.h5","SSc10.h5","SSc11.h5","SSc12.h5","SSc13.h5","SSc14.h5","SSc15.h5",'SSc16.h5',"SSc17.h5","SSc18.h5","SSc19.h5","SSc20.h5","SSc21.h5","HD1.h5","HD2.h5","HD3.h5","HD4.h5","HD5.h5","HD6.h5")),
  'sample' = c('SSc1_before',"SSc2","SSc3","SSc4","SSc5",'SSc6',"SSc7","SSc8","SSc9","SSc10","SSc11","SSc12","SSc13","SSc14","SSc15",'SSc16',"SSc17","SSc18","SSc19","SSc20","SSc21","HD1","HD2","HD3","HD4","HD5","HD6"),
  'SRC' = c('yes',"no", "no",'yes',"no", 'no',"no", "no",'no','no',"no", 'no',"no","no", 'no', "no", 'no',"no","yes","yes","no","HD","HD","HD", "HD","HD", "HD"),
  'ILD' = c('yes',"no", 'yes',"yes", 'no',"yes", "yes",'yes','no',"yes", 'yes',"no", 'yes','yes','no',"yes", 'yes',"no", 'no', "no","no","HD","HD","HD","HD","HD","HD"),
  'ILDwithoutSRC' = c('no',"no", 'yes',"no", 'no',"yes", "yes",'yes','no',"yes", 'yes',"no", 'yes','yes','no',"yes", 'yes',"no", 'no', "no","no","HD","HD","HD","HD","HD","HD"))
sample_info

load_and_prep_RNA <- function(x){
  tmp.data <- Read10X_h5( x['path'])
  tmprna <- CreateSeuratObject(counts = tmp.data[["Gene Expression"]], assay = "RNA")
  rownames(x=tmp.data[["Antibody Capture"]]) <- paste0("ADT-", rownames(tmp.data[["Antibody Capture"]]))
  tmp.data <- tmp.data[["Antibody Capture"]]
  tmpadt<- CreateSeuratObject(counts = tmp.data, assay = "ADT")
  tmprna[["ADT"]]<-tmpadt[["ADT"]]
  tmprna <- PercentageFeatureSet(tmprna, pattern="^MT-", col.name="percent.mt")
  tmprna <- subset(tmprna, subset = 
                     percent.mt < 20 & 
                     nFeature_RNA < 5000 & 
                     nFeature_RNA > 200)
  tmprna <- AddMetaData(tmprna, metadata = x['disease'], col.name="disease")
  tmprna <- AddMetaData(tmprna, metadata = x['sample'], col.name="sample")
  tmprna <- AddMetaData(tmprna, metadata = x['SRC'], col.name="SRC")
  tmprna <- AddMetaData(tmprna, metadata = x['ILD'], col.name="ILD")
  tmprna <- SCTransform(tmprna, method="glmGamPoi", vars.to.regress = c("percent.mt"), seed.use=1448145)
  return(tmprna)}
samples_rna <- apply(sample_info, 1, load_and_prep_RNA)
gc()

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DefaultAssay(s.int) <- 'SCT'
anchors <- FindTransferAnchors(
  reference = reference,
  query = s.int,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
gc()

s.int <- MapQuery(
  anchorset = anchors,
  query = s.int,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

Idents(s.int)<-"predicted.celltype.l2"
DimPlot(s.int, split.by = NULL, reduction = 'ref.umap', label = TRUE, repel = TRUE, label.box=TRUE, label.size = 2) 

Mo<-subset(s.int, ident=c("CD14 Mono","CD16 Mono","cDC2","cDC1"))
CD4<-subset(s.int, ident=c("CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","Treg"))
CD8<-subset(s.int, ident=c("CD8 Naive","CD8 TCM","CD8 TEM"))
B<-subset(s.int, ident=c("B naive","B intermediate","B memory"))
NK<-subset(s.int, ident=c("NK","NK_CD56bright"))
pDC<-subset(s.int, ident=c("pDC"))

s.int <- RunUMAP(s.int, reduction = "ref.spca", dims = 1:30, seed.use = 42)
s.int <- FindNeighbors(s.int, reduction = "ref.spca", dims = 1:30)
s.int <- FindClusters(s.int, graph.name="SCT_snn",algorithm = 3,resolution =0.5,verbose = FALSE)
gc()
#Manually remove doublet, platelets, and red blood cells using gene expression data and surface antigen data (Omitted).
DimPlot(s.int, split.by = "disease", reduction = 'umap',label = TRUE, repel = FALSE, label.box=TRUE, label.size = 2) 
DimPlot(s.int, split.by = "sample", reduction = 'umap',label = TRUE, repel = FALSE, label.box=TRUE, label.size = 2) 
saveRDS(s.int, "s.int.rds")

Mo <- RunUMAP(Mo, reduction = "ref.spca", dims = 1:30, seed.use = 42)
Mo <- FindNeighbors(Mo, reduction = "ref.spca", dims = 1:30)
Mo <- FindClusters(Mo, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(Mo, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(Mo, cols= c("red","yellow","white"), col.max=2, col.min=-2, scale=TRUE, dot.scale=5,
        features = c("EGR1","NFKBIA","IER2","SGK1","CCR1","IL1B","MX1","HERC5","ISG15","IFI44L","IFI44","S100A8","PLBD1","S100A12","S100A9","PADI4","VCAN","CD14","LYZ","CD99","TAGLN2","HLA-DRB1","HLA-DMA","HLA-DMB","HLA-DRB5","HLA-DRA","APOBEC3A","PLAC8","PSME2","CD52","IFI30","CDKN1C","FCGR3A","LYPD2","RHOC","MS4A7","IFITM1","IFITM3","IFITM2","TCF7L2","HES4","CLEC9A","IRF8","CPVL","SNX3","WDFY4","CD1C","FCER1A","CLEC10A","ENHO","PLD4"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=12, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(Mo, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
Idents(Mo)<-"disease"
Mo_2 <- subset(Mo, ident= "SSc")
DimPlot(Mo_2, split.by = "SRC",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
saveRDS(Mo, "Mo.rds")
saveRDS(Mo_2, "Mo_2.rds")


CD4 <- RunUMAP(CD4, reduction = "ref.spca", dims = 1:30, seed.use = 42)
CD4 <- FindNeighbors(CD4, reduction = "ref.spca", dims = 1:30)
CD4 <- FindClusters(CD4, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(CD4, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(CD4, cols= c("red","yellow","white"), col.max=1.8, col.min=-1, scale=TRUE, dot.scale=5,
        features = c("CCR7","LEF1","SATB1","CD7","RPS13","IL7R","INPP4B","AP3M2","FXYD5","LTB","GATA3","S100A4","CRIP1","S100A10","ANXA1","FOXP3","CTLA4","TIGIT","IL10RA","IL32","GZMK","DUSP2","LYAR","TNFAIP3","CD74","NKG7","GZMH","PRF1","GZMB","GZMA","MX1","IFI6","IFI44L","ISG15","STAT1"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=9, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(CD4, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
Idents(CD4)<-"disease"
CD4_2 <- subset(CD4, ident= "SSc")
DimPlot(CD4_2, split.by = "ILD",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
saveRDS(CD4, "CD4.rds")
saveRDS(CD4_2, "CD4_2.rds")


CD8 <- RunUMAP(CD8, reduction = "ref.spca", dims = 1:30, seed.use = 42)
CD8 <- FindNeighbors(CD8, reduction = "ref.spca", dims = 1:30)
CD8 <- FindClusters(CD8, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(CD8, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(CD8_3, cols= c("red","yellow","white"), col.max=1.7, col.min=-1.5, scale=TRUE, dot.scale=8,
        features = c("CCR7","LEF1","TCF7","RPS13","SELL","IL7R","LTB","JUNB","TXNIP","GATA3","FXYD5","PLP2","FXYD7","CRIP2","GZMK","CMC1","DUSP2","EOMES","CCL5","GBP5","CXCR3","COTL1","HLA-DRB1","CD74","GZMB","GZMH","PRF1","GNLY","NKG7"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=12, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(CD4, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
Idents(CD8)<-"disease"
CD8_2 <- subset(CD8, ident= "SSc")
DimPlot(CD8_2, split.by = "ILD",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
saveRDS(CD8, "CD8.rds")
saveRDS(CD8_2, "CD8_2.rds")

B <- RunUMAP(B, reduction = "ref.spca", dims = 1:30, seed.use = 42)
B <- FindNeighbors(B, reduction = "ref.spca", dims = 1:30)
B <- FindClusters(B, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(B, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(B, cols= c("red","yellow","white"), col.max=2, col.min=-2, scale=TRUE, dot.scale=8,
        features = c("TCL1A","IGHD","BACH2","YBX3","IL4R","IGHM","GPR183","FCRL2","TNFRSF13B","MARCKS","CD27","COCH","ITGB1","AIM2","CRIP1","ITGAX","TBX21","TNFRSF1B","FGR","IGHG1"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=12, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(B, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
saveRDS(B, "B.rds")

NK <- RunUMAP(NK, reduction = "ref.spca", dims = 1:30, seed.use = 42)
NK <- FindNeighbors(NK, reduction = "ref.spca", dims = 1:30)
NK <- FindClusters(NK, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(NK, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(NK, cols= c("red","yellow","white"), col.max=2, col.min=-2, scale=TRUE, dot.scale=8,
        features = c("SELL","GZMK","TCF7","IL7R","NCAM1","ITGA6","PIK3R1","KLRG1","CEBPD","CXCR4","GZMH","GZMB","PRF1","FGFBP2","FCGR3A"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=12, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(NK, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
saveRDS(NK, "NK.rds")

pDC <- RunUMAP(pDC, reduction = "ref.spca", dims = 1:30, seed.use = 42)
pDC <- FindNeighbors(pDC, reduction = "ref.spca", dims = 1:30)
pDC <- FindClusters(pDC, graph.name="SCT_snn",algorithm = 3,resolution = 0.5,verbose = FALSE)
#Manually remove doublet and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(pDC, split.by = "disease",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2) 
DotPlot(pDC, cols= c("red","yellow","white"), col.max=1, col.min=-0.2, scale=TRUE, dot.scale=8,
        features = c("MX1","IFI44L","IRF8","IRF7","TLR7","TLR9","HLA-DRB1","TCF4"),
        group.by="seurat_clusters") + theme(axis.text.x=element_text(size=12, angle=90)) +
  scale_color_gradient2(low="blue", mid="white", high="red")
DimPlot(pDC, split.by = "sample",reduction = 'umap', label = TRUE, repel = FALSE, label.size = 2)
saveRDS(pDC, "pDC.rds")
