#Abseq

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(data.table)
library(tidyr)
library(glmGamPoi)

sample_info <- data.frame( 
  'disease' = c("SSc"),
  'path' = paste0(fileroot, c("AI122-1.txt")),
  'sample' = c("AI122-1"))
sample_info

sample_info2 <- data.frame( 
  'disease' = c("kidney"),
  'path' = paste0(fileroot, c("AI122-1_kidney.txt")),
  'sample' = c("AI122-1_kidney"))
sample_info2

load_and_prep_RNA <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  counts <- select(counts, -c(18,28,30))
  x.rn <- counts[,36:length(counts),]
  x.rna <- CreateSeuratObject(counts = t(x.rn))
  x.rna <- AddMetaData(x.rna, metadata = x['disease'], col.name="disease")
  x.rna <- AddMetaData(x.rna, metadata = x['sample'], col.name="sample")
  x.rna[["percent.mt"]] <- PercentageFeatureSet(x.rna, pattern = "^MT.")
  x.rna <- subset(x.rna, subset = 
                    percent.mt < 20 & 
                    nFeature_RNA < 5000 & 
                    nFeature_RNA > 200)
  x.rna<-SCTransform(x.rna, method="glmGamPoi",return.only.var.genes = FALSE,variable.features.n = 23000, seed.use=1448145)
  return(x.rna)}
samples_rna <- apply(sample_info, 1, load_and_prep_RNA)

gc()

load_and_prep_ADT <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  counts <- select(counts, -c(18,28,30))
  x.rn <- counts[,36:length(counts),]
  x.rna <- CreateSeuratObject(counts = t(x.rn), assay ="RNA")
  x.rna[["percent.mt"]] <- PercentageFeatureSet(x.rna, pattern = "^MT.")
  x.rna <- subset(x.rna, subset = 
                    percent.mt < 20 & 
                    nFeature_RNA < 5000 & 
                    nFeature_RNA > 200)
  x.ab <- counts[,1:35,]
  x.ab <- x.ab[colnames(x.rna),]
  x.adt <- CreateSeuratObject(counts = t(x.ab), assay ="ADT")
  x.adt <- AddMetaData(x.adt, metadata = x['disease'], col.name="disease")
  x.adt <- AddMetaData(x.adt, metadata = x['sample'], col.name="sample")
  x.adt <- NormalizeData(x.adt, normalization.method="CLR", margin = 2) ## normalization across cells 
  return(x.adt)}
samples_adt <- apply(sample_info, 1, load_and_prep_ADT)

gc()

load_and_prep_RNA2 <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  counts <- select(counts, -c(18,28,30))
  x.rn <- counts[,36:length(counts),]
  x.rna <- CreateSeuratObject(counts = t(x.rn))
  x.rna <- AddMetaData(x.rna, metadata = x['disease'], col.name="disease")
  x.rna <- AddMetaData(x.rna, metadata = x['sample'], col.name="sample")
  x.rna[["percent.mt"]] <- PercentageFeatureSet(x.rna, pattern = "^MT.")
  x.rna <- subset(x.rna, subset = 
                    percent.mt < 40 & 
                    nFeature_RNA < 5000 & 
                    nFeature_RNA > 200)
  x.rna<-SCTransform(x.rna, method="glmGamPoi",return.only.var.genes = FALSE,variable.features.n = 23000, seed.use=1448145)
  return(x.rna)}
samples_rna2 <- apply(sample_info2, 1, load_and_prep_RNA2)

gc()

load_and_prep_ADT2 <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  counts <- select(counts, -c(18,28,30))
  x.rn <- counts[,36:length(counts),]
  x.rna <- CreateSeuratObject(counts = t(x.rn), assay ="RNA")
  x.rna[["percent.mt"]] <- PercentageFeatureSet(x.rna, pattern = "^MT.")
  x.rna <- subset(x.rna, subset = 
                    percent.mt < 40 & 
                    nFeature_RNA < 5000 & 
                    nFeature_RNA > 200)
  x.ab <- counts[,1:35,]
  x.ab <- x.ab[colnames(x.rna),]
  x.adt <- CreateSeuratObject(counts = t(x.ab), assay ="ADT")
  x.adt <- AddMetaData(x.adt, metadata = x['disease'], col.name="disease")
  x.adt <- AddMetaData(x.adt, metadata = x['sample'], col.name="sample")
  x.adt <- NormalizeData(x.adt, normalization.method="CLR", margin = 2) ## normalization across cells 
  return(x.adt)}
samples_adt2 <- apply(sample_info2, 1, load_and_prep_ADT2)

gc()

samples_rna[[2]] <- samples_rna2[[1]]
samples_adt[[2]] <- samples_adt2[[1]]

features_ADT <- rownames(samples_adt[[1]])
samples_adt <- lapply(X=samples_adt, FUN=function(x){
  x <- ScaleData(x, features=features_ADT, verbose=FALSE)
  x <- RunPCA(x, features=features_ADT, verbose=FALSE)
})
gc()

adt.anchors <- FindIntegrationAnchors(object.list = samples_adt,reduction="rpca",k.anchor = 15)
s.int <- IntegrateData(anchorset = adt.anchors)
gc()

features_RNA <- SelectIntegrationFeatures(object.list = samples_rna, assay=rep('SCT', length(samples_rna)), nfeatures = 23000)
samples_rna <- PrepSCTIntegration(object.list = samples_rna, anchor.features = features_RNA)
samples_rna <- lapply(X=samples_rna, FUN = RunPCA, features=features_RNA)
rna.anchors <- FindIntegrationAnchors(object.list = samples_rna, normalization.method='SCT', anchor.features = features_RNA, reduction="rpca", k.anchor=15)
gc()
int.rna <- IntegrateData(anchorset = rna.anchors, normalization.method="SCT")
gc()

s.int[["RNA"]] <- int.rna[["RNA"]]
s.int[["nFeature_RNA"]] <- int.rna[["nFeature_RNA"]]
s.int[["nCount_RNA"]] <- int.rna[["nCount_RNA"]]
s.int[["percent.mt"]] <- int.rna[["percent.mt"]]
s.int[["SCT"]] <- int.rna[["integrated"]]
s.int[["nFeature_SCT"]] <- int.rna[["nFeature_SCT"]]
s.int[["nCount_SCT"]] <- int.rna[["nCount_SCT"]]

DefaultAssay(s.int) <- 'SCT'
s.int <- RunPCA(s.int, verbose=TRUE)
s.int <- FindNeighbors(s.int, reductiobn="pca")
s.int <- RunUMAP(s.int, reductiobn="pca", dims=1:30, , seed.use = 42)
s.int <- FindClusters(s.int, reductiobn="pca", graph.name="SCT_snn",resolution=0.5)
Mo <- subset(s.int, idents =c("4","5","14")) #subset monocytes, macrophages, and conventional dendritic cells from "s.int".
#Manually remove doublets, platelets, and red blood cells from "s.int" using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(s.int,split.by="sample", reduction = 'umap', label = F, repel = TRUE, label.box=T, label.size = 3) 
#Manually remove doublets and debris using gene expression data and surface antigen data. Re-cluster and annotate cell populations (Omitted).
DimPlot(Mo,split.by="sample", reduction = 'umap', label = F, repel = TRUE, label.box=T, label.size = 3)
