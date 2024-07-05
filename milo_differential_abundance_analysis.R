#Milo differential abundance analysis

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(data.table)
library(tidyr)
library(miloR)
library(patchwork)
library(scater)
library(scran)
library(ggraph)

#Differential abundance analysis of PBMCs, between SSc patients with scleroderma renal crisis (SRC) and those without SRC.
s.int[["seurat_clusters"]] <- Idents(s.int)
s.int$SRC <- as.character(s.int$SRC)
Idents(s.int) <- "disease"
s.int_SSc <- subset(s.int, ident = "SSc")
Idents(s.int_SSc) <- "SRC"
subsampled_cells <- s.int_SSc@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  sample_n(1000) %>%
  pull(cell_id)
s.int_SSc_1000 <- subset(s.int_SSc, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(s.int_SSc_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.5, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,56,82,72,64,60,80,58,47,76,81,57,46,61,50,62,77,67,53,74,73)) #age of all SSc patients
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+SRC, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results.rds")
saveRDS(age, file = "age.rds")
saveRDS(traj_design, file = "traj_design.rds")
saveRDS(traj_milo, file = "traj_milo.rds")

#Differential abundance analysis of PBMCs, between interstitial lung disease (ILD) without SRC patients and other SSc patients.
s.int_SSc$ILDwithoutSRC <- as.character(s.int_SSc$ILDwithoutSRC)
Idents(s.int_SSc) <- "ILDwithoutSRC"
subsampled_cells <- s.int_SSc@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  sample_n(1000) %>%
  pull(cell_id)
s.int_SSc_1000 <- subset(s.int_SSc, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(s.int_SSc_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.5, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILDwithoutSRC")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILDwithoutSRC")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,56,82,72,64,60,80,58,47,76,81,57,46,61,50,62,77,67,53,74,73)) #age of all SSc patients
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+ILDwithoutSRC, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results2.rds")
saveRDS(age, file = "age2.rds")
saveRDS(traj_design, file = "traj_design2.rds")
saveRDS(traj_milo, file = "traj_milo2.rds")


#Differential abundance analysis of monocytes, between SSc patients with SRC and those without SRC.
Mo[["seurat_clusters"]] <- Idents(Mo)
Mo$SRC <- as.character(Mo$SRC)
Idents(Mo) <- "disease"
Mo_SSc <- subset(Mo, ident = "SSc")
Idents(Mo_SSc) <- "SRC"
subsampled_cells <- Mo_SSc@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
Mo_SSc_1000 <- subset(Mo_SSc, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(Mo_SSc_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,56,82,72,64,60,80,58,47,76,81,57,46,61,50,62,77,67,53,74,73)) #age of all SSc patients
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+SRC, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results3.rds")
saveRDS(age, file = "age3.rds")
saveRDS(traj_design, file = "traj_design3.rds")
saveRDS(traj_milo, file = "traj_milo3.rds")

#Differential abundance analysis of monocytes, between SSc patients with SRC and healthy donors.
Idents(Mo) <- "SRC"
Mo_SRC_HD <- subset(Mo, ident = c("yes","HD"))
subsampled_cells <- Mo_SRC_HD@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
Mo_SRC_HD_1000 <- subset(Mo_SRC_HD, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(Mo_SRC_HD_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "SRC")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,72,53,74,64,62,70,67,59,84)) #age of SRC patients and healthy donors
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+SRC, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results4.rds")
saveRDS(age, file = "age4.rds")
saveRDS(traj_design, file = "traj_design4.rds")
saveRDS(traj_milo, file = "traj_milo4.rds")

#Differential abundance analysis of CD8+ T cells, between SSc patients with interstitial lung disease (ILD) and those without ILD.
CD8[["seurat_clusters"]] <- Idents(CD8)
CD8$ILD <- as.character(CD8$ILD)
Idents(CD8) <- "disease"
CD8_SSc <- subset(CD8, ident = "SSc")
Idents(CD8_SSc) <- "ILD"
subsampled_cells <- CD8_SSc@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
CD8_SSc_1000 <- subset(CD8_SSc, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(CD8_SSc_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,56,82,72,64,60,80,58,47,76,81,57,46,61,50,62,77,67,53,74,73)) #age of all SSc patients
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+ILD, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results5.rds")
saveRDS(age, file = "age5.rds")
saveRDS(traj_design, file = "traj_design5.rds")
saveRDS(traj_milo, file = "traj_milo5.rds")

#Differential abundance analysis of CD8+ T cells, between SSc patients with ILD and healthy donors.
Idents(CD8) <- "ILD"
CD8_ILD_HD <- subset(CD8, ident = c("yes","HD"))
subsampled_cells <- CD8_ILD_HD@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
CD8_ILD_HD_1000 <- subset(CD8_ILD_HD, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(CD8_ILD_HD_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,82,72,60,80,58,76,81,46,61,62,77,64,62,70,67,59,84)) #age of SSc-ILD patients and healthy donors
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+ILD, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results6.rds")
saveRDS(age, file = "age6.rds")
saveRDS(traj_design, file = "traj_design6.rds")
saveRDS(traj_milo, file = "traj_milo6.rds")

#Differential abundance analysis of CD4+ T cells, between SSc patients with interstitial lung disease (ILD) and those without ILD.
CD4[["seurat_clusters"]] <- Idents(CD4)
CD4$ILD <- as.character(CD4$ILD)
Idents(CD4) <- "disease"
CD4_SSc <- subset(CD4, ident = "SSc")
Idents(CD4_SSc) <- "ILD"
subsampled_cells <- CD4_SSc@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
CD4_SSc_1000 <- subset(CD4_SSc, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(CD4_SSc_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,56,82,72,64,60,80,58,47,76,81,57,46,61,50,62,77,67,53,74,73)) #age of all SSc patients
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+ILD, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results7.rds")
saveRDS(age, file = "age7.rds")
saveRDS(traj_design, file = "traj_design7.rds")
saveRDS(traj_milo, file = "traj_milo7.rds")

#Differential abundance analysis of CD4+ T cells, between SSc patients with ILD and healthy donors.
Idents(CD4) <- "ILD"
CD4_ILD_HD <- subset(CD4, ident = c("yes","HD"))
subsampled_cells <- CD4_ILD_HD@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  mutate(total_cells = n()) %>%
  group_modify(~ if (.x$total_cells[1] > 1000) sample_n(.x, 1000) else .x) %>%
  pull(cell_id)
CD4_ILD_HD_1000 <- subset(CD4_ILD_HD, cells = subsampled_cells)
pbmc_small_sce <- as.SingleCellExperiment(CD4_ILD_HD_1000)
traj_milo <- Milo(pbmc_small_sce)
gc()

traj_milo <- buildGraph(traj_milo, k = 5, d = 30, reduced.dim="REF.SPCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 5, d=30,reduced_dim="REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="REF.SPCA")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- data.frame(colData(traj_milo))[,c("sample", "ILD")]
traj_design <- distinct(traj_design)

age <- data.frame(c(67,82,72,60,80,58,76,81,46,61,62,77,64,62,70,67,59,84)) #age of SSc-ILD patients and healthy donors
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
gc()

da_results <- testNhoods(traj_milo,reduced.dim="REF.SPCA", design = ~ age+ILD, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)

nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,layout="UMAP",alpha=0.2) + scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))
nh_graph_pl

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters",alpha=0.2) + geom_boxplot(outlier.shape=NA, aes(group = seurat_clusters)) + scale_color_gradient2(midpoint=0,low="blue",mid="white",high="red",space="Lab",limits=c(-6,6)) + geom_hline(yintercept=0, linetype="dashed")
pDa

# Save milo data
saveRDS(da_results, file = "da_results8.rds")
saveRDS(age, file = "age8.rds")
saveRDS(traj_design, file = "traj_design8.rds")
saveRDS(traj_milo, file = "traj_milo8.rds")