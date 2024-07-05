#Monocle 3

library(Seurat)
library(Signac)
library(SeuratWrappers)
library(remotes)
library(monocle3)
library(ggplot2)
library(patchwork)

#"Mo" is a seurat object created by integration analysis of Abseq data derived from a SSc patient at the onset of SRC. Please refer to the script "Abseq_Kidney_Blood_cells_clustering.R".

DefaultAssay(Mo) <- "SCT"
pbmc.combined.cds <- as.cell_data_set(Mo)
pbmc.combined.cds <- estimate_size_factors(pbmc.combined.cds)
rowData(pbmc.combined.cds)$gene_short_name <- row.names(rowData(pbmc.combined.cds))
pbmc.combined.cds <- cluster_cells(cds = pbmc.combined.cds, reduction_method = "UMAP", resolution = 0.005)
pbmc.combined.cds@clusters@listData[["UMAP"]][["clusters"]]
pbmc.combined.cds <- learn_graph(pbmc.combined.cds, use_partition = TRUE)
plot_cells(pbmc.combined.cds, cell_size = 1)
# set root node.
pbmc.combined.cds <- order_cells(pbmc.combined.cds)
# result
plot_cells(pbmc.combined.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,trajectory_graph_color="red",trajectory_graph_segment_size = 1.5,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
