#Principal component analysis (PCA)

library(FactoMineR)
library(Factoshiny)
library(missMDA)
library(FactoInvestigate)
library("factoextra")
library(mice)
library(ggplot2)
library(ggrepel)

#In the factominer1.csv file, the ratio of each cell subpopulation to the total PBMC is included. The first column contains the sample ID and the first row contains the cell subpopulation name. The second column contains the name of the disease group to which each sample belongs (SRC, ILD without SRC, No SRC No ILD, or HD).
data <- read.csv("your_file_directory/factominer1.csv", header=TRUE, row.names=1)

# Perform PCA
pca_data <- data[, -1]
pca_result <- PCA(pca_data, scale.unit=TRUE, ncp=10, axes = c(1,2))
plot(pca_result, choix="ind", axes=c(1,2))
variable_groups <- rep("group1", ncol(pca_data))
variables_loadings <- get_pca_var(pca_result)$coord
variables_loadings <- as.data.frame(variables_loadings)
variables_loadings$Dim1_Dim2_length <- sqrt(variables_loadings[["Dim.1"]]^2 + variables_loadings[["Dim.2"]]^2)

#Myeloid cells, Lymphoid cells
wanted_elements <- c("CD14.VCAN", "CD14.PLBD1", "CD16", "CD14.ISG", "CD16.ISG", "cDC2", "intermediate", "CD14.HLA", "CD14.EGR1", "cDC1","pDC","pDC.ISG","ASDC")
indexes <- match(wanted_elements, colnames(pca_data))
indexes <- indexes[!is.na(indexes)]
variable_groups[indexes] <- "Myeloid"
wanted_elements2 <- c("CD4.CTL", "CD4.Naive", "CD4.TEM", "Treg", "CD4.ISG", "CD4.TCM_1", "CD4.TCM_2","CD8.CTL", "CD8.naive", "CD8.TEM", "CD8.TEM_T2ISG", "CD8.TCM", "CD8.TCM_GATA3","dnT","MAIT","gdT","B.Naive", "PB", "ABC", "B.Memory_post.switched", "B.Memory_pre.switched", "B.Naive_Activated", "Plasma.cell","NK.bright","NK.intermediate", "NK.dim", "ILC","Proliferating")
indexes2 <- match(wanted_elements2, colnames(pca_data))
indexes2 <- indexes2[!is.na(indexes2)]
variable_groups[indexes2] <- "Lymphoid"

arrow_size <- 0.8
fviz_pca_var(pca_result, col.var = variable_groups, 
             palette = c("Lymphoid" = "orange", "Myeloid" = "red"),
             arrowsize = arrow_size,
             legend.title = "Variable Groups",
             label = , addEllipses = FALSE) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +
  theme(
    panel.grid.major.x = element_line(size = 0.7, linetype = 'dotted', color = "black"),  # X軸の補助線
    panel.grid.major.y = element_line(size = 0.7, linetype = 'dotted', color = "black"),  # Y軸の補助線
    panel.grid.minor = element_blank()  # 細かい補助線を非表示
  )

groups <- data[, 1]
pca_data <- data[, -1]
pca_result <- PCA(pca_data, scale.unit=TRUE, ncp=5)
pca_scores <- as.data.frame(pca_result$ind$coord[, c(1, 2)])
names(pca_scores) <- c("Dim1", "Dim2")
pca_scores$Group <- groups

desired_order <- c("SRC", "ILD without SRC", "NO SRC No ILD", "HD")
pca_scores$Group <- factor(pca_scores$Group, levels = desired_order)
group_colors <- c("HD" = "black", "No SRC No ILD" = "deepskyblue3", "ILD without SRC" = "orange", "SRC" = "red")

contributions <- pca_result$eig[, 2]
xlab <- paste("PC1 (", round(contributions[1], 2), "%)", sep="")
ylab <- paste("PC2 (", round(contributions[2], 2), "%)", sep="")
group_means <- aggregate(cbind(Dim1, Dim2) ~ Group, data = pca_scores, FUN = mean)

ggplot(pca_scores, aes(x = Dim1, y = Dim2, color = Group)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.7) +  # 先にX軸の主軸を描画
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +  # 先にY軸の主軸を描画
  geom_point(size = 5) +  # 点をプロット
  geom_segment(data = group_means, aes(xend = Dim1, yend = Dim2, x = 0, y = 0, color = Group), arrow = arrow(), size = 1.5) +  # 矢印を追加
  scale_color_manual(values = group_colors) +
  labs(x = xlab, y = ylab, color = "Group") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(size = 0.7, linetype = 'dotted', colour = "black"), # X軸の補助線
    panel.grid.major.y = element_line(size = 0.7, linetype = 'dotted', colour = "black"), # Y軸の補助線
    panel.grid.minor.x = element_blank(), # X軸の細い補助線を消去
    panel.grid.minor.y = element_blank(), # Y軸の細い補助線を消去
    axis.text.x = element_text(size = 10, color = "black"), # X軸の数値のフォントサイズ
    axis.text.y = element_text(size = 10, color = "black") # Y軸の数値のフォントサイズ
  )

###############
#In the factominer2.csv file, the ratio of each cell subpopulation to the total PBMC is included. The first column contains the sample ID and the first row contains the cell subpopulation name. The second column contains the name of the disease group to which each sample belongs (SSc or HD).
# Read data from CSV file
data <- read.csv("your_file_directory/factominer2.csv", header=TRUE, row.names=1)

# Perform PCA
pca_data <- data[, -1]
pca_result <- PCA(pca_data, scale.unit=TRUE, ncp=10, axes = c(1,2))
plot(pca_result, choix="ind", axes=c(1,2))

variable_groups <- rep("group1", ncol(pca_data))

variables_loadings <- get_pca_var(pca_result)$coord
variables_loadings <- as.data.frame(variables_loadings)
variables_loadings$Dim1_Dim2_length <- sqrt(variables_loadings[["Dim.1"]]^2 + variables_loadings[["Dim.2"]]^2)

#Myeloid cells, Lymphoid cells
wanted_elements <- c("CD14.VCAN", "CD14.PLBD1", "CD16", "CD14.ISG", "CD16.ISG", "cDC2", "intermediate", "CD14.HLA", "CD14.EGR1", "cDC1","pDC","pDC.ISG","ASDC")
indexes <- match(wanted_elements, colnames(pca_data))
indexes <- indexes[!is.na(indexes)] # NA値を除外
variable_groups[indexes] <- "Myeloid"
wanted_elements2 <- c("CD4.CTL", "CD4.Naive", "CD4.TEM", "Treg", "CD4.ISG", "CD4.TCM_1", "CD4.TCM_2","CD8.CTL", "CD8.naive", "CD8.TEM", "CD8.TEM_T2ISG", "CD8.TCM", "CD8.TCM_GATA3","dnT","MAIT","gdT","B.Naive", "PB", "ABC", "B.Memory_post.switched", "B.Memory_pre.switched", "B.Naive_Activated", "Plasma.cell","NK.bright","NK.intermediate", "NK.dim", "ILC","Proliferating")
indexes2 <- match(wanted_elements2, colnames(pca_data))
indexes2 <- indexes2[!is.na(indexes2)] # NA値を除外
variable_groups[indexes2] <- "Lymphoid"

arrow_size <- 0.8
fviz_pca_var(pca_result, col.var = variable_groups, 
             palette = c("Lymphoid" = "orange", "Myeloid" = "red"),
             arrowsize = arrow_size,
             legend.title = "Variable Groups",
             label = , addEllipses = FALSE) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +
  theme(
    panel.grid.major.x = element_line(size = 0.7, linetype = 'dotted', color = "black"),  # X軸の補助線
    panel.grid.major.y = element_line(size = 0.7, linetype = 'dotted', color = "black"),  # Y軸の補助線
    panel.grid.minor = element_blank()  # 細かい補助線を非表示
  )

groups <- data[, 1]
pca_data <- data[, -1]
pca_result <- PCA(pca_data, scale.unit=TRUE, ncp=5)
pca_scores <- as.data.frame(pca_result$ind$coord[, c(1, 2)])
names(pca_scores) <- c("Dim1", "Dim2")
pca_scores$Group <- groups

desired_order <- c("SSc", "HD")
pca_scores$Group <- factor(pca_scores$Group, levels = desired_order)
group_colors <- c("HD" = "black", "SSc" = "mediumorchid3")

contributions <- pca_result$eig[, 2]
xlab <- paste("PC1 (", round(contributions[1], 2), "%)", sep="")
ylab <- paste("PC2 (", round(contributions[2], 2), "%)", sep="")
group_means <- aggregate(cbind(Dim1, Dim2) ~ Group, data = pca_scores, FUN = mean)

ggplot(pca_scores, aes(x = Dim1, y = Dim2, color = Group)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.7) +  # 先にX軸の主軸を描画
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.7) +  # 先にY軸の主軸を描画
  geom_point(size = 5) +  # 点をプロット
  geom_segment(data = group_means, aes(xend = Dim1, yend = Dim2, x = 0, y = 0, color = Group), arrow = arrow(), size = 1.5) +  # 矢印を追加
  scale_color_manual(values = group_colors) +
  labs(x = xlab, y = ylab, color = "Group") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(size = 0.7, linetype = 'dotted', colour = "black"), # X軸の補助線
    panel.grid.major.y = element_line(size = 0.7, linetype = 'dotted', colour = "black"), # Y軸の補助線
    panel.grid.minor.x = element_blank(), # X軸の細い補助線を消去
    panel.grid.minor.y = element_blank(), # Y軸の細い補助線を消去
    axis.text.x = element_text(size = 10, color = "black"), # X軸の数値のフォントサイズ
    axis.text.y = element_text(size = 10, color = "black") # Y軸の数値のフォントサイズ
  )