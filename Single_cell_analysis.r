library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data = Read10X(data.dir = "C:/Users/HP/OneDrive/Documents/New folder")
pbmc.data=pbmc.data$`Gene Expression`
# Initialize the Seurat object with the raw (non-normalized data).
pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc
pbmc.data[1:50, 1:10]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^mt-")
head(pbmc@meta.data)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
pbmc = subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 1350 & percent.mt < 5)
pbmc
pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1350)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA@scale.data[1:50, 1:5]
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)
#pbmc = FindClusters(pbmc, resolution = 0)
#pbmc = FindClusters(pbmc, resolution = 1)

head(pbmc@meta.data)
pbmc = RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", label = T)
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)
#a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#a
filtered_markers = pbmc.markers %>% 
  filter(p_val_adj < 0.01, avg_log2FC > 0.1)
cluster_genes_list <- vector("list", length(unique(filtered_markers$cluster)))

# Populate the list with gene names for each cluster
for (i in seq_along(cluster_genes_list)) {
  cluster <- unique(filtered_markers$cluster)[i]
  genes <- filtered_markers$gene[filtered_markers$cluster == cluster]
  cluster_genes_list[[i]] <- genes
}

# Determine the maximum number of genes in a cluster
max_genes <- max(sapply(cluster_genes_list, length))

# Fill in shorter lists with NA to match max_genes
cluster_genes_list <- lapply(cluster_genes_list, function(x) {
  length_diff <- max_genes - length(x)
  c(x, rep(NA, length_diff))
})

# Create a data frame with gene names for each cluster
cluster_genes_df <- as.data.frame(do.call(cbind, cluster_genes_list))

cluster_to_celltype <- c("0" = "macrophages",
           "1" = "macrophages",
           "2" = "macrophages",
           "3" = "macrophages",
           "4" = "T-cells",
           "5" = "Erythroid cells",
           "6" = "Erythroid cells",
           "7" = "Erythroid cells",
           "8" = "T-cells",
           "9" = "T-cells")

# Assign cell types to clusters and add a new column 'cell_type' to metadata
pbmc@meta.data$cell_type <- cluster_to_celltype[as.character(pbmc@meta.data$seurat_clusters)]

# Now the 'cell_type' column is added to the metadata with assigned cell types
print(pbmc@meta.data)

# Create a UMAP plot grouped by cell type
DimPlot(pbmc, group.by = "cell_type")

