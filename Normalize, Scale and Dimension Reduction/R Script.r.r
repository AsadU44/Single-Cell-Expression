library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(reticulate)
reticulate::py_install(packages =
                          'umap-learn')

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "E:/scRNA-seq/Clustering and Marker Identification")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = '#C71000', plot.cor = T, jitter = T, pt.size = 1) 
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = '#00468B', plot.cor = T, jitter = T, pt.size = 1)
plot1 + plot2

#Discard cells having unique feature counts over 2,500 and less than 200
#Also filter out cells having >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Let's visualize again after QC 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalizing data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Alternative way
pbmc <- NormalizeData(pbmc)

#Identify top variable features (2000 default for a single dataset)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the dataset
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# Scaling on large dataset may require extra time. A better way to 
#avoid the time lag by applying scaling only to the variable Features
#identified previously (i.e., 2000 in this case) by the following way
pbmc <- ScaleData(pbmc)

#Performing dimensionality reduction for PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize different dimensions
VizDimLoadings(pbmc, 
               dims = 1:3, #Number of dimensions to display
               reduction = "pca",
               col = '#5773CC',
               nfeatures = 10, #Number of features to display
               combine = T, #Combine all PCs
               ncol = 3, #Number of columns in figure
               balanced = T # Split total n (features) to + and -
                )
#Visualize in the form of dimension plot
DimPlot(pbmc, reduction = "pca",
        pt.size = 2, #Size of the dots
        #group.by = 'orig.ident' # Group by class
        cols = '#9632B8')

# Investigate the heterogeneity among cells
DimHeatmap(pbmc, dims = 1,
           nfeatures = 30, #Number of top genes based on PCA scores 
           cells = 500, #Number of top cells based on PCA scores
           balanced = TRUE, #Equal number of top genes based on both +&- PCA scores
           #disp.min = -10.5,
           #disp.max = 10.5
           )
# Investigate the heterogeneity among different number of dimensions
DimHeatmap(pbmc, dims = 1:6, nfeatures = 20, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = c(1,2,3,7:9), nfeatures = 20, cells = 500, balanced = TRUE)

# Interpreting heterogeneity among cells across different PCs is not feasible. The
# number of significant dimensions could be identified with an elbowplot
# Elbowplot utilizes a heuristic method to calculate SD among the PCs. Commonly used
ElbowPlot(pbmc)

#Elbowplot suggests that at around 10th PC, the hetereogeneity among the dimensions 
# identified drops significantly

#Another alternative way
#JackStraw method takes random identified dimendions and utilizes statistical model
#to assign p-value
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
#This usually takes time for larger dataset. Elbowplot is advised to use.

#Forming clusters
pbmc<- FindNeighbors(pbmc, dims = 1:15)

#Setting resolution
pbmc <- FindClusters(pbmc, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(pbmc@meta.data)
DimPlot(pbmc, group.by = "RNA_snn_res.1", label = TRUE)

#Umap plot
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = T)

# Save RDS
saveRDS(pbmc, file = "../pbmc_tutorial.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
