library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Rfast2)
library(celldex)
library(SingleR)

setwd('E:/scRNA-seq/Spatial scRNAseq/Dataset')

### Load dataset
lung.data<- Load10X_Spatial(data.dir = 'E:/scRNA-seq/Spatial scRNAseq/Dataset',
                            filename = 'raw_feature_bc_matrix.h5',
                            assay = "Spatial", # Name of the assay
                            slice = "Lung", # Name of the test image
                            filter.matrix = TRUE, 
                            to.upper = FALSE)
### Initial Assessment
#nFeature_Spatial: the number of unique genes in each sample
#nCount_Spatial: the total number of detected molecules in each sample
lung.data
#Inspect metadata
View(lung.data@meta.data)
#Inspect dimension
dim(lung.data) # 33601 features across 4992 samples
#Check feature names
head(rownames(lung.data), n = 5)
lung.data@assays$Spatial@counts[5:10, 1:3]

# Plot 
plot1 <- VlnPlot(lung.data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(lung.data, features = "nCount_Spatial") + 
        theme(legend.position = "right")

wrap_plots(plot1, plot2)

#Check mitochondrial RNA percentage
lung.data[["percent.mt"]] <- PercentageFeatureSet(lung.data, pattern = "^MT-")

#Check ribosomal RNA percentage
lung.data[["percent.rb"]] <- PercentageFeatureSet(lung.data, pattern = "^RP[SL]")
View(lung.data@meta.data)

## Assessing correlation between assay variables
count.v.feature <- FeatureScatter(lung.data, feature1 = "nCount_Spatial", 
                  cols = 'darkgreen', feature2 = "nFeature_Spatial") +NoLegend()
cout.v.mt <- FeatureScatter(lung.data, feature1 = "nCount_Spatial", 
             feature2 = "percent.mt", cols = 'darkred')+ NoLegend()

count.v.feature+cout.v.mt 

### Visual inspection of the distribution of the assay data before QC
raw.vp<- VlnPlot(lung.data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", 
                           'percent.rb'), ncol = 4,
                             cols = '#C71000')
raw.vp
# The % of Rb genes looks fine (i.e., below 45%)
raw.fp<- SpatialFeaturePlot( lung.data, features = c("nFeature_Spatial", "nCount_Spatial", 
                 "percent.mt")) & theme(legend.position = "bottom")  
raw.fp


### Performing QC
lung.qc <- subset(lung.data, subset = nFeature_Spatial < 5000 & percent.mt < 20
               & nCount_Spatial < 20000)


### Visual inspection of the distribution of the assay data before normalization
qc.vp<- VlnPlot(lung.qc, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", 
                  'percent.rb'), ncol = 4, cols = '#C71000')
qc.vp

#Combine and compare
wrap_plots(raw.vp, qc.vp, nrow = 2)

qc.fp<- SpatialFeaturePlot( lung.qc, features = c("nFeature_Spatial", "nCount_Spatial", 
                           "percent.mt")) & theme(legend.position = "bottom")  
qc.fp

wrap_plots(raw.fp, qc.fp, nrow = 2)

### Perform normalization
#In spatial scRNA-seq the heterogeneity across cells can't only be attributed to
# technical issues. Rather a lot of dissimilarities also account for the specific
# position of the cells. The developer suggests that using SCTransform rather than
# logNormalize performs well in terms of retaining location(tissue)-specific heterogeneity 
# across the samples
lung.norm<- SCTransform(lung.qc, assay = "Spatial", verbose = FALSE)
names(lung.norm)
rm(lung.data)
rm(lung.qc)

### Downstream Analysis
# Run dimensionality reduction with PCA
lung.norm <- RunPCA(lung.norm, assay = "SCT", verbose = FALSE)

# Select number of dimensions
ElbowPlot(lung.norm)

# Compute Shared nearest neighbors (SNN)
lung.norm <- FindNeighbors(lung.norm, reduction = "pca", dims = 1:15)

# Leiden algorithm for community detection
lung.norm <- FindClusters(lung.norm, verbose = FALSE)
View(lung.norm@meta.data)

# RUN UMAP, PCA result is the default UMAP input
lung.final <- RunUMAP(lung.norm, reduction = "pca", dims = 1:15)
rm(lung.norm)
#Plotting
plot3 <- DimPlot(lung.final, reduction = "umap", group.by = 'seurat_clusters',
                 label = TRUE, pt.size = 1) 
plot3<- plot3+theme_bw()+labs(x='UMAP 1', y='UMAP 2') + ggtitle('Seurat Clusters')+
        theme(axis.text = element_text(size = 12))+ NoLegend()

plot4 <- SpatialDimPlot(lung.final, label = TRUE, label.size = 3) 
plot4<- plot4+theme_bw()+labs(x='', y='')+ NoLegend()+ theme(axis.text 
        = element_blank(), axis.ticks = element_blank())
plot3 + plot4


#### Subset and plot
subset <- subset(lung.final, idents = c(1, 2, 7, 8, 12))

SpatialDimPlot(subset, pt.size.factor = 2)+theme_bw()+ labs(x='', y='')+  
      theme(axis.text = element_blank(), axis.ticks = element_blank())

### Find Markers
#find all markers of cluster 3
cluster3_markers <- FindMarkers(lung.final, ident.1 = 3, min.pct = 0.25)
head(cluster3_markers, n = 5)

#### Visualize markers in the form of violin plot
VlnPlot(lung.final, features = c('SCGB3A1', 'SLPI', 'BPIFA1'))

#### Visualize markers in the form of feature plot
SpatialFeaturePlot(lung.final, features = c('SCGB3A1', 'SLPI', 'RPS14'), pt.size.factor = 2)

### Find All Markers
lung.markers <- FindAllMarkers(lung.final, only.pos = TRUE, min.pct = 0.20, 
                                    logfc.threshold = 0.20)
lung.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top10 <- lung.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(lung.final, features = top10$gene) + NoLegend() +  
           scale_fill_gradient(high = 'darkblue', low = '#FFB900')

### Identify spatially variable genes
#### Method I: Moran's I

#Perform Moran's I
lung.moransi <- FindSpatiallyVariableFeatures(
  lung.final, assay = "SCT", 
  features = VariableFeatures(lung.final)[1:10], #We are considering only 10 genes
  #Since running against 17000 features is time consuming
  selection.method = "moransi") 

#Get dataframe and inspect
lung.moransi.result <- lung.moransi@assays$SCT@meta.features %>%
   na.exclude #NA value due to calculating moransi only for 10 genes
head(lung.moransi.result[order(lung.moransi.result$MoransI_observed, decreasing = T),])

### Plot 3 most variable genes according to Moran's I experiment
#According to higher Moran's I Value
Best_Moransi_Marker <- head(
  SpatiallyVariableFeatures(lung.moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(lung.moransi,  features = Best_Moransi_Marker, ncol = 3, alpha 
                   = c(0.1, 1), pt.size.factor = 2) + 
  plot_annotation(
    title = "Moran's I: Best 3 Markers That Are Spatially Distinct",
    subtitle = "Accoding to Higher Moran's I Score Rank")

#According to lower Moran's I Value
Least_Moransi_Marker <- tail(
  SpatiallyVariableFeatures(lung.moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(lung.moransi,  features = Least_Moransi_Marker, ncol = 3, alpha 
                   = c(0.1, 1), pt.size.factor = 2) + 
  plot_annotation(
    title = "Moran's I: Best 3 Markers That Are Spatially Distinct",
    subtitle = "Accoding to Lower Moran's I Score Rank")

#### Method II: Variogram
#Perform variogram
lung.variogram <- FindSpatiallyVariableFeatures(
  lung.final, assay = "SCT", 
  features = VariableFeatures(lung.final)[1:10], #Same 10 marker as MoransI
  selection.method = "markvariogram") 

#Get result
lung.variogram.reuslt<- lung.variogram @assays$SCT@meta.features %>%
  na.exclude 
head(lung.variogram.reuslt[order(lung.variogram.reuslt$r.metric.5), ])

#According to higher Moran's I Value
Best_Variogram_Marker <- head(
  SpatiallyVariableFeatures(lung.variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(lung.variogram,  features = Best_Variogram_Marker, ncol = 3, alpha 
                   = c(0.1, 1), pt.size.factor = 2) + 
  plot_annotation(
    title = "Variogram: Best 3 Markers That Are Spatially Distinct",
    subtitle = "Accoding to Higher Variogram Score Rank")

#According to lower Variogram Value
Least_Ranked_Marker <- tail(
  SpatiallyVariableFeatures(lung.variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(lung.variogram,  features = Least_Ranked_Marker, ncol = 3, alpha 
                   = c(0.1, 1), pt.size.factor = 2) + 
  plot_annotation(
    title = "Variogram: Best 3 Markers That Are Spatially Distinct",
    subtitle = "Accoding to Lower Variogram Score Rank")

### Annotation and Visualization
#We performed manual annotation from PanglaoDB and CellMarker 2.0 due to 
#shortage of memory and now renaming clusters
new.names<- c('0'='Immune Cell', '1'='Basal Cell', '2'='Immune Cell', 
              '3'='Airway Secretory Cell', '4'='Fibroblast', '5'='Immune Cell',
              '6'='Smooth Muscle Cell', '7'='Airway Secretory Cell', '8'=
                'Fibroblast', '9'='Smooth Muscle Cell', '10'='Fibroblast',
              '11'='Airway Secretory Cell', '12'='Epithelial Cell')
#Rename
lung.final@meta.data$Annotation<- new.names[lung.final@meta.data$seurat_clusters]
table(lung.final@meta.data$Annotation)
#Final plotting
cols2 <- c('Airway Secretory Cell'='#F68282','Basal Cell'='#1FA195','Fibroblast'='#B95FBB',
           'Immune Cell'='#ff9a36','Smooth Muscle Cell'='#4B4BF7','Epithelial Cell'='darkred')

#Final plotting
#Dimension plot
jpeg('Annotated.jpg', width = 20, height = 12, units = 'cm', res = 600)
Final.dim<-DimPlot(lung.final, reduction = 'umap', group.by = 'Annotation', label = F, cols = cols2, pt.size = 2)

Final.dim<- Final.dim+ theme_bw()+labs(x='UMAP 1', y='UMAP 2') + ggtitle('Annotated Clusters')+
  theme(axis.text = element_text(size = 12), legend.text = element_text(size=12))

#Spatial feature plot
Final.fet <- SpatialDimPlot(lung.final, group.by = 'Annotation', label = F,
                            pt.size.factor  = 3, cols = cols2) 
Final.fet<- Final.fet+theme_bw()+labs(x='', y='')+ theme(axis.text 
                       = element_blank(), axis.ticks = element_blank())+ NoLegend()

Final.dim+Final.fet

dev.off()
