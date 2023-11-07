library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
library(Azimuth)
library(patchwork)
library(CellChat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(reticulate)
library(viridis)
library(ggmin)
library(DESeq2)
library(colorspace)
library(RColorBrewer)
library(ggpubr)
library(SingleR)
library(celldex)
library(pheatmap)
library(DoubletFinder)

#Set working directory
setwd('D:/CancerData/PDAC_scRNAseq')

#Read  patient data and add group
Patient1_Normal<- Read10X(data.dir = "D:/CancerData/PDAC_scRNAseq/Patient 1_Normal")
Patient1_Normal <- CreateSeuratObject(counts = Patient1_Normal, project = 'Normal')
Patient1_Normal$Group <- "Normal"

Patient1_PDAC<- Read10X(data.dir = "D:/CancerData/PDAC_scRNAseq/Patient 1_PDAC")
Patient1_PDAC <- CreateSeuratObject(counts = Patient1_PDAC, project = 'Cancer')
Patient1_PDAC$Group <- "Cancer"

#Combine dataset
PDAC.combined <- merge(Patient1_Normal, y = Patient1_PDAC, add.cell.ids = c("N", "C"), project = "PDAC")
PDAC.combined
head(colnames(PDAC.combined))

View(PDAC.combined@meta.data)

# create a sample column
PDAC.combined$sample <- rownames(PDAC.combined@meta.data)


# calculate mitochondrial and ribosomal percentage genes
PDAC.combined$mitoPercent <- PercentageFeatureSet(PDAC.combined, pattern='^MT-')
PDAC.combined$riboPercent <- PercentageFeatureSet(PDAC.combined, pattern='^RP[SL]')

# explore QC
VlnPlot(PDAC.combined, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent", 'riboPercent'), ncol = 4)


# filtering
merged_seurat_filtered <- subset(PDAC.combined, subset = nCount_RNA < 8000 &
                                   nFeature_RNA < 1500 &
                                   mitoPercent < 50)
#Plot after filtering
plot1 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = c('#C71000','#00468B'), plot.cor = T, jitter = T, pt.size = 1)
plot2 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "mitoPercent",
                        cols = c('#C71000','#00468B'), plot.cor = T, jitter = T, pt.size = 1)
plot3 <- FeatureScatter(merged_seurat_filtered, feature1 = "mitoPercent", feature2 = "riboPercent",
                        cols = c('#C71000','#00468B'), plot.cor = T, jitter = T, pt.size = 1)
plot1 + plot2 + plot3



# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:15)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:15)

#Show variable features
Variable_PDAC <- FindVariableFeatures(merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(Variable_PDAC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

###############################################################################
#Finding doublets
## pK Identification (no ground-truth)
sweep.res.list_PDAC <- paramSweep_v3(merged_seurat_filtered, PCs = 1:15, sct = FALSE)
sweep.stats_PDAC <- summarizeSweep(sweep.res.list_PDAC, GT = FALSE)
bcmvn_PDAC <- find.pK(sweep.stats_PDAC)

ggplot(bcmvn_PDAC, aes(pK, BCmetric, group = 1)) +
  geom_point() + theme_bw()+
  geom_line()

#The highest BCmetric corresponding to pK value represents the optimum pk value 
# (0.220 in this case)

## select the pK that corresponds to max bcmvn to perform doublet detection
pK <- bcmvn_PDAC %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate 
annotations <- merged_seurat_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

#Additional step for homotypic doublets since DoubletFinder is less sensitive
#to homotypic doublets
nExp_doub <- round(0.015*nrow(merged_seurat_filtered@meta.data))  

## Assuming 1.5% doublet in our dataset given 2700 cells (follow 10x genomics 
#user guide or DoubletFinder instructions)
nExp_doub.adj <- round(nExp_doub*(1-homotypic.prop))

# Run doubletFinder 
merged_seurat_filtered<- doubletFinder_v3(merged_seurat_filtered, 
                        PCs = 1:20, 
                        pN = 0.25, #default doesn't affect DblFnr performance
                        pK = pK, 
                        nExp = nExp_doub.adj,
                        reuse.pANN = FALSE, sct = FALSE) #No sct trans.

View(merged_seurat_filtered@meta.data)

#Numver of doublets and singlets
table(merged_seurat_filtered@meta.data$DF.classifications_0.25_0.24_59)

# visualize doublets
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.24_59",
        pt.size = 1)

VlnPlot(merged_seurat_filtered, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.24_59")

#Remove doublets
merged_seurat_filtered = merged_seurat_filtered[, merged_seurat_filtered@meta.data[, "DF.classifications_0.25_0.24_59"] == "Singlet"]
###############################################################################

# Examine and visualize PCA results a few different ways
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize different dimensions
jpeg('Plot_3.jpeg', width = 10, height = 7, units = 'in', res = 600)
VizDimLoadings(merged_seurat_filtered, 
               dims = 1:3, #Number of dimensions to display
               reduction = "pca",
               col = '#5773CC',
               nfeatures = 10, #Number of features to display
               combine = T, #Combine all PCs
               ncol = 3, #Number of columns in figure
               balanced = T # Split total n (features) to + and -
)
dev.off()

# Investigate the heterogeneity among cells
DimHeatmap(merged_seurat_filtered, dims = 1,
           nfeatures = 30, #Number of top genes based on PCA scores 
           cells = 500, #Number of top cells based on PCA scores
           balanced = TRUE, #Equal number of top genes based on both +&- PCA scores
           #disp.min = -10.5, #Negative cutoff
           #disp.max = 10.5 #Positive cutoff
           fast = FALSE #Add legend (PC score scale)
)


# Investigate the heterogeneity among different number of dimensions
DimHeatmap(merged_seurat_filtered, dims = 1:12, nfeatures = 20, cells = 500, balanced = TRUE,fast = F)



# plot
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Group', pt.size = 1)


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Group')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:30)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)
#Plot
DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Group', pt.size = 1)


# Visualization
p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Group", pt.size = 1)
p2 <- DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
p1 + p2


#Split Cluster According to Group
DimPlot(seurat.integrated, reduction = "umap", split.by = "Group", pt.size = 1)

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(seurat.integrated) <- "RNA"
nk.markers <- FindConservedMarkers(seurat.integrated, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers)

write.csv(nk.markers, file='Markers_NormalvsCancer.csv')

#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(seurat.integrated, features = c("PLAUR", "EREG", "VCAN", "TYMP", "NAMPT", "CD14" ), min.cutoff = "q9")



#The DotPlot() function with the split.by parameter can be useful for viewing 
#conserved cell type markers across conditions, showing both the expression 
#level and the percentage of cells in a cluster expressing any given gene. 
#Here we plot 2-3 strong marker genes for each of our 14 clusters.
markers.to.plot <- c("PLAUR", "EREG", "VCAN", "TYMP", "NAMPT", "CD14",'PHGDH', "PSAT1", 'NRF2', "ATF4", "KEAP1","PSPH")

DotPlot(seurat.integrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "Group") +
  RotatedAxis()


#Now that we've aligned the stimulated and control cells, we can start to do comparative 
#analyses and look at the differences induced by stimulation. One way to look broadly 
#at these changes is to plot the average expression of both the stimulated and control 
#cells and look for genes that are visual outliers on a scatter plot.

#Feature Plot
FeaturePlot(seurat.integrated, features = c("CD3D", "GNLY", "IFI6", 'ATF4'), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

#################################################
#Set reference dataset
ref<- celldex::HumanPrimaryCellAtlasData() 
View(as.data.frame(colData(ref)))

#Get expression marix
PDAC.counts<- GetAssayData(seurat.integrated, slot = 'counts')

#Get prediction matrix from SingleR
pred <- SingleR(test = PDAC.counts,
                ref = ref,
                labels = ref$label.main)
# Tidy and inspect
pred$pruned.labels<- gsub('_'," ", as.character(pred$pruned.labels))
table(pred$pruned.labels)

# Assign labels to count matrix
seurat.integrated$annotation <- pred$pruned.labels[match(rownames(seurat.integrated@meta.data), rownames(pred))]
DimPlot(seurat.integrated, reduction = 'umap', group.by = 'annotation', label = F, pt.size = 1)
View(seurat.integrated@meta.data)

#####################################################################
plots <- VlnPlot(seurat.integrated, features = c("PLAUR", 'ATF4'), split.by = "Group", group.by = "annotation",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)




