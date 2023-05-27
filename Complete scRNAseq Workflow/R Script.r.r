library(Seurat)
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

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "E:/scRNA-seq/Clustering and Marker Identification")
str(pbmc.data)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
str(pbmc)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#Check mitochondrial RNA percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Check ribosomal RNA percentage
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
View(pbmc@meta.data)

#Note MT:Bar-codes with a low count depth, few detected genes, and a 
#high fraction of mitochondrial counts are indicative of cells whose cytoplasmic 
#mRNA has leaked out through a broken membrane, and thus, only mRNA located in the 
#mitochondria is still conserved

#Note RB: Here RB means ribosomal protein and RNA coding genes. High number of
# RB genes reflect metabolically active cells. Based on cell types 15-45% RB genes 
# can be mapped. Based on requirement these genes can be removed (i.e., these
#genes don't contribute to heterogeneity)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                           'percent.rb'), ncol = 4,
        cols = '#C71000')


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = '#00468B', plot.cor = T, jitter = T, pt.size = 1)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = '#C71000', plot.cor = T, jitter = T, pt.size = 1) 
plot3 <- FeatureScatter(pbmc, feature1 = "percent.mt", feature2 = "percent.rb",
                        cols = 'darkgreen', plot.cor = T, jitter = T, pt.size = 1)

plot1 + plot2 + plot3 #Result is in Pearson-rank

#Discard cells having unique feature counts over 2,500 and less than 200
#Also filter out cells having >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Let's visualize again after QC 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                           "percent.rb"), ncol = 4, cols = '#C71000')

#Normalizing data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Alternative way that achives the same attribute as above method
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
# Scaling on large dataset may require extra time. An alternative way to 
#avoid the time lag by applying scaling only to the variable Features
#identified previously (i.e., 2000 in this case) by the following way
pbmc <- ScaleData(pbmc)

#Performing dimensionality reduction for PCA
pbmc <- RunPCA(object = pbmc)

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
           #disp.min = -10.5, #Negative cutoff
           #disp.max = 10.5 #Positive cutoff
           fast = F #Add legend (PC score scale)
           )
# Investigate the heterogeneity among different number of dimensions
DimHeatmap(pbmc, dims = 1:12, nfeatures = 20, cells = 500, balanced = TRUE)

#In particular DimHeatmap() allows for easy exploration of the primary sources 
#of heterogeneity in a dataset, and can be useful when trying to decide which PCs 
#to include for further downstream analyses. The impression from this dim plot
#is that PC1 and PC2 separetes the cell population into half whereas other PCs
#are not that comprehendible to understand heterogeneity like PC1 and PC2.

#Heatmap with manual settings
DimHeatmap(pbmc, dims = c(1,2,3,7:9), nfeatures = 20, cells = 500, balanced = TRUE)


# Alternative: Number of significant dimensions could be identified with an elbowplot
# Elbowplot utilizes a heuristic method to calculate SD among the PCs. Commonly used
El_plot<-ElbowPlot(pbmc)
El_plot<- El_plot+theme_bw()+theme(axis.text = element_text(size = 12))
El_plot
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

#Run Umap 
pbmc <- RunUMAP(pbmc, dims = 1:15)
###############################################################################
#Finding doublets
## pK Identification (no ground-truth)
sweep.res.list_pbmc <- paramSweep_v3(pbmc, PCs = 1:15, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() + theme_bw()+
  geom_line()
#The highest BCmetric corresponding to pK value represents the optimum pk value 
# (0.02 in this case)

## select the pK that corresponds to max bcmvn to perform doublet detection
pK <- bcmvn_pbmc %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate 
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
#Additional step for homotypic doublets since DoubletFinder is less sensitive
#to homotypic doublets
nExp_doub <- round(0.015*nrow(pbmc@meta.data))  

## Assuming 1.5% doublet in our dataset given 2700 cells (follow 10x genomics 
#user guide or DoubletFinder instructions)
nExp_doub.adj <- round(nExp_doub*(1-homotypic.prop))

# Run doubletFinder 
pbmc<- doubletFinder_v3(pbmc, 
                        PCs = 1:20, 
                        pN = 0.25, #default doesn't affect DblFnr performance
                        pK = pK, 
                        nExp = nExp_doub.adj,
                        reuse.pANN = FALSE, sct = FALSE) #No sct trans.

View(pbmc@meta.data)

#Numver of doublets and singlets
table(pbmc@meta.data$DF.classifications_0.25_0.02_33)

# visualize doublets
DimPlot(pbmc, reduction = 'umap', group.by = "DF.classifications_0.25_0.02_33",
        pt.size = 1) 
VlnPlot(pbmc, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.02_33")

#Remove doublets
pbmc = pbmc[, pbmc@meta.data[, "DF.classifications_0.25_0.02_33"] == "Singlet"]
###############################################################################

DimPlot(pbmc, group.by = "RNA_snn_res.1", reduction = "umap", label = T , pt.size=1)

#Customization
# Setting colors
color.length <- subset(pbmc, idents = 0:10)
levels(Idents(color.length ))
cols <- c('3'='#F68282','5'='#1FA195','1'='#B95FBB',
          '9'='#ff9a36','8'='#2FF18B','6'='#faf4cf',
          '2'='#CCB1F1','7'='#A4DFF2', '0'='#28CECA',
          '4'='#4B4BF7','10'='#E6C122')
cols <- cols[order(as.integer(names(cols)))]
scales::show_col(cols)

# Plot and finalize
umap<-DimPlot(pbmc, group.by="RNA_snn_res.1", reduction="umap", label=T , 
        pt.size=1, cols = cols)
umap<- umap+ theme_bw()+labs(x='UMAP 1', y='UMAP 2') + ggtitle('Resolution 1')+
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
umap

#TSNE Plot
pbmc<-RunTSNE(pbmc, dims = 1:10)
tsne<-DimPlot(pbmc, group.by = "RNA_snn_res.1", reduction = 'tsne', label = TRUE,
              cols=cols, pt.size=1)

#Customization
tsne<- tsne+theme_bw()+labs(x='tSNE 1', y='tSNE 2') + ggtitle('Resolution 1')+
       theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
tsne

#We can also visualize the percent of MT and RB genes in our datasets
#Make some palletes
Mt.pal<- viridis(n = 5, option = "D") 
Rb.pal<- viridis(n = 60, option = "A") #N depends on percent allowed in dataset

#Plot and customize
Mt<-FeaturePlot(pbmc, reduction='umap', features='percent.mt', pt.size=1)#, cols = Mt.pal) 
Rb<-FeaturePlot(pbmc, reduction='umap', features='percent.rb', pt.size=1)#, cols = Rb.pal)

Mt<-Mt+theme_bw()+labs(x='UMAP 1', y='UMAP 2')+ ggtitle('Mitochondrial Gene Abundance (%)')+
    theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
Rb<-Rb+theme_bw()+labs(x='UMAP 1', y='UMAP 2')+ ggtitle('Ribosomal Gene Abundance (%)')+
    theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))

Mt+Rb

# Save RDS
saveRDS(pbmc, file = "PBMC.rds")

# Find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# Find markers for every cluster compared to all remaining cells, report only the positive
# ones
All.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Find markers with specific tests
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, 
                              test.use = "roc", only.pos = TRUE)

#Checking features across all clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), log = T)

# Checking raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#Check in the form of feature plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                               "FCGR3A", "LYZ", "PPBP", "CD8A"), pt.size = 1, ncol = 3) &
               scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Let's higlight top 10 genes found in all clusters
All.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene)+ NoLegend()+
         scale_fill_gradient(high = 'darkblue', low = '#FFB900') +
         theme(axis.text = element_text(size=6, face = 'bold'))

#Annotation cluster
pbmc<- readRDS('PBMC.rds')
View(pbmc@meta.data)
DimPlot(pbmc, reduction = 'umap')

#Set reference dataset
ref<- celldex::HumanPrimaryCellAtlasData() 
View(as.data.frame(colData(ref)))

#Get expression marix
pbmc.counts<- GetAssayData(pbmc, slot = 'counts')

#Get prediction matrix from SingleR
pred <- SingleR(test = pbmc.counts,
                ref = ref,
                labels = ref$label.main)
# Tidy and inspect
pred$pruned.labels<- gsub('_'," ", as.character(pred$pruned.labels))
table(pred$pruned.labels)

# Assign labels to count matrix
pbmc$annotation <- pred$pruned.labels[match(rownames(pbmc@meta.data), rownames(pred))]

#Final plotting
cols2 <- c('B cell'='#F68282','CMP'='#1FA195','Monocyte'='#B95FBB',
           'NK cell'='#ff9a36','Platelets'='#2FF18B','Pre-B cell CD34-'='#faf4cf',
           'Pro-B cell CD34+'='#4B4BF7','T cells'='darkred')

Final<-DimPlot(pbmc, reduction = 'umap', group.by = 'annotation', label = F, cols = cols2, pt.size = 1)
Final<- Final+ theme_bw()+labs(x='UMAP 1', y='UMAP 2') + ggtitle('Annotated Clusters')+
    theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))
Final

# Annotation diagnostics 
#Based on the scores within cells 
pred
pred$scores
plotScoreHeatmap(pred)

# Based on deltas across cells 
plotDeltaDistribution(pred)

# Comparing to unsupervised clustering
pheatmap.data <- table(Assigned=pred$labels, Clusters=pbmc$seurat_clusters)
pheatmap(log10(pheatmap.data+10), color = colorRampPalette(c('white','blue'))(10))


## Table for expected doublets

#Rate (%)|#Cells Loaded|Cells Recovered
0.40%-----|800-------|500-------
0.80%-----|1,600-----|1,000-------
1.60%-----|3,200-----|2,000-------
2.30%-----|4,800-----|3,000-------
3.10%-----|6,400-----|4,000-------
3.90%-----|8,000-----|5,000-------
4.60%-----|9,600-----|6,000-------
5.40%-----|11,200----|7,000-------
6.10%-----|12,800----|8,000-------
6.90%-----|14,400----|9,000-------
7.60%-----|16,000----|10,000------