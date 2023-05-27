library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)
library(DoubletFinder)

#Set directory and read data
setwd('E:/scRNA-seq/Clustering and Marker Identification')
pbmc<- readRDS('PBMC.rds')
ElbowPlot(pbmc)



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
table(pbmc@meta.data$DF.classifications_0.25_0.03_34)

# visualize doublets
DimPlot(pbmc, reduction = 'umap', group.by = "DF.classifications_0.25_0.03_34",
        pt.size = 1) 
VlnPlot(pbmc, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.03_34")

#Remove doublets
pbmc = pbmc[, pbmc@meta.data[, "DF.classifications_0.25_0.03_34"] == "Singlet"]

## Table for expected doublets

#Rate (%)|#Cells Loaded|Cells Recovered
0.40%-----|800-------|500-------
0.80%-----|1,600-----|1,000-------
1.60%-----|3,200-----|2,000-------
2.30%-----|4,800-----|3,000-------
3.10%-----|6,400-----|4,000-------
3.90%-----|8,000-----|5,000-------
4.60%-----|9,600-----|6,000-------
5.40%-----|11,200-----|7,000-------
6.10%-----|12,800-----|8,000-------
6.90%-----|14,400-----|9,000-------
7.60%-----|16,000-----|10,000------