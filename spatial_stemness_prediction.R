# Title: Spatial stemness prediction script
# Author: Maycon Marção, Felipe Segato
# Contact: marcao.legatum@gmail.com

# Data downloaded from https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard



### Stemness prediction on Spatial Omics data (Visium example) ###
setwd("/media/hd/maycon/Elena_stemness_request")

install.packages('Seurat')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DropletUtils")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
BiocManager::install("hdf5r")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")

library(Seurat)
library(glmGamPoi)
library(DropletUtils)
library(tidyr)

options(future.globals.maxSize = 8 * 1024^3)







# untaring 
untar("C:/Users/admin/Desktop/Bioinformatics/Project/Dataset/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz")

untar("C:/Users/admin/Desktop/Bioinformatics/Project/Dataset/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz")


# Create a .h5

feature_bc_matrix <- Read10X(data.dir = "C:/Users/admin/Desktop/Bioinformatics/Project/Dataset/filtered_feature_bc_matrix")
write10xCounts("C:/Users/admin/Desktop/Bioinformatics/Project/Dataset/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5", feature_bc_matrix, type = "HDF5",
               genome = "GRCh38",
               version = "3", overwrite = TRUE,
               gene.id = rownames(feature_bc_matrix),
               gene.symbol = rownames(feature_bc_matrix))

#Loading .h5 into R
data.dir <- "C:/Users/admin/Desktop/Bioinformatics/Project/Dataset"
list.files(data.dir, recursive = TRUE)
spData <- Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)
plot2 <- SpatialFeaturePlot(spData, features = "nCount_Spatial") #+ theme(legend.position = "right")
print(plot2)
#Transforming the data before applying stemness model
spData <- SCTransform(spData, assay = "Spatial", verbose = FALSE) 
SpatialFeaturePlot(spData, features = c("HPCA"))
spData@assays$Spatial@data
spData@assays$Spatial@counts
spData <- NormalizeData(spData, assay = "Spatial")
spData <- ScaleData(spData, assay = "Spatial")

# Extracting gene matrix 
gene_matrix <- as.matrix(spData@assays$SCT$data)
dim(gene_matrix)
gene_matrix[1:10, 1:10]

View(Assays(spData))
spData@assays$SCT
spData <- SCTransform(spData, assay = "Spatial", verbose = TRUE, method = "glmGamPoi")

# Loading Stemness model 
load("./model_RNA_MALTA.2018.Rda")
load("C:/Users/admin/Desktop/Bioinformatics/Project/model_RNA_MALTA.2018.Rda")
w = mm$w
w[1:5]
length(w) #12953 (number of genes on the model)
# Filtering matrix expression by the genes on stemness model
matrix = gene_matrix
length(intersect(rownames(matrix), names(w))) #12259 
predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
length(rownames(predict.DATA)) # 12259
w = w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA))) # 12259 
w[1:5]
length(names(w)) # 12259
is.vector(w) #TRUE

# Scoring the Matrix `X` using Spearman correlation
s = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scaling the scores to be between 0 and 1
s = s - min(s)
s = s / max(s)
s[1:5]
s = as.data.frame(s)
names(s) = "stemness"

# Checking cells order before add stemness index into the data
cellbarcodes_seq_ori <- rownames(spData@meta.data)
cellbarcodes_seq_s <- rownames(s)

all(cellbarcodes_seq_ori %in% cellbarcodes_seq_s) #TRUE
identical(cellbarcodes_seq_ori, cellbarcodes_seq_s) #TRUE

# Creating a stemness index column in metadata
spData@meta.data$stemness <- s$stemness
SpatialFeaturePlot(spData, features = c("stemness"))

#clustering
spData_clusterised <- spData
spData_clusterised <- RunPCA(object = spData_clusterised) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:15)
DimPlot(spData_clusterised, reduction = "umap")
spData_clusterised$seurat_clusters
SpatialDimPlot(spData_clusterised, interactive=TRUE) + SpatialFeaturePlot(spData_clusterised, features = c("stemness"))
SpatialDimPlot(spData_clusterised, interactive=TRUE)
?FindClusters
?SpatialDimPlot
Assays(spData_clusterised)
markers <- FindAllMarkers(
  object = spData_clusterised,
  only.pos = TRUE,       # only return upregulated genes (optional)
  min.pct = 0.25,        # only test genes expressed in at least 25% of cells
  logfc.threshold = 0.25 # minimum log fold-change
)
#Viewing the data
View(subset(markers, markers$cluster==13))
spData_clusterised$seurat_clusters[1:10]
View(spData_clusterised)

#average stemness
mean(spData$stemness)
#median stemness
median(spData$stemness)
#amount of cells with high stemness (above 0.75)
length(subset(spData$stemness, spData$stemness > 0.75))
#amount of total cells
length(spData$stemness)
#ratio of high stemness cells to total cells
length(subset(spData$stemness, spData$stemness > 0.75))/length(spData$stemness)


     