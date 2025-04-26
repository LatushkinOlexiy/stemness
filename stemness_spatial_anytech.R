### Applying Stemness on Spatial datasets (Seurat objects) ###

# AKOYA/COSMX/XENIUM/VISIUM ... ---------
library(Seurat)
library(dplyr)
library(data.table)

# Create stemness function
stemness_akoya_seurat <- function(obj, ML_model_weight) {
  w <- ML_model_weight$w
  
  matrix_exp <- obj@assays$Akoya@data # You may want to consider using scale.data if you are experiencing challenges visualizing the stemness in your dataset.
  
  # Check the gene label between the model and the matrix_prediction
  length(intersect(rownames(matrix_exp),names(w))) 
  
  # Filter matrix_prediction by stemness model genes
  matrix_exp <- matrix_exp[rownames(matrix_exp)%in%  names(w) ,]
  length(rownames(matrix_exp))  
  
  # Filter stemness model with 'predict.DATA'
  w <- w[ rownames(matrix_exp) ]
  length(intersect(names(w),rownames(matrix_exp)))  
  length(names(w)) 
  is.vector(w) #TRUE
  
  # Score the predict.DATA using Spearman correlation
  s <- apply( matrix_exp, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
  s[1:5]
  
  # Scale the scores to be between 0 and 1
  s <- s - min(s)
  s <- s / max(s)
  s[1:5]
  
  # Assign stemness to seurat object
  s <- as.data.frame(s)
  obj[["stemness"]] <- s$s
  return(obj)
}

# Load the model
load("./model_RNA_MALTA.2018.Rda") # mm = stemness model

# Load your seurat object
obj_s185464 <- akoya_qupath_seurat_workflow('./measurements_S185464.tsv', fov = 'S185464')

# Run stemness function 
obj_s185464 <- stemness_akoya_seurat(obj_s185464, mm)





