#scPred Script
devtools::install_github("immunogenomics/harmony", force = TRUE)
devtools::install_github("powellgenomicslab/scPred")
library(scPred)
library(Seurat)
library(magrittr)
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)


# Load in multimodal citeseq data
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
# We first perform pre-processing and dimensional reduction on both assays independently. We use standard normalization, but you can also use SCTransform or any alternative method.

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)



# Load in RNA data alone from CITEseq to train the data set, separate line so I dont get confused
cbmc.rna <- as.sparse(read.csv(file = "/STEMREM205/Project/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "/STEMREM205/Project/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# Change default assay
Assays(cbmc)
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(bm) <- "RNA"

# Load in the scPred provided PBMC dataset
provided_reference <- scPred::pbmc_1


# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(bm, label = TRUE)


####Workflow for setting up training model####
# This sets up using RNA object from CITEseq to be used as reference
reference <- cbmc
reference
# Similar to clustering in `Seurat`, `scPred` uses the cell embeddings from a 
# principal component analysis to make inferences about cell-type identity. 
# however —unlike clustering—, `scPred` trains classifiers for each cell type of 
# interest in a supervised manner by using the known cell identity from a 
# reference dataset to guide the classification of cells in a different data set.


# Normalization
# Step 1: Normalizes the gene expression data from the reference data set weighting by the counts by the "sequencing depth" of each cell and applying a natural logarithmic transformation
# Step 2 Finds a set of variable features by modeling the mean and variance log-transformed expression
# Step 3 Scales the expression of each gene bu subtracting the mean expression across all cells and dividing by the standard deviation
# Step 4 Runs a PCA
# Step 5 Runs a UMAP using the top 30 most variable PCs
bm <- bm %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
# The column `cell_type` contains the identity of each cell in the meta data slot.
# Let's plot the UMAP and grouping the cells by cell type.
bm$celltype.l2
DimPlot(bm, group.by = "celltype.l2", label = TRUE, repel = TRUE)
# Training classifiers with `scPred`
#Firstly, let's get the feature space to train the classifiers. By default, `scPred` 
#will use all principal components. The reference labels of each cells are 
#specified as the second parameter value of the function (in this case the cell_type` column. 

# Training on the multimodal clusters
reference <- getFeatureSpace(bm , "celltype.l2")

# Secondly, we train the classifiers for each cell using the trainModel function. By default, scPred will use a support vector machine with a radial kernel.
reference <- trainModel(reference)
# Training probabilities for each cell in the reference data can be accessed using the get_probabilities method:
get_probabilities(reference) %>% head()
#We can use the get_scpred method to retrieve the scPred object from the Seurat object. Printing a scPred object will show for each cell type:
get_scpred(reference)
# To visualize the performance for each cell type we can use the plot_probabilities function:
plot_probabilities(reference)


# An important requirement for classifying cells is using the same normalization method for both the reference and the query datasets.
# First, let’s normalize the query dataset (cells to be classfied).
bm <- NormalizeData(bm)

# This is if you want to recluster the reference provided using the CITEseq training
provided_reference <- NormalizeData(provided_reference)


# classify the cells from the query data using the scPredict function. The first argument corresponds to the query object and the second to the reference object (with a scPred model trained already).
# scPred now uses Harmony to align the query data onto the training low-dimensional space used as reference. Once the data is aligned, cells are classified using the pre-trained models.
bm <- scPredict(bm, reference)
provided_reference <- scPredict(provided_reference, reference)
provided_reference$cell_type


#scPred will store the final classifications in the scpred_prediction column of the Seurat meta data. 
DimPlot(bm, group.by = "scpred_prediction", reduction = "scpred")
DimPlot(provided_reference, group.by = "scpred_prediction", reduction = "scpred")
bm <- RunUMAP(bm, reduction = "scpred", dims = 1:30)
provided_reference <- RunUMAP(provided_reference, reduction = "scpred", dims = 1:30)


# and plot the predicted labels for each cell type over the UMAP:
DimPlot(bm, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(provided_reference, group.by = "scpred_prediction", label = TRUE, repel = TRUE)


# compare with original labels
DimPlot(provided_reference, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(provided_reference, group.by = "cell_type", label = TRUE, repel = TRUE)


Dimplot(bm, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(provided_reference, group.by = "cell_type", label = TRUE, repel = TRUE)

# To verify the performance of the models in the query dataset, we can use the crossTab to create a contingency table using two colums from the metadata. 
# In this example, the cell type info is sotred in the cell_type columns and the predicted labels for each cell in the scpred_prediction column.
crossTab(bm, "cell_type", "scpred_prediction")

get_classifiers(reference)


####END####