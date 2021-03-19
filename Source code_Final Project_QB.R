library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(scmap)
library(SAFEclustering)
library(cidr)
library(SingleCellExperiment)
library(SC3)
library(scater)
#library(ascend)
library(TSCAN)
library("scPred")
library(singleCellNet)
library(mclust)
library(Wind) #wRI

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")


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

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

#Proof on concept, should generate a sankey plot that maps onto itself. Make sure all relevant libraries are loaded!
testbarcodes <- bm@meta.data
Testlabels <- bm@meta.data[["celltype.l2"]]
table = data.frame() #this table is important! Will be referenced throughout the code, as it contains the ground truth information. 
for (i in 1:nrow(testbarcodes)){
  table<-rbind(table,c(rownames(testbarcodes[i,]),Testlabels[i]))
}
plot(getSankey(table[,2],table[,2]))

#Sankey attempt on same data with different labels. 
testbarcodes2 <- bm@meta.data
Testlabels2 <- bm@meta.data[["celltype.l1"]]
table2 = data.frame()
for (i in 1:nrow(testbarcodes2)){
  table2<-rbind(table2,c(rownames(testbarcodes2[i,]),Testlabels2[i]))
}
colnames(table2)<- c("Barcodes","cell_type1")

plot(getSankey(table[,2],table2[,2]))

#SAFEClustering, used to generate seurat plots. 
#generates count matrix from seurat data
counts.matrix <- bm@assays[["RNA"]]@data
#confirms that we have some MxN matrix
dim(counts.matrix) #SC3 is a wash
cluster.results <- individual_clustering(inputTags = counts.matrix, 
                                         mt_filter = FALSE, 
                                         nGene_filter = FALSE, 
                                         SC3 = FALSE, 
                                         gene_filter = FALSE, 
                                         CIDR = FALSE, 
                                         nPC.cidr = 25, 
                                         Seurat = TRUE, 
                                         nPC.seurat = 25, 
                                         resolution = 0.65, 
                                         tSNE = FALSE, 
                                         dimensions = 3, 
                                         perplexity = 30, 
                                         SEED = 123) # Seurat resolution parameter is nested in here. 

cluster.results <- data.frame(cluster.results)
colnames(cluster.results) <- colnames(counts.matrix)
cluster_comp_table<- t(cluster.results)
plot(getSankey(table[,2],cluster_comp_table[,1]))
adjustedRandIndex(table[,2],cluster_comp_table[,1]) #ARI generated here. Can be easily automated to generate multiple clusterings based on seurat resolution. 

#scmapcluster analysis begins here. 
colnames(table)<- c("Barcodes","cell_type1")

table <- as.data.frame(table)
#iterate through here
sce <- SingleCellExperiment(
  assays = list(
    normcounts = as.matrix(2^(counts.matrix-1))),
  colData=table)
logcounts(sce) = as.matrix(counts.matrix)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- selectFeatures(sce, suppress_plot = FALSE)
sce <- indexCluster(sce)
#Threshold parameter changes clustering fidelity.
scmapCluster_results <- scmapCluster(
  projection = sce, 
  index_list = list(
    clust = metadata(sce)$scmap_cluster_index
  ),threshold = 0.3
)
head(scmapCluster_results$scmap_cluster_labs)
scmap_results = data.frame()
for (i in 1:nrow(testbarcodes2)){
  scmap_results<-rbind(scmap_results,c(table2[i,1],scmapCluster_results$scmap_cluster_labs[i]))
}
plot(getSankey(table[,2],scmap_results[,2]))
adjustedRandIndex(table[,2],scmap_results[,2]) #CITE vs. scmap

#scmapcell analyses start here. Counts matrix is loaded into a single cell experiment. 
set.seed(1)
sce <- SingleCellExperiment(
  assays = list(
    normcounts = as.matrix(2^(counts.matrix-1))),
  colData=table)
logcounts(sce) = as.matrix(counts.matrix)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- selectFeatures(sce, suppress_plot = FALSE)
sce <- indexCell(sce)
names(metadata(sce)$scmap_cell_index)
dim(metadata(sce)$scmap_cell_index$subclusters)
metadata(sce)$scmap_cell_index$subclusters[1:5,1:5]
scmapCell_results <- scmapCell(
  sce, 
  list(
    clust = metadata(sce)$scmap_cell_index
  )
)
names(scmapCell_results)
scmapCell_results$clust$cells[,1:3]
ari_vals = vector()
thresholds = seq(0,1,0.01) #for generating sigmoid curve that shows ARI as a function of resolution. 
for (i in 1:length(thresholds)){
  scmapCell_clusters <- scmapCell2Cluster(
    scmapCell_results, 
    list(
      as.character(colData(sce)$cell_type1)
    ),threshold=thresholds[i]
  )
  head(scmapCell_clusters$scmap_cluster_labs)
  scmap_results = data.frame()
  ari_vals <- append(ari_vals,adjustedRandIndex(table[,2],scmapCell_clusters$scmap_cluster_labs)) #CITE vs. scmap
}
plot(getSankey(table[,2],scmapCell_clusters$scmap_cluster_labs))

#scmapcluster - General Schema for holding out data in scmap. 
set.seed(1)
#set holdout parameter, somewhere between 0 and 1
holdout_train = 1
holdout_test = 1
holdout_train = seq(.1,1,0.3) #alternatively set a range for iteration. If mapping single value, simply comment out the seq lines. 
holdout_test = seq(.1,1,0.3)
ARI_scmap_cluster = vector()
for (i in 1:length(holdout_train)){
  for (r in 1:length(holdout_test)){
    #transpose data such that rows are cells
    origmatrix <- as.data.frame(t(as.matrix(counts.matrix)))
    origmatrix <- cbind(origmatrix,table[,2])
    #randomly subset this matrix
    trainmatrix <- origmatrix[sample(nrow(origmatrix),holdout_train*nrow(origmatrix)), ] #put indeces back!
    set.seed(20)
    testmatrix <-  origmatrix[sample(nrow(origmatrix),holdout_test*nrow(origmatrix)), ]
    train_table <- cbind(rownames(trainmatrix),trainmatrix[,17010])
    colnames(train_table)<- c("Barcodes","cell_type1")
    test_table <- cbind(rownames(testmatrix),testmatrix[,17010])
    colnames(test_table)<- c("Barcodes","cell_type1")
    train_table <- as.data.frame(train_table)
    test_table <- as.data.frame(test_table)
    #keep the above matrices as they contain the original labels in the last column. Feed this to all subsequent analyses, including the sankey and ARI. 
    
    sce_train <- SingleCellExperiment(
      assays = list(
        normcounts = as.matrix(2^(t(trainmatrix[,1:17009])-1))),
      colData=train_table)
    logcounts(sce_train) = as.matrix(t(trainmatrix[,1:17009]))
    rowData(sce_train)$feature_symbol <- rownames(sce_train)
    
    sce_test <- SingleCellExperiment(
      assays = list(
        normcounts = as.matrix(2^(t(testmatrix[,1:17009])-1))),
      colData=test_table)
    logcounts(sce_test) = as.matrix(t(testmatrix[,1:17009]))
    rowData(sce_test)$feature_symbol <- rownames(sce_test)
    
    sce_train <- selectFeatures(sce_train, suppress_plot = FALSE)
    
    sce_train <- indexCluster(sce_train)
    
    scmapCluster_results <- scmapCluster(
      projection = sce_test, 
      index_list = list(
        clust = metadata(sce_train)$scmap_cluster_index
      ),threshold = 0.3
    )
    head(scmapCluster_results$scmap_cluster_labs)
    
    adjustedRandIndex(testmatrix[,17010],scmapCluster_results$scmap_cluster_labs)#CITE vs. scmap
    ARI_scmap_cluster <- append(ARI_scmap_cluster,adjustedRandIndex(testmatrix[,17010],scmapCluster_results$scmap_cluster_labs))
  }}

#For generating heatmap. Here, test is a matrix of all ARI values across the tested range. 
plot(as.matrix(test),breaks = range(test),axis.col=NULL,axis.row = NULL, xlab='',ylab='',key=list(side=2, cex.axis=1),digits = 4,text.cell=list(cex=3),fmt.cell='%.2f')
#to generate this matrix, an example of the code is shown below (ari_matrix):
holdout_train = seq(.05,1,0.1)
holdout_test = seq(.05,1,0.1)
ari_matrix <- matrix(,nrow=length(holdout_train),ncol=length(holdout_test))
ARI_scmap_cluster = vector()
for (i in 1:length(holdout_train)){
  for (r in 1:length(holdout_test)){
    print(i)
    print(r)
    #transpose data such that rows are cells
    origmatrix <- as.data.frame(t(as.matrix(counts.matrix)))
    origmatrix <- cbind(origmatrix,table[,2])
    #randomly subset this matrix
    trainmatrix <- origmatrix[sample(nrow(origmatrix),holdout_train[i]*nrow(origmatrix)), ]
    set.seed(20)
    testmatrix <-  origmatrix[sample(nrow(origmatrix),holdout_test[r]*nrow(origmatrix)), ]
    train_table <- cbind(rownames(trainmatrix),trainmatrix[,17010])
    colnames(train_table)<- c("Barcodes","cell_type1")
    test_table <- cbind(rownames(testmatrix),testmatrix[,17010])
    colnames(test_table)<- c("Barcodes","cell_type1")
    train_table <- as.data.frame(train_table)
    test_table <- as.data.frame(test_table)
    #keep the above matrices as they contain the original labels in the last column. Feed this to all subsequent analyses, including the sankey and ARI. 
    
    sce_train <- SingleCellExperiment(
      assays = list(
        normcounts = as.matrix(2^(t(trainmatrix[,1:17009])-1))),
      colData=train_table)
    logcounts(sce_train) = as.matrix(t(trainmatrix[,1:17009]))
    rowData(sce_train)$feature_symbol <- rownames(sce_train)
    
    sce_test <- SingleCellExperiment(
      assays = list(
        normcounts = as.matrix(2^(t(testmatrix[,1:17009])-1))),
      colData=test_table)
    logcounts(sce_test) = as.matrix(t(testmatrix[,1:17009]))
    rowData(sce_test)$feature_symbol <- rownames(sce_test)
    
    sce_train <- selectFeatures(sce_train, suppress_plot = FALSE)
    
    sce_train <- indexCluster(sce_train)
    
    scmapCluster_results <- scmapCluster(
      projection = sce_test, 
      index_list = list(
        clust = metadata(sce_train)$scmap_cluster_index
      ),threshold = 0.3
    )
    head(scmapCluster_results$scmap_cluster_labs)
    
    adjustedRandIndex(testmatrix[,17010],scmapCluster_results$scmap_cluster_labs)#CITE vs. scmap
    # ARI_scmap_cluster <- append(ARI_scmap_cluster,adjustedRandIndex(testmatrix[,17010],scmapCluster_results$scmap_cluster_labs))
    ari_matrix[i,r]= adjustedRandIndex(testmatrix[,17010],scmapCluster_results$scmap_cluster_labs)
  }}


