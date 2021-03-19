library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(scmap)
library(remotes)
devtools::install_github("VCCRI/CIDR")
devtools::install_github("yycunc/SAFEclustering")
library(SAFEclustering)
library(cidr)
library(SingleCellExperiment)
library(SC3)
library(scater)
devtools::install_github("zji90/TSCAN")
library(TSCAN)
library(scPred)
devtools::install_github("zji90/TSCAN")
devtools::install_github("pcahan1/singleCellNet") 


#Proof on concept, should generate a sankey plot that maps onto itself. 
#Uses metadata/barcodes from bone marrow multimodal dataset and maps to itself
testbarcodes <- bm@meta.data
Testlabels <- bm@meta.data[["celltype.l2"]]
table = data.frame()
for (i in 1:nrow(testbarcodes)){
  table<-rbind(table,c(rownames(testbarcodes[i,]),Testlabels[i]))
}
plot(getSankey(table[,2],table[,2]))

#Sankey attempt on same data with different labels. 
#Uses metadata/barcodes from bone marrow multimodal dataset and maps onto scpred predicted clusters
testbarcodes2 <- bm@meta.data
Testlabels2 <- bm@meta.data[["scpred_prediction"]]
table2 = data.frame()
for (i in 1:nrow(testbarcodes2)){
  table2<-rbind(table2,c(rownames(testbarcodes2[i,]),Testlabels2[i]))
}
colnames(table2)<- c("Barcodes","scpred_prediction")
plot(getSankey(table[,2],table2[,2]))


####Using reference data####
#Proof on concept, should generate a sankey plot that maps onto itself. 
testbarcodes <- provided_reference@meta.data
Testlabels <- provided_reference@meta.data[["cell_type"]]
table = data.frame()
for (i in 1:nrow(testbarcodes)){
  table<-rbind(table,c(rownames(testbarcodes[i,]),Testlabels[i]))
}
plot(getSankey(table[,2],table[,2]))
#Sankey attempt on same data with different labels. 
testbarcodes2 <- provided_reference@meta.data
Testlabels2 <- provided_reference@meta.data[["scpred_prediction"]]
table2 = data.frame()
for (i in 1:nrow(testbarcodes2)){
  table2<-rbind(table2,c(rownames(testbarcodes2[i,]),Testlabels2[i]))
}
colnames(table2)<- c("Barcodes","scpred_prediction")
plot(getSankey(table[,2],table2[,2]))



# lets find ARI
devtools::install_github("yycunc/SAMEclustering")
library("SAMEclustering")
dim(bm$scpred_prediction)
library(cidr)
# Cell labels of ground truth
head(bm$scpred_prediction)
head(bm$celltype.l2)
# Calculating ARI for cluster ensemble
adjustedRandIndex(bm$scpred_prediction, bm$celltype.l2)
adjustedRandIndex(bm$celltype.l2, bm$scpred_prediction)

provided_reference
adjustedRandIndex(provided_reference$cell_type, provided_reference$scpred_prediction)
adjustedRandIndex(provided_reference$celltype.l2, bm$scpred_prediction)

