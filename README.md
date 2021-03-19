# STEMREM205-Course-Project-QB-JB

The choice of methodology used to determine distinct cell clusters within a single cell RNA-seq dataset is an extremely critical part of the analytical process and could vastly alter any downstream conclusions or derived correlations. In this study, we address the growing need to compare and evaluate different algorithms used to classify undetermined cells within a query dataset based on a matching protein-expression-based ground truth.

#scPred 

scPred_script.R can be utilized to allow Alquicira-Hernandez et al.'s scPred training model to recluster a provided PBMC dataset using a reference of Satija et al.'s CITEseq dbone marrow mononuclear cell data.

#Seurat and scMapCluster

Scmap is a method developed for projecting cells from an scRNA-seq data set onto cell types or individual cells for other experiments (Kiselev, Yiu and Hemberg, 2018).


Source code_Final Project_QB.R can be utilized to recluster Satija et al.'s CITEseq dbone marrow mononuclear cell data using Seurat and scMapCluster

#Visualization

Sankey_ARI_Generator.R allows for the generation of river plots that can compare the inclusion or exclusion of cells in different clustering methods


#Data sets

Data sets used in this study are available and loadable onto local devices once the aforementioned have been installed (specific instances cited in code).
