#!/usr/bin/env Rscript
## This is the R script version of the R markdown file: predict_cell_type_using_pretrained_model.rmd
## Please modified the paths of the pretrained model and your query data
## Also set the file names of the output in the second code chunk below
## Then simply run this code by: ./predict_cell_types_using_pretrained_model.R

################################################################################################################################
# change logs:
#
# 04/24/2024: Xinlei Gao added the codes (two methods) to fix the Seurat version compitability issue
################################################################################################################################


## ------------------------------------------------------------------------------------------------
#set the R library path to the old one: 2023q1
set_lib_paths <- function(lib_vec) {
  lib_vec <- normalizePath(lib_vec, mustWork = TRUE)
  shim_fun <- .libPaths
  shim_env <- new.env(parent = environment(shim_fun))
  shim_env$.Library <- character()
  shim_env$.Library.site <- character()
  environment(shim_fun) <- shim_env
  shim_fun(lib_vec)
}

set_lib_paths(c("/nfs/apps/lib/R/4.2-focal/site-library.2023q1", "/opt/R/4.2.1/lib/R/library"))

# alternatively,
#library(irlba, lib.loc = "/nfs/apps/lib/R/4.2-focal/site-library.2023q4") # the 'irlba' package is required by Seurat for linear algebra
#library(Seurat, lib.loc = "/nfs/apps/lib/R/4.2-focal/site-library.2023q1")

## ------------------------------------------------------------------------------------------------
library(scPred)
library(Seurat)
library(cowplot)
library(ComplexHeatmap)


## ------------------------------------------------------------------------------------------------
# define your pretrained model, modified in your case
model_path <- "/nfs/BaRC_Public/BaRC_code/R/scRNAseq_cell_type_classification/reference_data/Tabula_sapiens/model/Tabula_sapiens_downsample100_scPred_default_model.rds"
# file containing the query data seurat object, here I use human lung sample as an example
query_path <- "/nfs/BaRC_Public/BaRC_code/R/scRNAseq_cell_type_classification/sample_dataset/human_lung_seurat_object_ENSEMBL_ID.rds"
# output file containing the prediction score matrix
output_matrix <- "./scPred_cell_type_prediction_score_matrix.rds"
# output UMAP plots with the original and predicted labels side-by-side
output_umap <- "./compare_human_lung_true_labels_with_scpred.pdf"
# output heatmap displaying confusion matrix
output_heatmap <- "./confusion_matrix_heatmap.pdf"
# output true labels and predicted labels
output_true_label <- "./human_lung_cell_type_true_labels.csv"
output_predicted_label <- "./human_lung_cell_type_scpred_prediction.csv"


## ------------------------------------------------------------------------------------------------
scpred <- readRDS(model_path)


## ------------------------------------------------------------------------------------------------
query <- readRDS(query_path)


## ------------------------------------------------------------------------------------------------
# Normalize the query data
query <- NormalizeData(query)
# Predict cell type labels using scPred
query <- scPredict(query, scpred) # default threshold = 0.55


## ------------------------------------------------------------------------------------------------
# Plot a histogram of the maximum probability
hist(query$scpred_max)


## ------------------------------------------------------------------------------------------------
#query <- scPredict(query, scpred, recompute_alignment = FALSE, threshold = 0.3)


## ------------------------------------------------------------------------------------------------
table(query$scpred_prediction)


## ------------------------------------------------------------------------------------------------
scpred_score_mtx <- query@meta.data[,grepl("scpred_", colnames(query@meta.data))]
# output the prediction score matrix
saveRDS(scpred_score_mtx, file = output_matrix)


## ------------------------------------------------------------------------------------------------
# Run UMAP using aligned data as input
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)

# Compare true labels with predicted labels
p1 <- DimPlot(query, group.by = "cell_type1", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(query, reduction = "umap", group.by = "scpred_prediction", label = TRUE, repel = TRUE) + NoLegend()

# Save the comparison plot
pdf(output_umap, 
    height = 8, width = 15)
options(ggrepel.max.overlaps = Inf)
plot_grid(p1, p2)
dev.off()


## ------------------------------------------------------------------------------------------------
confusion_matrix <- crossTab(query, "cell_type1", "scpred_prediction", output = "prop")
# draw heatmap
pdf(output_heatmap, height = 16, width = 16)
pheatmap(as.matrix(confusion_matrix), cluster_rows = TRUE, cluster_cols = TRUE, column_names_side = "top", fontsize_number = 5)
dev.off()


## ------------------------------------------------------------------------------------------------
write.csv(query$cell_type1, 
          file = output_true_label)
write.csv(query$scpred_prediction, 
          file = output_predicted_label)

