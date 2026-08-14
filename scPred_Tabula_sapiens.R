#!/usr/bin/env Rscript
#
################################################################################################################################
# 04/18/2024: The default version of Seurat was changed to version 5. However, scPred is only compatible with Seurat version 4. 
#  Xinlei Gao changed the command line to load Seurat version 4. 
# 04/24/2024: Xinlei Gao added an alternative way to fix the Seurat version compitability issue
################################################################################################################################
# use scPred to train a model on Tabula Sapiens data
args <- commandArgs(T)
reference_path <- args[1]
output_path <- args[2]
model_filename <- args[3]
downsample_size <- args[4] # set this option to sample n cells from each cell type; or set as NULL to use full dataset.
n_cpus <- args[5]

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
# import packages 
library("scPred")
library("Seurat")
library("magrittr")

# an alternative way to load Seurat version 4
#library(irlba, lib.loc = "/nfs/apps/lib/R/4.2-focal/site-library.2023q4") # the 'irlba' package is required by Seurat for linear algebra
#library(Seurat, lib.loc = "/nfs/apps/lib/R/4.2-focal/site-library.2023q1")

# print session info
sessionInfo()

# You can run this script to train other models
reference <- readRDS(reference_path)
Idents(reference) <- reference$cell_type
# (optional but recommended) subset reference to reduce size
# Downsample the number of cells per identity class
if(!is.null(downsample_size) & downsample_size!="NULL")
  reference <- subset(reference, downsample = as.numeric(downsample_size))
# standard preprocessing steps
reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

#DimPlot(reference, group.by = "cell_type", label = FALSE) + NoLegend()

# Training classifiers with scPred
reference <- getFeatureSpace(reference, "cell_type")

# train the classifiers for each cell using the trainModel function. 
# By default, scPred will use a support vector machine with a radial kernel.
# parallel training
library(doParallel)
cl <- makePSOCKcluster(as.numeric(n_cpus))
registerDoParallel(cl)
reference <- trainModel(reference, allowParallel = TRUE)
stopCluster(cl)

# to visualize the performance
pdf(paste0(output_path, "/", "scPred_default_model_performance.pdf"))
plot_probabilities(reference)
dev.off()

# We can use the get_scpred method to retrieve the scPred object from the Seurat object
# save the model without seurat object
scpred <- get_scpred(reference)
saveRDS(scpred, file = paste0(output_path, "/", model_filename))

# using command line to train the model
#$./train_model.sh -r ./reference_data/Tabula_murius/SMART-seq/full_seurat_object.rds -o ./reference_data/Tabula_murius/SMART-seq/ -m Tabula_murius_all_default_model.rds -n 32
