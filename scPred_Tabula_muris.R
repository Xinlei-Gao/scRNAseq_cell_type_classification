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


# 1. To train a model on full data of Tabula Muris (SMART-seq)
setwd("/nfs/BaRC_Public/BaRC_code/R/scRNAseq_cell_type_classification/")
count_files <- list.files("./reference_data/Tabula_muris/SMART-seq/FACS/")

reference_path <- "./reference_data/Tabula_muris/SMART-seq/FACS/"

# preprocess the SMART-seq data
# load the count matrix
count_files
tissue_names <- unlist(strsplit(count_files, split = "-counts.csv"))
count_matrix_list <- list()
for (i in 1:length(tissue_names)){
  count_matrix_list[[tissue_names[i]]] <- read.csv(paste0(reference_path, count_files[i]), row.names = 1)
}
# create seurat objects
seurat_objs <- list()
for (i in 1:length(tissue_names)){
  seurat_objs[[tissue_names[i]]] <- CreateSeuratObject(counts = count_matrix_list[[tissue_names[i]]],
                                                      project = tissue_names[i])
}

# save seurat object for each tissue separately
for (i in 1:length(tissue_names)){
  saveRDS(seurat_objs[[tissue_names[i]]], file = paste0("./reference_data/Tabula_muris/SMART-seq/", tissue_names[i], "_seurat_object.rds"))
}

# merge Seurat Object
seu_all <- merge(seurat_objs[["Aorta"]], y = c(seurat_objs[["Bladder"]], 
                                               seurat_objs[["Brain_Myeloid"]],
                                               seurat_objs[["Brain_Non-Myeloid"]],
                                               seurat_objs[["Diaphragm"]],
                                               seurat_objs[["Fat"]],
                                               seurat_objs[["Heart"]],
                                               seurat_objs[["Kidney"]],
                                               seurat_objs[["Large_Intestine"]],
                                               seurat_objs[["Limb_Muscle"]],
                                               seurat_objs[["Liver"]],
                                               seurat_objs[["Lung"]],
                                               seurat_objs[["Mammary_Gland"]],
                                               seurat_objs[["Marrow"]],
                                               seurat_objs[["Pancreas"]],
                                               seurat_objs[["Skin"]],
                                               seurat_objs[["Spleen"]],
                                               seurat_objs[["Thymus"]],
                                               seurat_objs[["Tongue"]],
                                               seurat_objs[["Trachea"]]
                                            ), project = "Tabula_murius")
# save the full seurat object
saveRDS(seu_all, file = "./reference_data/Tabula_muris/SMART-seq/full_seurat_object.rds")
# Tabula murius applied a QC cutoff of at least 500 genes and 50,000 reads
# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687
seu_all <- subset(seu_all, subset = nCount_RNA >= 50000 & nFeature_RNA >= 500)
# read metadata file
meta_data <- read.csv("./reference_data/Tabula_muris/SMART-seq/metadata_FACS.csv")
# load cell annotation 
cell_anno <- read.csv("./reference_data/Tabula_muris/SMART-seq/annotations_facs.csv")
# only keep the cells have annotation
cell_anno <- cell_anno[cell_anno$cell_ontology_class!="", ]
seu_all$cell.name <- rownames(seu_all@meta.data)
seu_all <- subset(seu_all, subset = cell.name %in% cell_anno$cell)
# sort the cell annotation table by the order of cell.name in seurat object metadata
cell_anno <- cell_anno[order(match(cell_anno$cell, seu_all$cell.name)),]
identical(cell_anno$cell, as.vector(seu_all$cell.name))
seu_all$cell_type <- cell_anno$cell_ontology_class
# save the seurat object with cell type annotation
saveRDS(seu_all, file = "./reference_data/Tabula_muris/SMART-seq/full_seurat_object.rds")

seu_all <- readRDS( "./reference_data/Tabula_muris/SMART-seq/full_seurat_object.rds")

Idents(seu_all) <- seu_all$cell_type

# standard preprocessing steps
seu_all <- seu_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(seu_all, group.by = "cell_type", label = FALSE) + NoLegend()

# Training classifiers with scPred
seu_all <- getFeatureSpace(seu_all, "cell_type")
saveRDS(seu_all, file = "./reference_data/Tabula_muris/SMART-seq/full_seurat_object.rds")

########################################################################
# the following lines can be run on slurm cluster
#./train_model.sh -r ./reference_data/Tabula_murius/SMART-seq/full_seurat_object.rds -o ./reference_data/Tabula_murius/SMART-seq/ -m Tabula_murius_all_default_model.rds -n 32 -mem 238G
########################################################################

# 2. train separate models for each tissue
setwd("/nfs/BaRC_Public/BaRC_code/R/scRNAseq_cell_type_classification/")

seu_all <- readRDS("./reference_data/Tabula_muris/SMART-seq/full_seurat_object.rds")

tissue_names <- unique(seu_all$orig.ident)
write.table(tissue_names, file = "./reference_data/Tabula_muris/SMART-seq/Tabula_Muris_SMART-seq_tissue_list.txt",
            quote=F, row.names = F, col.names = F)

for (i in 1:length(tissue_names)){
  seu <- subset(seu_all, subset = orig.ident == tissue_names[i])
  # save the seurat object of this tissue
  saveRDS(seu, file = paste0("./reference_data/Tabula_muris/SMART-seq/", tissue_names[i], "_seurat_object.rds"))
}

########################################################################
# the following lines can be run on slurm cluster
#while IFS=$'\t' read -r tissue_name
#do
#./train_model.sh -r ./reference_data/Tabula_muris/SMART-seq/${tissue_name}_seurat_object.rds -o ./reference_data/Tabula_muris/SMART-seq/ -m Tabula_murius_${tissue_name}_default_model.rds -n 32 -mem 238G
#done < ./reference_data/Tabula_muris/SMART-seq/Tabula_Muris_SMART-seq_tissue_list.txt
########################################################################
