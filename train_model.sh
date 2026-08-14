#!/bin/bash

# Default values
default_memory="64G"
default_cpus=16
default_downsample_size="NULL"

# Define usage and parameter descriptions
function show_usage() {
    echo "Description: this script is used to train a model on well annotated scRNA-seq reference data (such as Tabula sapiens) to predict cell types on new query data using the R package scPred."
    echo "Usage: $0 -r <reference_path> -o <output_path> -m <model_filename> [-d <downsample_size>] [-n <n_cpus>] [-mem <memory_size>]"
    echo "Options:"
    echo "  -h, --help                Show this help message"
    echo "  -r, --reference_path      Path to the reference data (Seurat object)"
    echo "  -o, --output_path         Output path for saving the model and results"
    echo "  -m, --model_filename      Filename for saving the scPred model"
    echo "  -d, --downsample_size     Downsample size (default: $default_downsample_size to use the full dataset)"
    echo "  -n, --n_cpus              Number of CPUs for parallel training (default: $default_cpus)"
    echo "  -mem, --memory_size       Required memory size (default: $default_memory)"
}

# Set default values
downsample_size="$default_downsample_size"
n_cpus="$default_cpus"
memory_size="$default_memory"

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -r|--reference_path)
            reference_path="$2"
            shift
            ;;
        -o|--output_path)
            output_path="$2"
            shift
            ;;
        -m|--model_filename)
            model_filename="$2"
            shift
            ;;
        -d|--downsample_size)
            downsample_size="$2"
            shift
            ;;
        -n|--n_cpus)
            n_cpus="$2"
            shift
            ;;
        -mem|--memory_size)
            memory_size="$2"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
    shift
done

# Check if required parameters are provided
if [ -z "$reference_path" ] || [ -z "$output_path" ] || [ -z "$model_filename" ]; then
    echo "Error: Missing required parameters."
    show_usage
    exit 1
fi

# check if the output path exists, otherwise create it
if [ ! -d "$output_path" ]; then
  mkdir -p $output_path
fi

# Run Rscript using sbatch
sbatch -p 20 --job-name=scPred --mem="$memory_size" --cpus-per-task="$n_cpus" --wrap "/nfs/BaRC_Public/BaRC_code/R/scRNAseq_cell_type_classification/scPred_Tabula_sapiens.R "$reference_path" "$output_path" "$model_filename" "$downsample_size" "$n_cpus""

echo "Job submitted successfully!"
