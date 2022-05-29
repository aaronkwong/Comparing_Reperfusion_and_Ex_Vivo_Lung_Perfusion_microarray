#!/bin/bash 

script_dir="/mnt/c/Users/Aaron Wong/OneDrive - UHN/Masters/Project #3/github/published_code"
root_directory="/mnt/c/Users/Aaron Wong/Desktop/New folder/"
sample_info_transplant_path="/mnt/c/Users/Aaron Wong/Desktop/New folder/data/transplant_sample.info.tab"
sample_info_evlp_path="/mnt/c/Users/Aaron Wong/Desktop/New folder/data/evlp_sample.info.tab"
transplant_cel_directory="/mnt/c/Users/Aaron Wong/Desktop/New folder/GSE127003_RAW/"
evlp_cel_directory="/mnt/c/Users/Aaron Wong/Desktop/New folder/GSE127055_RAW/"

# setup R virtual enviornmment for reproducibility
Rscript "${script_dir}/setup_renv.R"

# normalize transplant microarrays (pre and post) using RMA and output gene exprssion matrix
Rscript "${script_dir}/process_transplant_data.R" \
--root_directory="$root_directory" \
--sample_info_transplant_path="$sample_info_transplant_path" \
--transplant_cel_directory="$transplant_cel_directory"

# normalize EVLP microarrays (pre and post) using RMA and output gene exprssion matrix
Rscript "${script_dir}/process_evlp_dataset.R" \
--root_directory="$root_directory" \
--sample_info_evlp_path="$sample_info_evlp_path" \
--evlp_cel_directory="$evlp_cel_directory"

# conduct differential gene expression for post vs pre trandsplant
Rscript "${script_dir}/differential_gene_expression_transplant.R" \
--root_directory="$root_directory" \
--sample_info_transplant_path="$sample_info_transplant_path"

# conduct differential gene expression for post vs pre evlp
Rscript "${script_dir}/differential_gene_expression_evlp.R" \
--root_directory="$root_directory" \
--sample_info_evlp_path="$sample_info_evlp_path"

#take differentially expressed genes from transplant and EVLP, make ranked list for each using only intersecting genes
Rscript "${script_dir}/make_common_ranks.R" \
--root_directory="$root_directory"

# run gene set enrichment analysis on each dataset
Rscript "${script_dir}/GSEA.R" \
--root_directory="$root_directory"