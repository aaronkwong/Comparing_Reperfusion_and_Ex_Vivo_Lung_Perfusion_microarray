#!/bin/bash 

#set the path to appropriate R version
R_version="/opt/R/4.0.3/bin/Rscript"

# analysis is meant to be run in the unzipped repository
root_directory="."
script_dir="./scripts"
sample_info_transplant_path="./sample_info/transplant_sample.info.tab"
sample_info_evlp_path="./sample_info/evlp_sample.info.tab"
transplant_cel_directory="./GSE127003_RAW"
evlp_cel_directory="./GSE127055_RAW"

# setup R virtual enviornmment for reproducibility
${R_version} "${script_dir}/setup_renv.R"

# normalize transplant microarrays (pre and post) using RMA and output gene exprssion matrix
${R_version} "${script_dir}/process_transplant_data.R" \
--root_directory="$root_directory" \
--sample_info_transplant_path="$sample_info_transplant_path" \
--transplant_cel_directory="$transplant_cel_directory"

# normalize EVLP microarrays (pre and post) using RMA and output gene exprssion matrix
${R_version} "${script_dir}/process_evlp_dataset.R" \
--root_directory="$root_directory" \
--sample_info_evlp_path="$sample_info_evlp_path" \
--evlp_cel_directory="$evlp_cel_directory"

# conduct differential gene expression for post vs pre trandsplant
${R_version} "${script_dir}/differential_gene_expression_transplant.R" \
--root_directory="$root_directory" \
--sample_info_transplant_path="$sample_info_transplant_path"

# conduct differential gene expression for post vs pre evlp
${R_version} "${script_dir}/differential_gene_expression_evlp.R" \
--root_directory="$root_directory" \
--sample_info_evlp_path="$sample_info_evlp_path"

#take differentially expressed genes from transplant and EVLP, make ranked list for each using only intersecting genes
${R_version} "${script_dir}/make_common_ranks.R" \
--root_directory="$root_directory"

# run gene set enrichment analysis on each dataset
${R_version} "${script_dir}/run_GSEA.R" \
--root_directory="$root_directory"
