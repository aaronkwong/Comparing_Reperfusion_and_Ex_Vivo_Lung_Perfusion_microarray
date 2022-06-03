suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))

# source custom functino
source(paste0("custom_functions/GSEA_pipeline.R"))
source(paste0("custom_functions/make_rnk.R"))

option_list <- list(
    make_option(c("-a", "--root_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory

#################################
#Run GSEA for batch 0
#################################

s.broad_GSEA_preranked(
	s.gsea_soft=paste0(root_directory,"/gsea/gsea-3.0.jar"),
	s.gsea_memory=4096,
	s.gsea_nperm=1000,
	s.gsea_rnk=paste0(root_directory,"/results/batch_0_commongenes.rnk"),
	s.gsea_gmt=paste0(root_directory,"/gsea/pathways.gmt"), # here we renamed gmt file due to path lenght limitation of java
	s.gsea_output_name="batch0_enrichment",
	s.gsea_output_location=paste0(root_directory,"/results"),
	s.timestamp=1499973691546
)



#################################
#Run GSEA for batch 2
#################################

s.broad_GSEA_preranked(
	s.gsea_soft=paste0(root_directory,"/gsea/gsea-3.0.jar"),
	s.gsea_memory=4096,
	s.gsea_nperm=1000,
	s.gsea_rnk=paste0(root_directory,"/results/batch_2_commongenes.rnk"),
	s.gsea_gmt=paste0(root_directory,"/gsea/pathways.gmt"),# here we renamed gmt file due to path lenght limitation of java
	s.gsea_output_name="batch2_enrichment",
	s.gsea_output_location=paste0(root_directory,"/results"),
	s.timestamp=1234567890123
)
