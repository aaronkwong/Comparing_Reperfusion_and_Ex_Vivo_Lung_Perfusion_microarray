source(paste0("custom_functions/GSEA_pipeline.R"))
source(paste0("custom_functions/make_rnk.R"))


suppressPackageStartupMessages(library(optparse))
library(renv)
renv::activate()
option_list <- list(
    # make_option(c("-h", "--help"), action="store_true", default=FALSE, 
    #              help="Show this help message and exit"),

    make_option(c("-a", "--root_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory

#################################
#Run GSEA for batch 0
#################################

s.broad_GSEA_preranked(
	s.gsea_soft=paste0(root_directory,"/data/gsea-3.0.jar"),
	s.gsea_memory=4096,
	s.gsea_nperm=1000,
	s.gsea_rnk=paste0(root_directory,"/results/batch_0_commongenes.rnk"),
	s.gsea_gmt=paste0(root_directory,"/data/Human_GOBP_AllPathways_no_GO_iea_June_01_2017_entrezgene.gmt"),
	s.gsea_output_name="batch0_enrichment",
	s.gsea_output_location=paste0(root_directory,"/results"),
	s.timestamp=1499973691546
)



#################################
#Run GSEA for batch 2
#################################

s.broad_GSEA_preranked(
	s.gsea_soft=paste0(root_directory,"/data/gsea-3.0.jar"),
	s.gsea_memory=4096,
	s.gsea_nperm=1000,
	s.gsea_rnk=paste0(root_directory,"/results/batch_2_commongenes.rnk"),
	s.gsea_gmt=paste0(root_directory,"/data/Human_GOBP_AllPathways_no_GO_iea_June_01_2017_entrezgene.gmt"),
	s.gsea_output_name="batch2_enrichment",
	s.gsea_output_location=paste0(root_directory,"/results"),
	s.timestamp=1234567890123
)
