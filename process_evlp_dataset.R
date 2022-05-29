#move from v3 to v4 is changing all writes and saves to use gdpath "root" functionality, and hasing outputs to be sure of reproduceability


source("gdpath.R")
#custom functions
source(gdpath("Masters/R Script Folder/lib/make_rnk.R"))
source(gdpath("Masters/R Script Folder/lib/GSEA_pipeline.R"))
# source(gdpath("Masters/R Script Folder/lib/win_slash.R"))
# source(gdpath("Masters/R Script Folder/lib/short_path.R"))
# source(gdpath("Masters/R Script Folder/lib/Ricardo_functions.R"))
# source(gdpath("Masters/R Script Folder/lib/script_library_tools.R"))
# source(gdpath("Masters/R Script Folder/lib/hash_functions.R"))

suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(limma))

# setwd("C:/Users/Aaron Wong/OneDrive - UHN/Masters/Project #3/github/published_code")
suppressPackageStartupMessages(library(optparse))
library(renv)
renv::restore()

option_list <- list(
    # make_option(c("-h", "--help"), action="store_true", default=FALSE, 
    #              help="Show this help message and exit"),

    make_option(c("-a", "--root_directory"), default="NA",
        help=""),
    make_option(c("-b", "--sample_info_evlp_path"), default="",
        help=""),
    make_option(c("-c", "--evlp_cel_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory
sample_info_path<-opt$sample_info_evlp_path
cel_directory<-opt$evlp_cel_directory


# root_directory<-opt$root_directory
# sample_info_transplant_path<-opt$sample_info_transplant_path
# transplant_cel_directory<-opt$transplant_cel_directory

# root_directory<-"C:/Users/Aaron Wong/Desktop/New folder/"
# sample_info<-"C:/Users/Aaron Wong/Desktop/New folder/data/evlp_sample.info.tab"
# cel_directory<-"C:/Users/Aaron Wong/Desktop/New folder/GSE127055_RAW/"


###################
#Batch 2 (EVLP Processing)
###################

#lets load in the batch 2 (EVLP) CEL files
batch_2_location<-cel_directory
#well have to sort files, since some are not CEL files
files_in_batch_2_folder<-list.files(gdpath(batch_2_location))
batch2_cel_files<-files_in_batch_2_folder[grep(".CEL",files_in_batch_2_folder)]
print(paste0(length(batch2_cel_files)," CEL files found for batch 2 (EVLP)"))
#lets load in the sample info from Ricardo so we can match the data tables
batch2_sample_info<-read.delim(sample_info_path,row.names=1)
#lets install the annotation and probe set database
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/clariomdhumanhsentrezg.db_22.0.0.tar.gz",repos=NULL,type="source")
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.clariomdhuman.hs.entrezg_22.0.0.tar.gz",repos=NULL,type="source")#mount our libraries
library(clariomdhumanhsentrezg.db)
library(pd.clariomdhuman.hs.entrezg)
#lets load in the batch 2 files, rma, summarize
raw.data <- read.celfiles(filenames = paste0(cel_directory,batch2_sample_info$GEO_filename), sampleNames = rownames(batch2_sample_info), pkgname = "pd.clariomdhuman.hs.entrezg")
rma.data <- rma(raw.data)
rma.expr <- exprs(rma.data)
#create annotation set
annotations<-data.frame(select(clariomdhumanhsentrezg.db,keys=keys(clariomdhumanhsentrezg.db),columns=c("ENTREZID","SYMBOL","GENENAME")))
#write.table(annotations,file=(gdpath("Masters/Project #3/Comparisons/PreRanked Batch0 vs Batch2 Human_GOBP_AllPathways_no_GO_iea_June_01_2017-common genes_pipeline_improved/batch2_annotations.tab")),sep="\t")
rownames(annotations)<-annotations[,1]
#unmount probe sets and annotations
remove.packages(c("clariomdhumanhsentrezg.db","pd.clariomdhuman.hs.entrezg"))
#combine annotations and gene expressionn matrix
batch2_eset<-rma.expr[rownames(rma.expr) %in% rownames(annotations),]
#save expression set
saveRDS(batch2_eset, file =paste0(root_directory,"/results/batch2_eset.rds"))
