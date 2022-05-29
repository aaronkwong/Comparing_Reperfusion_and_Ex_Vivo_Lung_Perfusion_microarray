#move from v3 to v4 is changing all writes and saves to use gdpath "root" functionality, and hasing outputs to be sure of reproduceability

suppressPackageStartupMessages(library(optparse))
option_list <- list(
    # make_option(c("-h", "--help"), action="store_true", default=FALSE, 
    #              help="Show this help message and exit"),

    make_option(c("-a", "--root_directory"), default="NA",
        help=""),
    make_option(c("-b", "--sample_info_transplant_path"), default="NA",
        help=""),
    make_option(c("-c", "--transplant_cel_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory
sample_info_transplant_path<-opt$sample_info_transplant_path
transplant_cel_directory<-opt$transplant_cel_directory

# root_directory<-"C:/Users/Aaron Wong/Desktop/New folder/"
# sample_info_transplant_path<-"C:/Users/Aaron Wong/Desktop/New folder/data/transplant_sample.info.tab"
# transplant_cel_directory<-"C:/Users/Aaron Wong/Desktop/New folder/GSE127003_RAW/"

#custom functions
# source(gdpath("Masters/R Script Folder/lib/make_rnk.R"))
# source(gdpath("Masters/R Script Folder/lib/Aaron_GSEA_pipeline_utility.R"))
# source(gdpath("Masters/R Script Folder/lib/win_slash.R"))
# source(gdpath("Masters/R Script Folder/lib/short_path.R"))
# source(gdpath("Masters/R Script Folder/lib/Ricardo_functions.R"))
# source(gdpath("Masters/R Script Folder/lib/script_library_tools.R"))
# source(gdpath("Masters/R Script Folder/lib/hash_functions.R"))


#libraries
# setwd("C:/Users/Aaron Wong/OneDrive - UHN/Masters/Project #3/github/published_code")
library(renv)
renv::activate()
library(oligo)
library(limma)

###################
#Batch 0 (LTx Processing)
###################

#we have a table with CEl filenames and sample id so we cant just grab from the table
batch0_sample_info<-read.delim(sample_info_transplant_path,row.names=1)
#lets load in the batch 0 (EVLP) CEL files
batch_0_location<-transplant_cel_directory
#well have to sort files, since some are not CEL files
files_in_batch_0_folder<-list.files(batch_0_location)
batch0_cel_files<-files_in_batch_0_folder[grep(".CEL",files_in_batch_0_folder)]
print(paste0(length(batch0_cel_files)," CEL files found for batch 0 (LTx)"))

#lets install the annotation and probe set database (these are binary packages)
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.hgu133plus2.hs.entrezg_22.0.0.tar.gz",repos=NULL,type="source")
install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133plus2hsentrezg.db_22.0.0.tar.gz",repos=NULL,type="source")
#mount our libraries
library(pd.hgu133plus2.hs.entrezg)
library(hgu133plus2hsentrezg.db)
#lets double check that the column of filename in the sample info table matches the names of the CEL files we actually have
sum(batch0_sample_info[,"GEO_filename"] %in% batch0_cel_files)==92
#lets load in the batch 2 files, rma, summarize. Well grab these from the table
raw.data <- read.celfiles(filenames = paste0(transplant_cel_directory,batch0_sample_info[,"GEO_filename"]), sampleNames = rownames(batch0_sample_info), pkgname = "pd.hgu133plus2.hs.entrezg")
rma.data <- rma(raw.data)
rma.expr <- exprs(rma.data)
#create annotation set
annotations<-data.frame(select(hgu133plus2hsentrezg.db,keys=keys(hgu133plus2hsentrezg.db),columns=c("ENTREZID","SYMBOL","GENENAME")))
#write.table(annotations,file=(gdpath("Masters/Project #3/Comparisons/PreRanked Batch0 vs Batch2 Human_GOBP_AllPathways_no_GO_iea_June_01_2017-common genes_pipeline_improved/batch0_annotations.tab")),sep="\t")
rownames(annotations)<-annotations[,1]
#unmount probe sets and annotations
remove.packages(c("hgu133plus2hsentrezg.db","pd.hgu133plus2.hs.entrezg"))
#combine annotations and gene expressionn matrix
batch0_eset<-rma.expr[rownames(rma.expr) %in% rownames(annotations),]
#save expression set
saveRDS(batch0_eset, file =paste0(root_directory,"results/batch0_eset.rds"))