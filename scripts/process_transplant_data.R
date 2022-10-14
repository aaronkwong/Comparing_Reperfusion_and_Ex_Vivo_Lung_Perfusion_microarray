# R package management
library(renv)
renv::activate()

suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-a", "--root_directory"), default="NA",
        help=""),
    make_option(c("-b", "--sample_info_transplant_path"), default="NA",
        help=""),
    make_option(c("-c", "--transplant_cel_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# root directory
root_directory<-opt$root_directory
# path to transplant sample info
sample_info_transplant_path<-opt$sample_info_transplant_path
# directory of transplant cel files
cel_directory<-opt$transplant_cel_directory

############
#Batch 0 (LTx Processing)
###################

# get CEl filenames from sample file
batch0_sample_info<-read.delim(sample_info_transplant_path,row.names=1)
# we will have to sort files, since some are not CEL files
files_in_batch_0_folder<-list.files(cel_directory)
batch0_cel_files<-files_in_batch_0_folder[grep(".CEL",files_in_batch_0_folder)]
print(paste0(length(batch0_cel_files)," CEL files found for batch 0 (LTx)"))

# lets install the annotation and probe set database (these are binary packages)
install.packages("./brain_array/hgu133plus2/hgu133plus2hsentrezg.db_22.0.0.tar.gz",repos=NULL,type="source")
install.packages("./brain_array/hgu133plus2/pd.hgu133plus2.hs.entrezg_22.0.0.tar.gz",repos=NULL,type="source")
# mount our libraries
suppressPackageStartupMessages(library(pd.hgu133plus2.hs.entrezg))
suppressPackageStartupMessages(library(hgu133plus2hsentrezg.db))
# lets double check that the column of filename in the sample info table matches the names of the CEL files we actually have
if(!sum(batch0_sample_info[,"GEO_filename"] %in% batch0_cel_files)==92){
	cat("Error. 92 Samples were expected but this number were not found.")
	quit()
}

# lets load in the batch 2 files, rma, summarize. Well grab these from the table
raw.data <- read.celfiles(filenames = paste0(cel_directory,"/",batch0_sample_info[,"GEO_filename"]), sampleNames = rownames(batch0_sample_info), pkgname = "pd.hgu133plus2.hs.entrezg")
# RMA normalization
rma.data <- rma(raw.data)
rma.expr <- exprs(rma.data)
#create annotation set
annotations<-data.frame(select(hgu133plus2hsentrezg.db,keys=keys(hgu133plus2hsentrezg.db),columns=c("ENTREZID","SYMBOL","GENENAME")))
# change rownames from probeset id to entrezid
rownames(annotations)<-annotations[,1]
#unmount probe sets and annotations
remove.packages(c("hgu133plus2hsentrezg.db","pd.hgu133plus2.hs.entrezg"))
#combine annotations and gene expressionn matrix
batch0_eset<-rma.expr[rownames(rma.expr) %in% rownames(annotations),]
#create results file
dir.create("results")
#save results
saveRDS(batch0_eset, file =paste0(root_directory,"/results/batch0_eset.rds"))
saveRDS(annotations, file =paste0(root_directory,"/results/batch0_annotations.rds"))