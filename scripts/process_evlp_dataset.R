# R package management
library(renv)
renv::activate()

suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-a", "--root_directory"), default="NA",
        help=""),
    make_option(c("-b", "--sample_info_evlp_path"), default="",
        help=""),
    make_option(c("-c", "--evlp_cel_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# the path to the root directory
root_directory<-opt$root_directory
# path to evlp sample info
sample_info_path<-opt$sample_info_evlp_path
# directory of evlp CEL files
cel_directory<-opt$evlp_cel_directory

###################
#Batch 2 (EVLP Processing)
###################

# pull CEL filenames from our sample info table
batch2_sample_info<-read.delim(sample_info_path,row.names=1)
# we will have to sort files, since some are not CEL files
files_in_batch_2_folder<-list.files(cel_directory)
batch2_cel_files<-files_in_batch_2_folder[grep(".CEL",files_in_batch_2_folder)]
print(paste0(length(batch2_cel_files)," CEL files found for batch 2 (EVLP)"))
# lets install the annotation and probe set database
install.packages("./brain_array/clariomd/clariomdhumanhsentrezg.db_22.0.0.tar.gz",repos=NULL,type="source")
install.packages("./brain_array/clariomd/pd.clariomdhuman.hs.entrezg_22.0.0.tar.gz",repos=NULL,type="source")#mount our libraries
library(clariomdhumanhsentrezg.db)
library(pd.clariomdhuman.hs.entrezg)
#lets load in the batch 2 files, rma, summarize
raw.data <- read.celfiles(filenames = paste0(cel_directory,"/",batch2_sample_info$GEO_filename), sampleNames = rownames(batch2_sample_info), pkgname = "pd.clariomdhuman.hs.entrezg")
# RMA
rma.data <- rma(raw.data)
# extract expression
rma.expr <- exprs(rma.data)
# create annotation set
annotations<-data.frame(select(clariomdhumanhsentrezg.db,keys=keys(clariomdhumanhsentrezg.db),columns=c("ENTREZID","SYMBOL","GENENAME")))
# change rownames from probeset id to entrezid
rownames(annotations)<-annotations[,1]
#unmount probe sets and annotations
remove.packages(c("clariomdhumanhsentrezg.db","pd.clariomdhuman.hs.entrezg"))
#combine annotations and gene expressionn matrix
batch2_eset<-rma.expr[rownames(rma.expr) %in% rownames(annotations),]
# create results file
dir.create("results")
saveRDS(batch2_eset, file =paste0(root_directory,"/results/batch2_eset.rds"))
saveRDS(annotations, file =paste0(root_directory,"/results/batch2_annotations.rds"))
