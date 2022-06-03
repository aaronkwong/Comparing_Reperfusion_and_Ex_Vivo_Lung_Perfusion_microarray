suppressPackageStartupMessages(library(renv))
renv::activate()
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(limma))

source(paste0("custom_functions/GSEA_pipeline.R"))
source(paste0("custom_functions/make_rnk.R"))



option_list <- list(
    make_option(c("-a", "--root_directory"), default="NA",
        help=""),
    make_option(c("-b", "--sample_info_evlp_path"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory
sample_info_path<-opt$sample_info_evlp_path


#################################
#Limma Differential Gene Expression Batch_2
#################################

# load files
batch2_sample_info<-read.delim(sample_info_path,row.names=1)
batch2_eset<-readRDS(file =paste0(root_directory,"/results/batch2_eset.rds"))
sum(batch2_sample_info[,"filename"]==colnames(batch2_eset))

# create design matrix
batch_2_patient<-factor(batch2_sample_info[,"sample.id"])
batch_2_treatment<-factor(batch2_sample_info[,"timepoint"],levels=c("cit1","cit2"))
batch_2_design<-model.matrix(~batch_2_patient+batch_2_treatment)

# fit model and empirical Bayes 
batch2_fit<-lmFit(batch2_eset,batch_2_design)
batch2_fit<-eBayes(batch2_fit)
p.value<-batch2_fit$p.value[,"batch_2_treatmentcit2"]

# summary output
batch2_summ<-data.frame(
	p<-p.value,
	fdr=p.adjust(p.value, method="fdr")
	)
	
# Fold Change
batch2_log2.pair<-vector(mode="numeric", length=nrow(batch2_eset))
names(batch2_log2.pair)<- rownames(batch2_eset)
for (i in 1:nrow(batch2_eset)){
	x<-batch2_eset[i,(batch2_sample_info[,"timepoint"]=="cit1")]
	y<-batch2_eset[i,(batch2_sample_info[,"timepoint"]=="cit2")]
	log2<-y-x
	batch2_log2.pair[i]<-mean(as.numeric(log2))
}

# create a table to output
batch2_summ$log2ratio.paired<-batch2_log2.pair
batch2_summ$fold.paired<-log2ratio2fold(batch2_log2.pair)
colnames(batch2_summ)<-c("p.value","fdr","log2ratio.paired","fold.paired")
rownames(batch2_summ)<-gsub("_at","",rownames(batch2_summ))


# make summary results
write.table(batch2_summ, file =paste0(root_directory,"/results/batch_2_improved_summ.tab"),row.names=TRUE,sep="\t",col.names=NA)

# make rnk file
make_rnk(x=batch2_summ,output_path=paste0(root_directory,"/results/batch2.rnk"),clean=T)

