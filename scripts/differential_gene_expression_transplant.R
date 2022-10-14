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
    make_option(c("-b", "--sample_info_transplant_path"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# root directory
root_directory<-opt$root_directory
# transplant sample info
sample_info_path<-opt$sample_info_transplant_path

#################################
#Limma Differential Gene Expression Batch_0
#################################

# load files
batch0_sample_info<-read.delim(sample_info_path,row.names=1)
batch0_eset<-readRDS(file =paste0(root_directory,"/results/batch0_eset.rds"))

# create design matrix
batch_0_patient<-factor(batch0_sample_info[,"sample.id"])
batch_0_treatment<-factor(batch0_sample_info[,"timepoint"],levels=c("cit","rep.2h"))
batch_0_design<-model.matrix(~batch_0_patient+batch_0_treatment)
batch0_fit<-lmFit(batch0_eset,batch_0_design)

# fit model and empirical Bayes 
batch0_fit<-eBayes(batch0_fit)
p.value<-batch0_fit$p.value[,"batch_0_treatmentrep.2h"]

# summary output
batch0_summ<-data.frame(
	p<-p.value,
	fdr=p.adjust(p.value, method="fdr")
	)
	
# Fold Change
batch0_log2.pair<-vector(mode="numeric", length=nrow(batch0_eset))
names(batch0_log2.pair)<- rownames(batch0_eset)
for (i in 1:nrow(batch0_eset)){
	x<-batch0_eset[i,(batch0_sample_info[,"timepoint"]=="cit")]
	y<-batch0_eset[i,(batch0_sample_info[,"timepoint"]=="rep.2h")]
	log2<-y-x
	batch0_log2.pair[i]<-mean(as.numeric(log2))
}

# create a table to output
batch0_summ$log2ratio.paired<-batch0_log2.pair
batch0_summ$fold.paired<-log2ratio2fold(batch0_log2.pair)
colnames(batch0_summ)<-c("p.value","fdr","log2ratio.paired","fold.paired")
rownames(batch0_summ)<-gsub("_at","",rownames(batch0_summ))

# make summary results
write.table(batch0_summ, file =paste0(root_directory,"/results/batch_0_improved_summ.tab"),row.names=TRUE,sep="\t",col.names=NA)

# make summary results with annotations
annotations<-readRDS(file=paste0(root_directory,"/results/batch0_annotations.rds"))
batch0_summ_with_annotations<-cbind(batch0_summ,annotations[match(rownames(batch0_summ),annotations$ENTREZID),])
write.table(batch0_summ_with_annotations, file =paste0(root_directory,"/results/batch_0_improved_summ_with_anno.tab"),row.names=TRUE,sep="\t",col.names=NA)


# make rnk file
make_rnk(x=batch0_summ,output_path=paste0(root_directory,"/results/batch0.rnk"),clean=T)

