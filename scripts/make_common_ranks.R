suppressPackageStartupMessages(library(optparse))
renv::activate()
suppressPackageStartupMessages(library(renv))

# source custom function
source(paste0("custom_functions/GSEA_pipeline.R"))
source(paste0("custom_functions/make_rnk.R"))

option_list <- list(
    make_option(c("-a", "--root_directory"), default="NA",
        help="")
    )

opt <- parse_args(OptionParser(option_list=option_list))

root_directory<-opt$root_directory

#lets cutoff the genes so that only common ones between the platforms are used
batch_0<-read.table(paste0(root_directory,"/results/batch0.rnk"),sep="\t")
batch_2<-read.table(paste0(root_directory,"/results/batch2.rnk"),sep="\t")
common_genes<-intersect(batch_0[,1],batch_2[,1])

# pull out only intersecting genes
batch_0_common<-batch_0[batch_0[,1] %in% common_genes,]
batch_2_common<-batch_2[batch_2[,1] %in% common_genes,]

# remame columns
colnames(batch_0_common)=c("entrez.id","formula")
colnames(batch_2_common)=c("entrez.id","formula")

#write new ranked files with the interesecting genes
write.table(batch_0_common,file=paste0(root_directory,"/results/batch_0_commongenes.rnk"),row.names=F,sep="\t",col.names=T,quote=F)
write.table(batch_2_common,file=paste0(root_directory,"/results/batch_2_commongenes.rnk"),row.names=F,sep="\t",col.names=T,quote=F)