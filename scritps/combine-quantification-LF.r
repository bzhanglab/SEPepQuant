#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 2) {
  cat("Usage: combine-quantification-LF input_dir output_dir\n")
  q(status = 1)
}

sample_list=dir(argv[1])
gene_list="XXXXXXXX"
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  expression_table=read.table(paste(argv[1],"/",tmt_name,"/sepep-level-quantification.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE)
  gene_list=union(gene_list,rownames(expression_table))
}
gene_list=gene_list[-1]

combined_results=matrix(data=NA,ncol = 1,nrow=length(gene_list))
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  expression_table=read.table(paste(argv[1],"/",tmt_name,"/sepep-level-quantification.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE)
  expression_table[setdiff(gene_list,rownames(expression_table)),]=NA
  expression_table=expression_table[gene_list,]
  colnames(expression_table)=paste(tmt_name,colnames(expression_table),sep="_")
  combined_results=cbind(combined_results,expression_table)
  print(ncol(combined_results))
}
combined_results=combined_results[,-1]


write.table(rbind(idx=colnames(combined_results),combined_results),file=paste(argv[2],"/","sepep_quant_matrix.txt",sep=""),sep="\t",row.names=TRUE,quote=FALSE,col.names = FALSE)

###############################################################################################################################


