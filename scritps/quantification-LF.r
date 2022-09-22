#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 2) {
  cat("Usage: quantification-LF input_dir output_dir\n")
  q(status = 1)
}

sample_list=dir(argv[1])
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  peptide = read.table(paste(argv[1],"/",tmt_name,"/sepepe-with-class.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
  peptide=peptide[!is.na(peptide[,ncol(peptide)]),]
  expression_table=as.data.frame(peptide[,c("SEPEP_Intensity","SEPEP_Spectral_Count")])  
  by_list=list(peptide$SEPEP_label)  
  expression_table=as.matrix(expression_table)
  expression_by_median=aggregate(x = expression_table, by = by_list, FUN = "sum",na.rm=TRUE, na.action=NULL)
  rownames(expression_by_median)=expression_by_median$Group.1
  expression_by_median=expression_by_median[,-1]  
  group_size=as.data.frame(table(peptide$sorted_protein_list))
  rownames(group_size)=group_size$Var1  
  group_size=group_size[rownames(expression_by_median),]
  colnames(group_size)=c("group","peptides") 
  
  write.table(rbind(idx=colnames(expression_by_median),expression_by_median),file = paste(argv[1],"/",tmt_name,"/sepep-level-quantification.txt",sep=""),row.names = TRUE,col.names = FALSE,sep="\t",quote = FALSE)
}

