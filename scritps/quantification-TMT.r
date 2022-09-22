#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 4) {
  cat("Usage: quantification input_dir output_dir tmt_plex RefTag\n")
  q(status = 1)
}
tmt_plex=as.numeric(argv[3])

sample_list=dir(argv[1])
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  peptide = read.table(paste(argv[1],"/",tmt_name,"/sepepe-with-class.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
  peptide=peptide[!is.na(peptide[,ncol(peptide)]),]
  expression_table=as.data.frame(peptide[,c((ncol(peptide)-tmt_plex+1):ncol(peptide))])  
  by_list=list(peptide$SEPEP_label)  
  expression_table=as.matrix(expression_table)
  expression_by_median=aggregate(x = expression_table, by = by_list, FUN = "median",na.rm=TRUE, na.action=NULL)
  rownames(expression_by_median)=expression_by_median$Group.1
  expression_by_median=expression_by_median[,-1]  
  group_size=as.data.frame(table(peptide$sorted_protein_list))
  rownames(group_size)=group_size$Var1  
  group_size=group_size[rownames(expression_by_median),]
  colnames(group_size)=c("group","peptides") 
  
  expression_by_median=expression_by_median[,!grepl(argv[4],colnames(expression_by_median))]
  write.table(rbind(idx=colnames(expression_by_median),expression_by_median),file = paste(argv[1],"/",tmt_name,"/sepep-level-quantification.txt",sep=""),row.names = TRUE,col.names = FALSE,sep="\t",quote = FALSE)
}

for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  peptide = read.table(paste(argv[1],"/",tmt_name,"/processed-protein-and-ratio.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE)
  high_quality_biclique = read.table(paste(argv[1],"/",tmt_name,"/filtered-sepep-by-fdr.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
  peptide=peptide[!is.na(peptide[,ncol(peptide)]),]        ###remove na
  peptide=peptide[peptide$`Mapped Genes`=="",] ###remove mapped to multiple genes
  peptide=peptide[peptide$sorted_protein_list %in% high_quality_biclique$sorted_protein_list,] ###remove peptide from low quality biclique
  #peptide = peptide[!is.na(peptide[,22]),]
  expression_table=as.data.frame(peptide[,c((ncol(peptide)-tmt_plex+1):ncol(peptide))])  
  by_list=list(peptide$Gene)  
  expression_table=as.matrix(expression_table)
  expression_by_median=aggregate(x = expression_table, by = by_list, FUN = "median",na.rm=TRUE, na.action=NULL)
  rownames(expression_by_median)=expression_by_median$Group.1
  expression_by_median=expression_by_median[,-1]  
  group_size=as.data.frame(table(peptide$sorted_protein_list))
  rownames(group_size)=group_size$Var1  
  group_size=group_size[rownames(expression_by_median),]
  colnames(group_size)=c("group","peptides") 

  expression_by_median=expression_by_median[,!grepl(argv[4],colnames(expression_by_median))]
  write.table(rbind(idx=colnames(expression_by_median),expression_by_median),file = paste(argv[1],"/",tmt_name,"/gene-level-quantification.txt",sep=""),row.names = TRUE,col.names = FALSE,sep="\t",quote = FALSE)
  
}

