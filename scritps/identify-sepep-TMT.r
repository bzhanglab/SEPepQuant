#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 2) {
  cat("Usage: identify-sepep input_dir tmt_plex\n")
  q(status = 1)
}
tmt_plex=as.numeric(argv[2])

peptide_table = read.table(paste(argv[1],"/processed-protein-and-ratio.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE)
expression_table=as.data.frame(peptide_table[,c((ncol(peptide_table)-tmt_plex+1):ncol(peptide_table))])  
by_list=list(peptide_table[,"sorted_protein_list"])  
expression_table=as.matrix(expression_table)
expression_by_median=aggregate(x = expression_table, by = by_list, FUN = "median",na.rm=TRUE, na.action=NULL)
rownames(expression_by_median)=expression_by_median$Group.1
expression_by_median=expression_by_median[,-1]  
group_size=as.data.frame(table(peptide_table$sorted_protein_list))
rownames(group_size)=group_size$Var1  
group_size=group_size[rownames(expression_by_median),]
colnames(group_size)=c("group","peptides") 
  
peptide_probability=as.data.frame(peptide_table[,c(1:(ncol(peptide_table)-tmt_plex))])
peptide_probability=peptide_probability[order(peptide_probability$Probability,peptide_probability$sorted_protein_list,decreasing = TRUE),]
peptide_probability=peptide_probability[!duplicated(peptide_probability$sorted_protein_list),]
rownames(peptide_probability)=peptide_probability$sorted_protein_list
peptide_probability=peptide_probability[rownames(expression_by_median),]
  
expression_by_median=cbind(peptide_probability,group_size,expression_by_median)
expression_by_median=expression_by_median[,-(ncol(peptide_table)-tmt_plex+1)]
write.table(expression_by_median,file = paste(argv[1],"/identify-sepep-and-quantified-by-median.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)

