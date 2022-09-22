#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 1) {
  cat("Usage: identify-sepep-LF input_dir \n")
  q(status = 1)
}

peptide_table = read.table(paste(argv[1],"/processed-protein-group.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE) 
by_list=list(peptide_table[,"sorted_protein_list"])  
expression_table=as.matrix(peptide_table[,c("Intensity","Spectral Count")])
expression_by_median=aggregate(x = expression_table, by = by_list, FUN = "sum",na.rm=TRUE, na.action=NULL)
rownames(expression_by_median)=expression_by_median$Group.1
expression_by_median=expression_by_median[,-1] 
colnames(expression_by_median)=c("SEPEP_Intensity","SEPEP_Spectral_Count")  
group_size=as.data.frame(table(peptide_table$sorted_protein_list))
rownames(group_size)=group_size$Var1  
group_size=group_size[rownames(expression_by_median),]
colnames(group_size)=c("group","peptides") 
  
peptide_probability=peptide_table # as.data.frame(peptide_table[,c(1:(ncol(peptide_table)-2))])
peptide_probability=peptide_probability[order(peptide_probability$Probability,peptide_probability$sorted_protein_list,decreasing = TRUE),]
peptide_probability=peptide_probability[!duplicated(peptide_probability$sorted_protein_list),]
rownames(peptide_probability)=peptide_probability$sorted_protein_list
peptide_probability=peptide_probability[rownames(expression_by_median),]
  
expression_by_median=cbind(peptide_probability,group_size,expression_by_median)
#expression_by_median=expression_by_median[,-(ncol(peptide_table)-tmt_plex+1)]
write.table(expression_by_median,file = paste(argv[1],"/identify-sepep-and-quantified-by-sum.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)

