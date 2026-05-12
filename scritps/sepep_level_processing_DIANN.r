#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 2) {
  cat("Usage:sepep_level_processing_DIANN input_dir output_dir\n")
  q(status = 1)
}

peptide_table = read.table(paste(argv[1],"/processed-report.pr_matrix.tsv",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE,na.strings = "NA")
expression_table=as.data.frame(peptide_table[,c(17:ncol(peptide_table))]) 
inf_table=as.data.frame(peptide_table[,c(1:16)])  
by_list=list(peptide_table$sorted_protein_list) ####for sepep quantification 
by_list_gene=list(peptide_table$Genes_refined)  ####for gene quantification 
expression_table = sapply(expression_table, as.numeric)
expression_by_sum=aggregate(x = expression_table, by = by_list, FUN = "sum",na.rm=TRUE, na.action=NULL)
expression_by_sum_gene=aggregate(x = expression_table, by = by_list_gene, FUN = "sum",na.rm=TRUE, na.action=NULL)
rownames(expression_by_sum)=expression_by_sum$Group.1
rownames(expression_by_sum_gene)=expression_by_sum_gene$Group.1
expression_by_sum=expression_by_sum[,-1]  
expression_by_sum_gene=expression_by_sum_gene[,-1] 
expression_by_sum_gene=expression_by_sum_gene[!grepl(",",rownames(expression_by_sum_gene)),]
group_size=as.data.frame(table(peptide_table$sorted_protein_list))
rownames(group_size)=group_size$Var1  
group_size=group_size[rownames(expression_by_sum),]
colnames(group_size)=c("group","number_of_peptides") 

inf_table=inf_table[!duplicated(inf_table$sorted_protein_list),]
rownames(inf_table)=inf_table$sorted_protein_list
inf_table=inf_table[rownames(expression_by_sum),]
combined_table=cbind(inf_table,group_size,expression_by_sum)

################sepep identification

gene_clique_number=as.data.frame(table(combined_table$sorted_gene_list))
rownames(gene_clique_number)=gene_clique_number$Var1
gene_clique_number$Freq=2

combined_table=combined_table[!duplicated(combined_table$sorted_protein_list),]

combined_table_multiple_gene=combined_table[combined_table$number_of_mapped_gene>1,]
combined_table_whole_gene=combined_table[combined_table$number_of_mapped_gene==1&(combined_table$number_of_mapped_protein==combined_table$protein_number_of_gene),]
combined_table_single_protein=combined_table[combined_table$number_of_mapped_gene==1&combined_table$number_of_mapped_protein==1&combined_table$protein_number_of_gene>1,]
combined_table_multiple_protein=combined_table[combined_table$number_of_mapped_gene==1&combined_table$number_of_mapped_protein>1&(combined_table$number_of_mapped_protein<combined_table$protein_number_of_gene),]

#setdiff(combined_table$sorted_protein_list,c(combined_table_multiple_gene$sorted_protein_list,combined_table_whole_gene$sorted_protein_list,combined_table_single_protein$sorted_protein_list,combined_table_multiple_protein$sorted_protein_list))
combined_table_multiple_gene[,"SEPEP_label"]="Multiple_SEPEP"
for(i in 1:nrow(combined_table_multiple_gene))
{
  combined_table_multiple_gene$SEPEP_label[i]=paste(combined_table_multiple_gene$SEPEP_label[i],".",i,"_C5",sep="")
  #print(combined_table_multiple_gene$SEPEP_label[i])
}

combined_table_whole_gene[,"SEPEP_label"]=combined_table_whole_gene$Genes_refined
combined_table_whole_gene=combined_table_whole_gene[!is.na(combined_table_whole_gene$number_of_mapped_protein),]
for(i in 1:nrow(combined_table_whole_gene))
{
  #print(combined_table_whole_gene$number_of_mapped_protein[i])
  if(combined_table_whole_gene$number_of_mapped_protein[i]==1)
  {
    combined_table_whole_gene[i,"SEPEP_label"]=paste(combined_table_whole_gene$Genes_refined[i],"_SEPEP.1_C1",sep="")
  } else {
    combined_table_whole_gene[i,"SEPEP_label"]=paste(combined_table_whole_gene$Genes_refined[i],"_SEPEP.1_C4",sep="")
  }
}

combined_table_single_protein[,"SEPEP_label"]=combined_table_single_protein$Genes_refined
combined_table_single_protein=combined_table_single_protein[!is.na(combined_table_single_protein$number_of_mapped_protein),]
#combined_table_single_protein$SEPEP_label=make.unique(combined_table_single_protein$SEPEP_label)
for(i in 1:nrow(combined_table_single_protein))
{
  combined_table_single_protein[i,"SEPEP_label"]=paste(combined_table_single_protein[i,"SEPEP_label"],"_SEPEP.",gene_clique_number[combined_table_single_protein$Genes_refined[i],"Freq"],"_C2",sep="")
  gene_clique_number[combined_table_single_protein$Genes_refined[i],"Freq"]=gene_clique_number[combined_table_single_protein$Genes_refined[i],"Freq"]+1
}

combined_table_multiple_protein[,"SEPEP_label"]=combined_table_multiple_protein$Genes_refined
combined_table_multiple_protein=combined_table_multiple_protein[!is.na(combined_table_multiple_protein$number_of_mapped_protein),]
#combined_table_multiple_protein$SEPEP_label=make.unique(combined_table_multiple_protein$SEPEP_label)
for(i in 1:nrow(combined_table_multiple_protein))
{
  combined_table_multiple_protein[i,"SEPEP_label"]=paste(combined_table_multiple_protein[i,"SEPEP_label"],"_SEPEP.",gene_clique_number[combined_table_multiple_protein$Genes_refined[i],"Freq"],"_C3",sep="")
  gene_clique_number[combined_table_multiple_protein$Genes_refined[i],"Freq"]=gene_clique_number[combined_table_multiple_protein$Genes_refined[i],"Freq"]+1
}

combined_table_with_sepep=rbind(combined_table_multiple_gene,combined_table_whole_gene,combined_table_single_protein,combined_table_multiple_protein)
mapping_table=combined_table_with_sepep[,c(ncol(combined_table_with_sepep),12:16,18)]

rownames(combined_table_with_sepep)=combined_table_with_sepep$SEPEP_label
sepep_expression_table=combined_table_with_sepep[,-ncol(combined_table_with_sepep)]
sepep_expression_table=sepep_expression_table[,c(-1:-18)]
sepep_expression_table=log2(sepep_expression_table)
sepep_expression_table[sepep_expression_table=="-Inf"]=NA

expression_by_sum_gene=log2(expression_by_sum_gene)
expression_by_sum_gene[expression_by_sum_gene=="-Inf"]=NA

write.table(rbind(idx=colnames(mapping_table),mapping_table),file=paste(argv[2],"/sepep_mapping_table.txt",sep=""),col.names = FALSE,row.names = FALSE,sep="\t",quote = FALSE)
write.table(rbind(idx=colnames(sepep_expression_table),sepep_expression_table),file=paste(argv[2],"/sepep_matrix_raw_log2.txt",sep=""),col.names = FALSE,row.names = TRUE,sep="\t",quote = FALSE)
write.table(rbind(idx=colnames(expression_by_sum_gene),expression_by_sum_gene),file=paste(argv[2],"/gene_matrix_raw_log2.txt",sep=""),col.names = FALSE,row.names = TRUE,sep="\t",quote = FALSE)


sample_median_gene_level=apply(expression_by_sum_gene,2,median,na.rm=TRUE)
sample_median_diff=sample_median_gene_level-median(sample_median_gene_level)

sepep_expression_table=t(t(sepep_expression_table)-sample_median_diff)
expression_by_sum_gene=t(t(expression_by_sum_gene)-sample_median_diff)
write.table(rbind(idx=colnames(sepep_expression_table),sepep_expression_table),file=paste(argv[2],"/sepep_matrix_MD_log2.txt",sep=""),col.names = FALSE,row.names = TRUE,sep="\t",quote = FALSE)
write.table(rbind(idx=colnames(expression_by_sum_gene),expression_by_sum_gene),file=paste(argv[2],"/gene_matrix_MD_log2.txt",sep=""),col.names = FALSE,row.names = TRUE,sep="\t",quote = FALSE)


