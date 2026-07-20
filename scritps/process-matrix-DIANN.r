#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 3) {
  cat("Usage: process-matrix-DIANN DIA-NN_output_folder gene-protein-statistic-table protein-gene-mapping-table\n")
  q(status = 1)
}

lines <- readLines(paste(argv[1],"/report.pr_matrix.tsv",sep=""), warn = FALSE)
lines <- gsub('"', "", lines)
lines <- gsub("\'", "", lines)
writeLines(lines, paste(argv[1],"/report.pr_matrix.tsv",sep=""))

peptide_table = read.table(paste(argv[1],"/report.pr_matrix.tsv",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
peptide_table=peptide_table[!grepl("^Cont",peptide_table$Protein.Ids,perl=TRUE),] ##remove decoy
gene_protein_number=read.table(argv[2],sep="\t",row.names = 1,header = FALSE)
protein_gene_mapping=read.table(argv[3],sep="\t",header = FALSE)
rownames(protein_gene_mapping)=protein_gene_mapping$V1

peptide_table$Genes_refined=peptide_table$Genes
for(j in 1:nrow(peptide_table))
{
  if(grepl(";",peptide_table$Protein.Ids[j]))
  {
    protein_list=as.vector(unlist(strsplit(peptide_table$Protein.Ids[j],"\\;")))
    protein_list=protein_list[!(grepl("^rev",protein_list,perl = TRUE)|grepl("^Cont",protein_list,perl = TRUE))] ##remove reversion and cont
    gene_list=protein_gene_mapping[protein_list,]$V2
    peptide_table$Genes_refined[j]=paste(unique(gene_list),sep=",",collapse =",")
  }
}
peptide_table=peptide_table[,c(1:4,ncol(peptide_table),5:(ncol(peptide_table)-1))]


peptide_inf=peptide_table[,c(1:11)]   ###peptide information
peptide_intensity=peptide_table[,c(12:ncol(peptide_table))]  ###peptide intensity

peptide_inf[,"sorted_protein_list"]=NA  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"sorted_gene_list"]=NA  ####remove gene name. Genes from decoy are excluded
peptide_inf[,"number_of_mapped_protein"]=0  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"number_of_mapped_gene"]=0  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"protein_number_of_gene"]=0  ####remove gene name. Remove decoy if has both target ad decoy

for(j in 1:nrow(peptide_inf))
{
  protein_list=peptide_inf$Protein.Ids[j]
  protein_list=gsub(" ","",protein_list)
  protein_list=as.vector(unlist(strsplit(protein_list,"\\;")))
  protein_list=protein_list[!(grepl("^rev",protein_list,perl = TRUE)|grepl("^Cont",protein_list,perl = TRUE))] ##remove reversion and cont
  peptide_inf[j,"sorted_protein_list"]=paste(sort(protein_list),sep=",",collapse =",")
  total_number_of_proteins=length(protein_list)
  peptide_inf[j,"number_of_mapped_protein"]=total_number_of_proteins
  #print(total_number_of_proteins)
  gene_list=peptide_inf$Genes_refined[j]
  gene_list=sub(" ","",gene_list)
  gene_list=as.vector(unlist(strsplit(gene_list,"\\,")))
  peptide_inf[j,"sorted_gene_list"]=paste(sort(gene_list),sep=",",collapse =",")
  peptide_inf[j,"number_of_mapped_gene"]=length(gene_list)
  if(length(gene_list)==1)
  {
    gene_name=gene_list[1]
    peptide_inf[j,"protein_number_of_gene"]=gene_protein_number[gene_name,1]
  }
}

peptide_intensity=format(peptide_intensity,scientific = TRUE)
combined_results=cbind(peptide_inf,peptide_intensity)
write.table(combined_results,file=paste(argv[1],"/processed-report.pr_matrix.tsv",sep=""),sep="\t",row.names=FALSE,col.names = TRUE,quote=FALSE)



