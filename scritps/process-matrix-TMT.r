#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 4) {
  cat("Usage: process-FragPipe-output fragpipe_output_folder_for_one_sample tmt_plex Ref_Tag gene-protein-statistic-table\n")
  q(status = 1)
}

tmt_plex=as.numeric(argv[2])
gene_protein_number=read.table(argv[4],sep="\t",row.names = 1,header = TRUE)
peptide_table = read.table(paste(argv[1],"/peptide.tsv",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",row.names=1,check.name=FALSE)
peptide_table=peptide_table[!grepl("Cont",peptide_table$Protein),] ##remove decoy

reference_location=grep(argv[3],colnames(peptide_table))[1]
print(paste("Location of TMT reference channel:",reference_location,sep=" "))

peptide_inf=peptide_table[,c(1:(ncol(peptide_table)-tmt_plex))]   ###peptide information
peptide_intensity=peptide_table[,c((ncol(peptide_table)-tmt_plex+1):ncol(peptide_table))]  ###peptide intensity
peptide_ratio=log2(peptide_intensity/peptide_table[,reference_location])

peptide_inf[,"sorted_protein_list"]=NA  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"sorted_gene_list"]=NA  ####remove gene name. Genes from decoy are excluded
peptide_inf[,"number_of_mapped_protein"]=0  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"number_of_mapped_gene"]=0  ####remove gene name. Remove decoy if has both target ad decoy
peptide_inf[,"protein_number_of_gene"]=0  ####remove gene name. Remove decoy if has both target ad decoy

for(j in 1:nrow(peptide_inf))
{
  protein_list=paste(peptide_inf[j,"Protein"],peptide_inf[j,"Mapped Proteins"],sep=",")
  protein_list=gsub(" ","",protein_list)
  peptide_inf[j,"sorted_protein_list"]=paste(sort(strsplit(protein_list,",")[[1]]),sep=",",collapse =",")
  total_number_of_proteins=strsplit(protein_list,",")[[1]]
  total_number_of_proteins=total_number_of_proteins[!grepl("rev_sp",total_number_of_proteins)&!grepl("Cont",total_number_of_proteins)] ##remove decoy and cont proteins for a peptide
  total_number_of_proteins=length(total_number_of_proteins)
  peptide_inf[j,"number_of_mapped_protein"]=total_number_of_proteins
  #print(total_number_of_proteins)
  gene_list=paste(peptide_inf[j,"Gene"],peptide_inf[j,"Mapped Genes"],sep=",")
  gene_list=sub(" ","",gene_list)
  peptide_inf[j,"sorted_gene_list"]=paste(sort(strsplit(gene_list,",")[[1]]),sep=",",collapse =",")
  gene_list_detail=strsplit(gene_list,",")[[1]]
  peptide_inf[j,"number_of_mapped_gene"]=length(gene_list_detail)
  if(length(gene_list_detail)==1)
  {
    gene_name=gene_list_detail[1]
    peptide_inf[j,"protein_number_of_gene"]=gene_protein_number[gene_name,1]
  }
}

peptide_ratio=format(peptide_ratio,scientific = TRUE)
combined_results=cbind(peptide_inf,peptide_ratio)
write.table(rbind(peptide=colnames(combined_results),combined_results),file=paste(argv[1],"/processed-protein-and-ratio.txt",sep=""),sep="\t",row.names=TRUE,col.names = FALSE,quote=FALSE)



