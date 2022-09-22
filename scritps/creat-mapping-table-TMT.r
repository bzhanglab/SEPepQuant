#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 3) {
  cat("Usage:creat-sepep-mapping-table input_dir output_dir tmt_plex\n")
  q(status = 1)
}
tmt_plex=as.numeric(argv[3])

sample_list=dir(argv[1])
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  print(paste(i,tmt_name),sep="\t")
  peptide = read.table(paste(argv[1],tmt_name,"filtered-sepep-by-fdr.txt",sep="/"),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE,na.strings = "NA")
  peptide = peptide[!is.na(peptide$peptides),]  ##remove na
  peptide = peptide[!is.na(peptide$sorted_protein_list),]  ##remove na
  peptide = peptide[!grepl("rev_",peptide$Protein),]  ##remove decoy only
  peptide = peptide[,c(1:(ncol(peptide)-tmt_plex))]
  if(i==1)
  {
    cliques_for_whole_cohort=peptide
  } else {
    cliques_for_whole_cohort=rbind(cliques_for_whole_cohort,peptide)
  }
}

cliques_for_whole_cohort=cliques_for_whole_cohort[!is.na(cliques_for_whole_cohort$sorted_protein_list),]
cliques_for_whole_cohort=cliques_for_whole_cohort[!is.na(cliques_for_whole_cohort$sorted_gene_list),]
cliques_for_whole_cohort=cliques_for_whole_cohort[!is.na(cliques_for_whole_cohort$Gene),]
cliques_for_whole_cohort=cliques_for_whole_cohort[!is.na(cliques_for_whole_cohort$number_of_mapped_protein),]

gene_clique_number=as.data.frame(table(cliques_for_whole_cohort$sorted_gene_list))
rownames(gene_clique_number)=gene_clique_number$Var1
gene_clique_number$Freq=2

cliques_for_whole_cohort=cliques_for_whole_cohort[!duplicated(cliques_for_whole_cohort$sorted_protein_list),]

cliques_for_whole_cohort_multiple_gene=cliques_for_whole_cohort[cliques_for_whole_cohort$number_of_mapped_gene>1,]
cliques_for_whole_cohort_whole_gene=cliques_for_whole_cohort[cliques_for_whole_cohort$number_of_mapped_gene==1&(cliques_for_whole_cohort$number_of_mapped_protein==cliques_for_whole_cohort$protein_number_of_gene),]
cliques_for_whole_cohort_single_protein=cliques_for_whole_cohort[cliques_for_whole_cohort$number_of_mapped_gene==1&cliques_for_whole_cohort$number_of_mapped_protein==1&cliques_for_whole_cohort$protein_number_of_gene>1,]
cliques_for_whole_cohort_multiple_protein=cliques_for_whole_cohort[cliques_for_whole_cohort$number_of_mapped_gene==1&cliques_for_whole_cohort$number_of_mapped_protein>1&(cliques_for_whole_cohort$number_of_mapped_protein<cliques_for_whole_cohort$protein_number_of_gene),]

#setdiff(cliques_for_whole_cohort$sorted_protein_list,c(cliques_for_whole_cohort_multiple_gene$sorted_protein_list,cliques_for_whole_cohort_whole_gene$sorted_protein_list,cliques_for_whole_cohort_single_protein$sorted_protein_list,cliques_for_whole_cohort_multiple_protein$sorted_protein_list))
cliques_for_whole_cohort_multiple_gene[,"SEPEP_label"]="Multiple_SEPEP"
for(i in 1:nrow(cliques_for_whole_cohort_multiple_gene))
{
  cliques_for_whole_cohort_multiple_gene$SEPEP_label[i]=paste(cliques_for_whole_cohort_multiple_gene$SEPEP_label[i],".",i,"_C5",sep="")
  #print(cliques_for_whole_cohort_multiple_gene$SEPEP_label[i])
}

cliques_for_whole_cohort_whole_gene[,"SEPEP_label"]=cliques_for_whole_cohort_whole_gene$Gene
cliques_for_whole_cohort_whole_gene=cliques_for_whole_cohort_whole_gene[!is.na(cliques_for_whole_cohort_whole_gene$number_of_mapped_protein),]
for(i in 1:nrow(cliques_for_whole_cohort_whole_gene))
{
  #print(cliques_for_whole_cohort_whole_gene$number_of_mapped_protein[i])
  if(cliques_for_whole_cohort_whole_gene$number_of_mapped_protein[i]==1)
  {
    cliques_for_whole_cohort_whole_gene[i,"SEPEP_label"]=paste(cliques_for_whole_cohort_whole_gene$Gene[i],"_SEPEP.1_C1",sep="")
  } else {
    cliques_for_whole_cohort_whole_gene[i,"SEPEP_label"]=paste(cliques_for_whole_cohort_whole_gene$Gene[i],"_SEPEP.1_C4",sep="")
  }
}

cliques_for_whole_cohort_single_protein[,"SEPEP_label"]=cliques_for_whole_cohort_single_protein$Gene
cliques_for_whole_cohort_single_protein=cliques_for_whole_cohort_single_protein[!is.na(cliques_for_whole_cohort_single_protein$number_of_mapped_protein),]
#cliques_for_whole_cohort_single_protein$SEPEP_label=make.unique(cliques_for_whole_cohort_single_protein$SEPEP_label)
for(i in 1:nrow(cliques_for_whole_cohort_single_protein))
{
  cliques_for_whole_cohort_single_protein[i,"SEPEP_label"]=paste(cliques_for_whole_cohort_single_protein[i,"SEPEP_label"],"_SEPEP.",gene_clique_number[cliques_for_whole_cohort_single_protein$Gene[i],"Freq"],"_C2",sep="")
  gene_clique_number[cliques_for_whole_cohort_single_protein$Gene[i],"Freq"]=gene_clique_number[cliques_for_whole_cohort_single_protein$Gene[i],"Freq"]+1
}

cliques_for_whole_cohort_multiple_protein[,"SEPEP_label"]=cliques_for_whole_cohort_multiple_protein$Gene
cliques_for_whole_cohort_multiple_protein=cliques_for_whole_cohort_multiple_protein[!is.na(cliques_for_whole_cohort_multiple_protein$number_of_mapped_protein),]
#cliques_for_whole_cohort_multiple_protein$SEPEP_label=make.unique(cliques_for_whole_cohort_multiple_protein$SEPEP_label)
for(i in 1:nrow(cliques_for_whole_cohort_multiple_protein))
{
  cliques_for_whole_cohort_multiple_protein[i,"SEPEP_label"]=paste(cliques_for_whole_cohort_multiple_protein[i,"SEPEP_label"],"_SEPEP.",gene_clique_number[cliques_for_whole_cohort_multiple_protein$Gene[i],"Freq"],"_C3",sep="")
  gene_clique_number[cliques_for_whole_cohort_multiple_protein$Gene[i],"Freq"]=gene_clique_number[cliques_for_whole_cohort_multiple_protein$Gene[i],"Freq"]+1
}

whole_mapping_table=rbind(cliques_for_whole_cohort_multiple_gene,cliques_for_whole_cohort_whole_gene,cliques_for_whole_cohort_single_protein,cliques_for_whole_cohort_multiple_protein)

mapping_table_columns=ncol(whole_mapping_table)
whole_mapping_table=whole_mapping_table[,c(mapping_table_columns,mapping_table_columns-4,mapping_table_columns-3,mapping_table_columns-2,mapping_table_columns-5,mapping_table_columns-6)]

write.table(whole_mapping_table,file=paste(argv[2],"/sepep_mapping_table.txt",sep=""),col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE)

whole_mapping_table=whole_mapping_table[,c(1,6)]
rownames(whole_mapping_table)=whole_mapping_table$sorted_protein_list

sample_list=dir(argv[1])
for(i in 1:length(sample_list))
{
  tmt_name=sample_list[i]
  #print(paste(i,tmt_name),sep="\t")
  peptide = read.table(paste(argv[1],tmt_name,"/filtered-sepep-by-fdr.txt",sep="/"),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
  sub_mapping_table=whole_mapping_table[peptide$sorted_protein_list,]
  peptide_part_1=peptide[,c(1:(ncol(peptide)-tmt_plex))]
  peptide_part_2=peptide[,c((ncol(peptide)-tmt_plex+1):ncol(peptide))]
  combined_table=cbind(peptide_part_1,sub_mapping_table,peptide_part_2)
  #combined_table=combined_table[,-24]
  write.table(combined_table,file=paste(argv[1],tmt_name,"/sepepe-with-class.txt",sep="/"),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
  
}















