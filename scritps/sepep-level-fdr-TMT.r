#!/usr/bin/env Rscript

argv <- commandArgs(TRUE)
if (length(argv) != 2) {
  cat("Usage: identify-sepep input_dir SEPEP_fdr\n")
  q(status = 1)
}
SEPEP_fdr=as.numeric(argv[2])

peptide = read.table(paste(argv[1],"/identify-sepep-and-quantified-by-median.txt",sep=""),header = TRUE,stringsAsFactors=FALSE, sep = "\t",check.name=FALSE)
peptide=peptide[order(peptide$Probability,decreasing = TRUE),]
peptide[,"FDR"]=0
  
true_positive=0;
false_positive=0;
for(j in 1:nrow(peptide))
{
   if(grepl("rev_",peptide$Protein[j]))
  {
    false_positive=false_positive+1
  }
  peptide[j,"FDR"]=false_positive/j
}
peptide=peptide[peptide$FDR<SEPEP_fdr,]
peptide=peptide[peptide$Gene %in% peptide$Gene,]
peptide=peptide[!grepl("rev_",peptide$Protein),]
peptide=peptide[,-ncol(peptide)]
peptide=peptide[!grepl("rev_sp",peptide$Protein),] ##remove decoy after FDR control
write.table(peptide,file=paste(argv[1],"/filtered-sepep-by-fdr.txt",sep=""),sep="\t",row.names=FALSE,col.names = TRUE,quote=FALSE)
  


