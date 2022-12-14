# Introduction
### SEPepQuant: a graph theory-based approach enables comprehensive protein isoform characterization in shotgun proteomics

### View our manuscript on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.03.515027v1)


 [<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/protein-and-peptide-distribution.jpg" width=500 class="center">](https://github.com/bzhanglab/SEPEPquant)

Among the 19449 protein coding genes annotated in a RefSeq database, 14698 (75.6%) have more than one protein isoforms, and 3409 (17.5%) have 10 or more protein isoforms (Fig. a). Most of isoforms from the same gene have very high sequence similarity (>90%, Fig. b). However, among the 11809 genes with three or more protein isoforms, 6165 (52.2%) have at least one pair of isoforms with a sequence similarity lower than 90%, or an average of one amino acid difference in every 10 amino acids, suggesting the possibility to identify isoform-discriminating peptide sequences for a substantial number of genes.
To further assess the challenge and opportunities of isoform characterization using shotgun proteomics, we performed in silico trypsin digestion of the RefSeq protein database to generate fully tryptic peptides with length 7 to 50 and no missed cleavage. Among the 1,883,206 resulting peptide sequences, 2.8% could be associated to multiple genes (i.e., multi-genes peptides), 13.6% to genes with a single protein isoform (i.e., single isoform peptides), and 83.5% to genes with more than one isoforms (i.e., multi-isoforms peptides). Within the group of multi-isoforms peptides, around half could be mapped to all protein isoforms of a gene and thus providing no information for isoform discrimination (i.e., non-discriminative peptides); however, another half, or 246,615 peptides, could be uniquely mapped to one isoform (i.e., fully discriminative peptides) or a subset of isoforms (i.e., partially discriminative peptides) (Fig. c).


 [<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/parsimony-selection.jpg" width=500 class="center">](https://github.com/bzhanglab/SEPEPquant)

Peptides shared by multiple genes or multiple protein isoforms of the same gene complicate protein inference and quantification. Based on Occam???s razor or the principle of parsimony, the current best practice in the proteomics field is to collapse proteins with the same or subset of supporting peptides into a minimal list of protein groups, and for quantitative rollup, peptides shared by multiple proteins are assigned only to the ones with the most identification evidence. Although practically useful, this parsimonious approach greatly limits the potential for protein isoform characterization.

 [<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/sepep-quantification.jpg" width=700 class="center">](https://github.com/bzhanglab/SEPEPquant)

The TGM approach involves four major steps (Fig. a-d). First, a tripartite graph is built with three vertices sets representing all peptides identified in a study (Pep1-Pep12), proteins of the peptides can be mapped to (Pro1.1-Pro3.1), and host genes of the proteins (Gene1-Gene3), respectively, and the vertices are connected by edges indicating their mapping relationships (Fig. a). Second, Peptides that are connected to exactly the same set of protein vertices are grouped together (Pep2 and Pep4; Pep 6 and Pep7; Pep9 and Pep9; Pep11 and Pep12), and each set of peptide vertices is defined as a Structurally Equivalent PEPtide set (SEPEP) (Fig. b). Third, the picked FDR approach29 is used to estimate target-decoy based FDR at the SEPEP level. Specifically, a SEPEP is considered a target hit if the highest scoring peptide in the SEPEP is from a forward protein sequence, and a decoy hit if the highest scoring peptide is from a decoy protein sequence. SEPEPs with FDR >0.01 are excluded from further analysis. Finally, the remaining SEPEPs are classified into 5 classes based on their patterns of connections to source proteins and genes in the tripartite graph (Fig. d). Class 1 through 5 correspond to single isoform peptide sets, fully discriminative peptide sets, partially discriminative peptide sets, non-discriminative peptide sets, and multi-genes peptide sets, respectively. Class 1 to 4 SEPEPs are labeled with gene name followed by SEPEP order within the gene and SEPEP class type, e.g., Gene1_SEPEP.1_C1. Class 5 SEPEPs are labeled with ???Multiple??? followed by SEPEP order across the whole study and C5, e.g., Multiple_SEPEP.1_C5. 
In contrast to existing methods that make protein inference and then use protein groups or genes as the reporting and quantification unit, the TGM approach uses SEPEP as the reporting and quantification unit (Fig. e). Both methods share the same database searching, PSM FDR control, and peptide FDR control protocols. SEPEP identification and SEPEP level FDR control are performed in parallel to the standard protein inference and protein level FDR control. Finally, the same algorithm can be used to report quantification at SEPEP, gene, and protein group levels. 


# Installation

SEPepQuant does not need setup or installation. While it required Rscript in your PATH. It has been test on Linux and MacOS.

Download SEPepQuant:
```sh
git clone https://github.com/bzhanglab/SEPepQuant.git
```

# Running SEPepQuant

SEPepQuant runs from a command line as follows:

```sh

Usage: perl SEPepQuant.pl [Options]

Parameters:
--database  protein database used for database searching
--quant LF or TMT
--plex  TMT plex, required if "--quant" is TMT
--RefTag    Tag of TMT reference channel, required if "--quant" is TMT
--input folder of FragPipe output
--output    output folder
--help  this information
```
### Notes of parameters:

--database:

The exact protein database used for FragPipe is required for SEPEPquant. The database should follow the [philosopher format rules](https://github.com/Nesvilab/philosopher/wiki/How-to-Prepare-a-Protein-Database#header-formatting) and in [UniProt](https://www.uniprot.org/help/fasta-headers) format. An example can be found in the [protein_database](https://github.com/bzhanglab/SEPEPquant/protein_database) folder.

--input:

Curently, SEPEPquant only support [FragPipe](https://fragpipe.nesvilab.org/) processed Label free and TMT data. The input files should be organized as the following format:

```shell
FragPipe_output
????????? Sample1/TMT1
??????? ????????? peptide.tsv
???
????????? Sample2/TMT2
??????? ????????? peptide.tsv
???
????????? Sample3/TMT3
??????? ????????? peptide.tsv

```

# Outputs

### 1. gene-protein-statistic.txt
This file contains numbers of protein(s) of each gene in the protein database.
* Column 1: gene name
* Column 2: # of proteins of a gene

### 2. sepep_mapping_table.txt
This file descripts mappings between SEPEPs, host gene(s) and host protein(s).
* SEPEP_label: name and class of a SEPEP
* number_of_mapped_protein: Number of proteins belong to this SEPEP
* number_of_mapped_gene: Number of genes belong to this SEPEP
* protein_number_of_gene: Number of proteins belongs to host gene(s), 0 for multiple gene SEPEPs
* sorted_gene_list: Sorted gene list
* sorted_protein_list: Sorted protein list

### 3. sepep quantification files
TMT data: It will report both raw (sepep_matrix_raw.txt) and median centered (sepep_matrix_MD.txt) Log2(TMT ratio) for TMT data.
 
Label free data: It will report both MS1 intensity and spectral count in a same file.


# Example usage

Enter the directory where you downloaded SEPepQuant and type the following command lines in your terminal:


### [Example 1](https://github.com/bzhanglab/SEPepQuant/tree/main/testing_data/Label_free): 3 samples from the [HCC-label free](https://www.nature.com/articles/s41586-019-0987-8) dataset
```r
perl SEPepQuant.pl --database protein_database/GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants.fa --fdr 0.01 --quant LF --input testing_data/Label_free --output testing_data_output_LF
```

### [Example 2](https://github.com/bzhanglab/SEPepQuant/tree/main/testing_data/TMT): 3 TMT plex from the [TCC-TMT](https://www.sciencedirect.com/science/article/pii/S0092867419310037) dataset
```r
perl SEPepQuant.pl --database protein_database/GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants.fa --fdr 0.01 --quant TMT --plex 11 --RefTag Mix --input testing_data/TMT --output testing_data_output_TMT
```

By above commond lines, SEPepQuant will identify and quantify SEPEPs for the example label free and TMT data in the [testing_data](https://github.com/bzhanglab/SEPepQuant/tree/main/testing_data/) folder. SEPepQuant will creat output folders testing_data_output_LF and testing_data_output_TMT for label free and TMT data, respectively. The running times for label free and TMT data are about 5 and 30 minutes. The same outputs are expected which are shown in the [testing_data_output](https://github.com/bzhanglab/SEPepQuant/tree/main/testing_data_output/) folder. 






