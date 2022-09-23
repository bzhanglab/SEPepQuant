# SEPEPquant
 Tripartite graph modeling enables comprehensive protein isoform characterization in shotgun proteomics

[<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/protein-and-peptide-distribution.jpg" width=400 class="center">](https://github.com/bzhanglab/SEPEPquant)

Among the 19449 protein coding genes annotated in a RefSeq database, 14698 (75.6%) have more than one protein isoforms, and 3409 (17.5%) have 10 or more protein isoforms (Fig. 1a). Most of isoforms from the same gene have very high sequence similarity (>90%, Fig. 1b). However, among the 11809 genes with three or more protein isoforms, 6165 (52.2%) have at least one pair of isoforms with a sequence similarity lower than 90%, or an average of one amino acid difference in every 10 amino acids, suggesting the possibility to identify isoform-discriminating peptide sequences for a substantial number of genes.
To further assess the challenge and opportunities of isoform characterization using shotgun proteomics, we performed in silico trypsin digestion of the RefSeq protein database to generate fully tryptic peptides with length 7 to 50 and no missed cleavage. Among the 1,883,206 resulting peptide sequences, 2.8% could be associated to multiple genes (i.e., multi-genes peptides), 13.6% to genes with a single protein isoform (i.e., single isoform peptides), and 83.5% to genes with more than one isoforms (i.e., multi-isoforms peptides). Within the group of multi-isoforms peptides, around half could be mapped to all protein isoforms of a gene and thus providing no information for isoform discrimination (i.e., non-discriminative peptides); however, another half, or 246,615 peptides, could be uniquely mapped to one isoform (i.e., fully discriminative peptides) or a subset of isoforms (i.e., partially discriminative peptides) (Fig. 1c).


[<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/parsimony-selection.jpg" width=400 class="center">](https://github.com/bzhanglab/SEPEPquant)

[<img src="https://github.com/bzhanglab/SEPEPquant/blob/main/doc/sepep-quantification.jpg" width=400 class="center">](https://github.com/bzhanglab/SEPEPquant)

# Usage 
Curently, SEPEPquant only support [FragPipe](https://fragpipe.nesvilab.org/) processed Label free and TMT data.  

```shell
FragPipe_output
├── Sample1/TMT1
│   ├── peptide.tsv
│
├── Sample2/TMT2
│   ├── peptide.tsv
│
├── Sample3/TMT3
│   ├── peptide.tsv
│
├── Sample4/TMT4
│   ├── peptide.tsv

```

## Installation

SEPEPquant does not need setup or installation. It has been test on Linux and Windows.

## Protein database
The exact protein database used for FragPipe is required for SEPEPquant. The database should follow the [philosopher format rules](https://github.com/Nesvilab/philosopher/wiki/How-to-Prepare-a-Protein-Database#header-formatting) and in [UniProt](https://www.uniprot.org/help/fasta-headers) format.


## Testing data set

### TMT11   data set from 
```r
perl SEPEPquant.pl --database protein_database\GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants_with_decoys.fa --fdr 0.01 --quant TMT --plex 11 --RefTag Mix --input testing_data\TMT --output testing_data_output_TMT
```
### Label free data set from 
```r
perl SEPEPquant.pl --database protein_database\GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants_with_decoys.fa --fdr 0.01 --quant LF --input testing_data\Label_free --output testing_data_output_LF
```

### Parameters

--database: protein database used for database searching

--quant: LF or TMT

--plex: TMT plex, required if --quant is TMT

--RefTag: Tag of TMT reference channel, required if --quant is TMT

--input: folder of FragPipe output

--output: output filder

--help: this information




