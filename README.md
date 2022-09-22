# SEPEPquant
 Tripartite graph modeling enables comprehensive protein isoform characterization in shotgun proteomics

# Usage 
Curently, SEPEPquant only support FragPipe processed Label free and TMT data.  

## Installation

SEPEPquant does not need setup or installation. It has been test on Linux and Windows.

## Testing data set

### TMT11   data set from 
```r
perl SEPEPquant.pl --database protein_database\GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants_with_decoys.fa --fdr 0.01 --quant TMT --plex 11 --RefTag Mix --input testing_data\TMT --output testing_data_output_TMT
```
### Label free data set from 
```r
perl SEPEPquant.pl --database protein_database\GRCh38_latest_protein_NP_YP_XP.changeHeaderFormatUniprot.maxquant_contaminants_with_decoys.fa --fdr 0.01 --quant LF --input testing_data\Label_free --output testing_data_output_LF
```


