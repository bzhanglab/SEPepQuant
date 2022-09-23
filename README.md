# SEPEPquant
 Tripartite graph modeling enables comprehensive protein isoform characterization in shotgun proteomics

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
The exact protein database used for FragPipe is required for SEPEPquant.


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




