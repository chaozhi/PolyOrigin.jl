# Prepare input

PolyOrigin requires two input files: genofile and pedfile. 

## Input genofile

The genofile must be in VCF format with extension ".vcf" or ".vcf.gz". 

A vcf file contains the physical map in the first three columns: #CHROM, POS, ID. A genetic map can be provided by adding keywords `PHYSCHROM` and `PHYSPOSBP` in the INFO column. Here is an exmaple of info field: "PHYSCHROM=1;PHYSPOSBP=1000". 

## Input pedfile

The pedfile must be in CSV format with extension ".csv".   

The pedfile stores pedigree information. It must contain at least four columns with names: individual, mother, father, ploidy.  The column names are case-insensitive. 



