# Prepare input

PolyOrigin requires two input files: genofile and pedfile. 

### Input genofile

The genofile should be in VCF format with extension ".vcf" or ".vcf.gz". A vcf file contains the physical map in the first three columns: #CHROM, POS, ID. A genetic map can be provided by adding keywords `PHYSCHROM` and `PHYSPOSBP` in the INFO column. Here is an exmaple of info field: "PHYSCHROM=1;PHYSPOSBP=1000". 

To be compatible with the previous versions, The genofile may also be in CSV format with extension ".csv", satisfying 

1. The first three columns must be `marker`,`chromosome`, and `position`. The aliases of `marker` include `markers` and `id`. The aliases of `chromosome` inclue `chr`, `chrom` and `linkagegroup`. The aliases of `position` include `pos`, `poscm`, and `cm`. They are case-insenstive. 
2. The rest columns are sampled members that are in pedfile. 

The position is parsed as physical map position in bp if option isphysmap = true, and otherwise it is parsed as genetic postion in cM. 

### Input pedfile

The pedfile must be in CSV format with extension ".csv".   

The pedfile stores pedigree information. It must contain at least four columns with names: individual, mother, father, ploidy.  The column names are case-insensitive. 



