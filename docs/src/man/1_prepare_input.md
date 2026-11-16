# Prepare input

PolyOrigin requires two input files: genofile and pedfile. 

### Input genofile

The genofile should be in VCF format with extension ".vcf" or ".vcf.gz". 

A vcf file contains the physical map in the first three columns: #CHROM, POS, ID. A genetic map can be provided by adding keywords `PHYSCHROM` and `PHYSPOSBP` in the INFO column. Here is an example of info field: "PHYSCHROM=1;PHYSPOSBP=1000". 

!!! note

    To be compatible with the previous versions, The genofile may also be in CSV format with extension ".csv", satisfying 
    1. The first three columns must be `marker`,`chromosome`, and `position`. The aliases of `marker` include `markers` and `id`. The aliases of `chromosome` include `chr`, `chrom` and `linkagegroup`. The aliases of `position` include `pos`, `poscm`, and `cm`. They are case-insenstive. 
    2. The rest columns are sampled members that are in pedfile. 

    The position is parsed as physical map position in bp if option isphysmap = true, and otherwise it is parsed as genetic postion in cM. 

    There are four possible formats of genotypes:
    * DS: dosages range from 0, 1, ..., ploidy, and NA for missing dosage
    * GT: unphased genotypes in format of g1/g2/.../g(ploidy), where g(i)= 0 or 1 for i=1, ..., ploidy. Delimiter '|' is used instead for phased genotypes. 
    * AD: allelic depths in form of c1&c2, where c1 and c2 are the number of reads for alleles 1 and 2, respectively. Missing genotypes are denoted by 0&0. Delimiter `|` is depreciated. 
    * GP: probability vector in form of p(0)&p(1)&...&p(ploidy), where p(i) denotes the probability of observed data given dosage i = 0, ..., ploidy, and the probabilities are normalized so that their sum is 1. Delimiter `|` is depreciated. 

### Input pedfile

The pedfile must be in CSV format with extension ".csv".   

The pedfile stores pedigree information. It must contain at least four columns with names: individual, mother, father, ploidy.  The column names are case-insensitive. 


