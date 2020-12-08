# Prepare input

The input files of function polyOrigin are two text files: genofile and pedfile. The input files are  in CSV format; set the option `delimchar` for a different delimiter.

## Input genofile

The genofile stores genetic map, parent genotypes, and offspring genotypes. Click to download [zipped example genofile](geno.csv.zip).
For example, a genofile looks like

marker | chromosome | position | ind1 | ind2 | ind3 | ind4 | ...
--- | --- | --- | --- | --- | --- | --- | ---
snp1 | 1 | 0.14 | 0 | 2 | 1 | 4 | ...
snp2 | 1 | 0.16 | 4 | 0 | NA | 2 | ...
snp3 | 1 | 0.21 | NA | 3 | 0 | 1 | ...

The genofile must pass the following check list
* Column 1: SNP IDs must be unique.
* Column 2: SNPs must be grouped by chromosome IDs.
* Column 3: Positions of markers within a
  chromosome must be non-decreasing. The position unit must be base-pair for physical map and centiMorgan for genetic map.
* Row 1, Col 4-end: Individual IDs must be unique, and must be in *pedfile*.
* Row 2:end, Col 4-end: all genotypes of parents must take on one of the following formats 1-4, and all genotypes of offspring must take one of the following formats 1-3.

List of four possible format of genotypes:
1. `dosage`: ranges from 0, 1, ..., ploidy, and NA for missing dosage;
2. `readcount`: c1|c2, where c1 and c2 are the number of reads for alleles 1 and 2, respectively. Missing genotypes are denoted by 0|0
3. `probability`: p(0)|p(1)|...|p(ploidy), where p(i) denotes the probability of observed data given dosage i = 0, ..., ploidy, and the probabilities are normalized so that their sum is 1.
4. `phasedgeno`: g1|g2|...|g(ploidy),  where g(i)=1 or 2 for i=1, ..., ploidy.

## Input pedfile

The pedfile stores pedigree information. Click to download [zipped example pedfile](ped.csv.zip). For example, a pedfile looks like

individual | population | motherid | fatherid | ploidy
--- | --- | --- | --- | ---
P1 | 0 | 0 | 0 | 4
P3 | 0 | 0 | 0 | 4
P3 | 0 | 0 | 0 | 4
offspring1 | pop1 | P1 | P2 | 4
offspring2 | pop1 | P1 | P2 | 4
offspring3 | pop2 | P1 | P3 | 4
offspring4 | pop2 | P1 | P3 | 4
offspring5 | pop3 | P2 | P3 | 4
offspring6 | pop4 | P3 | P3 | 4

The pedigree contains three founders (parents), two offspring from the cross
between parents 1 and 2, two offspring from the cross between parents 1 and 3,
one offspring from the cross between parents 2 and 3, and one offspring from the selfing of
parent 3.

The pedfile must pass the following check list
* Column 1: individual IDs must be unique, and must be in *genofile*.
* Column 2: Unique ID for each sub-population (F1 cross or selfing). The sub-population of founders must be the same.
* Column 3-4: parentID of founders must be 0.
* Column 5: ploidy must be 4. TODO for n=2, 4, 6, and 8.
