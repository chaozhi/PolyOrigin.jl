# Understand output

## Overview

PolyOrigin produces in the following files
* `oustem.log`: log file containing input argument values and computational time. 
* `outstem_plots`: a folder collecting all produced plots if isplot = true (false by default). 
* `outstem_parentphased.csv`: same as input genofile except that parents are phased
* `outstem_parentphased_corrected.csv`: same as input genofile except that parents are phased and corrected
* `outstem_maprefined.csv`: same as the above file except that map is refined
* `outstem_postdoseprob.csv`: same as the above file except that offspring genotypes are provided in terms of posterior dosage probabilites. 
* `outstem_valentsummary.csv`: estimation results on chromosomal pairing 
* `outstem_polyancestry.csv`: main results from ancestral inference
* `outstem_genoprob.csv`: a dataframe for posterior ancestral genotype probabilityes, which is also included in the above polyancestry file. 

We will focus on describing the polyancestry file, since it is the input file for downstream QTL mapping. 

## polyacnestry file

The output `polyancestryfile` is a CSV file containing a list of dataframes:

1. `designinfo`: parents for each subpopulation (F1 cross or selfing)
2. `parentinfo`: ploidy for each parent
3. `offspringinfo`: information for each offspring such as ploidy level, isoutlier (too many breakpoints) , and the posterior number of breakpoins in each linkage group
4. `delmarker`: list of deleted markers
5. `correction`: list of founder dosages being corrected in the step of parental phasing
6. `valentlist`: list of valent configurations for each full-sib family
7. `valentprob`: posterior probabilities of valent configurations for each offspring
8. `parentgeno`: phased genotpyes for all parents
9. `ancestralgenotype`: list of possible unphased origin-genotypes for each family. 
10. `genoprob`: marginal posterior probaiblities of unphased origin-genotypes at each marker for each individual

The unphased origin-genotypes of an offspring are defined in terms of the parental origins of its immediate parents, where the parental origins act as alleles. For example, if an offspring is produced by selfing a tetraploid parent and qudrivalent formation is allowed, there are 35 possible origin-genotypes: 1-1-1-1, 1-1-1-2, ..., 4-4-4-4, where 1, 2, 3, and 4 denote the parent origins (homologs) of the parent. If an offspring is produced by crossing tetraploid $P1$ with tetraploid $P2$ and qudrivalent formation is allowed, there are 100 possible origin-genotypes: 1-1-5-5, 1-1-5-6, ..., 4-4-8-8, where 1, 2, 3, and 4 denote the parent origins (homologs) of the 1st parent $P1$, and 5, 6, 7, and 8 denote the parent origins (homologs) of the 2nd parent $P2$. These possible origin-genotypes are listed in the dataframe `ancestralgenotype` for each subpopulation. 

In the `genoprob` dataframe, each cell denotes the posterior probabilities for a given offspring at a given marker. The sparse probability vector is represented in form of origin-genotype-indices => posterior-probabilites.  For example, 54|59|84=>0.995|0.005|0.001, the left part "54|59|84" denotes the indices 54, 59, and 84 of origin-genotypes 2-3-5-8, 2-3-7-8, and 3-4-5-8, respectively; and the right part denotes the corresonding posterior probabilities 0.995, 0.005, and 0.001, respectively. Each probability is round to three digits after the dicimal place and only non-zeros probabilities are listed. 