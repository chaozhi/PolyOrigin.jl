# Run `polyorigin` in cmd

First change into a working directory where genofile and pedfile are.  Then run the following command line: 

```
>polyorigin -g sim_4x_geno.vcf.gz -p sim_4x_ped.csv ^    
    --outstem sim_4x_out --nworker 3
```
where the caret(`^`) is a multi-line character for Windows CMD, and the backslash (`\`) should be used instead for Linux/macOS shells. 

By default, `polyorigin` analyzes all chromosomes (only one chromosome was simulted).  Here parallel computation launches workers; by default, 3 phasing runs for each chromosome. If `nworker` == 1, non-parallel computation. 


For help, run
```
>polyorigin --help
```
