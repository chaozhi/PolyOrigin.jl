
# Examples

See [`PolyOrigin.jl`](https://github.com/chaozhi/PolyOrigin.jl). 

See [`examples`](https://github.com/chaozhi/PolyOrigin2_Examples).


[Haplotype reconstruction in a real 3x3 half-diallel potato population](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato4x_3F1/real_potato4x_F1.md)

[Haplotype reconstruction for a real diploid potato pedigree](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato2x_ped/run_polyorigin.md)


[Haplotype reconstruction for a real tetraploid potato pedigree](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato4x_ped/real_potato4x.md)


### Real potato

### Load `PolyOrigin` and set inputs

````
using PolyOrigin 
cd(@__DIR__) 
dataid = "potato4x" # or "potato2x"
genofile = PolyOrigin.get_datafile("potato4x_geno.vcf.gz")
pedfile = PolyOrigin.get_datafile("potato4x_ped.csv")
````

Load the package and set the directory of this file as the working directory when the output files will be saved.  `PolyOrigin` includes two example datasets: `potato2x` and `potato4x`.

### Run `polyorigin`

````
polyorigin(genofile,pedfile;
    isphysmap = true,
    recomrate = 1.25,
    chrsubset = 1:1, # anaylze only the first chromosome
    maxrun = 1, # hide
    isplot = true, 
    nplot_subpop = 1,
    outstem = dataid*"_out"
)
````

The options `isphysmap` and `recomrate` transform the input physical map (positions in base pairs) into genetic map with specified constant recombination rate. If `isplot` is true, plot pedigree, map comparison, and founder-origin genotype probabilities, and `nplot_subpop` specifies plotting for up to 1 offspring in each subpopulation or family. `outstem` specifies the stem of output files, and by default, `workdir` is given by the current directory.