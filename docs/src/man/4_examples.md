
# Examples

See [`PolyOrigin.jl`](https://github.com/chaozhi/PolyOrigin.jl). 

See [`examples`](https://github.com/chaozhi/PolyOrigin2_Examples).


[Haplotype reconstruction in a real 3x3 half-diallel potato population](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato4x_3F1/real_potato4x_F1.md)

[Haplotype reconstruction for a real diploid potato pedigree](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato2x_ped/run_polyorigin.md)


[Haplotype reconstruction for a real tetraploid potato pedigree](https://github.com/chaozhi/PolyOrigin2_Examples/tree/main/real_potato4x_ped/real_potato4x.md)


### Example: real potato

### Load `PolyOrigin`


Load `PolyOrigin` for non-parallel computation
````
using PolyOrigin
````

Or for multicore parallel computation
````
using Distributed
procs() < 4 && addprocs(4-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using PolyOrigin
````

where the processes consist of 1 main and 3 workers: by default, `polyorigin` performs 3 repeated phasing runs for each chromosome. 


### Set inputs

````
cd(@__DIR__) 
dataid = "potato4x" # or "potato2x"
genofile = PolyOrigin.get_datafile("potato4x_geno.vcf.gz")
pedfile = PolyOrigin.get_datafile("potato4x_ped.csv")
````

Set the directory of this file as the working directory when the output files will be saved.  `PolyOrigin` includes two example datasets: `potato2x` and `potato4x`.

### Run `polyorigin`

````
polyorigin(genofile,pedfile;
    isphysmap = true,
    recomrate = 1.25,
    chrsubset = 1:1, # anaylze only the first chromosome    
    isplot = true, 
    nplot_subpop = 1,
    outstem = dataid*"_out"
)
````

The options `isphysmap` and `recomrate` transform the input physical map (positions in base pairs) into genetic map with specified constant recombination rate. If `isplot` is true, plot pedigree, map comparison, and founder-origin genotype probabilities, and `nplot_subpop` specifies plotting for up to 1 offspring in each subpopulation or family. `outstem` specifies the stem of output files, and by default, `workdir` is given by the current directory.