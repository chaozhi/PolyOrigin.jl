# Run polyorigin

Run function [`polyorigin`](@ref) for haplotype reconstruction in diploid or tetraploid multi-generation pedigrees. 

```@example
using PolyOrigin
workdir = abspath(pkgdir(PolyOrigin),"docs","run_polyorigin")
genofile = "geno.csv"
pedfile = "ped.csv"
outstem = "outstem"
polyorigin(genofile,pedfile;
    chrpairing_phase=22,    
    chrpairing=44,    
    isphysmap = false,
    isrefine = true,    
    workdir,
    outstem,    
)
filter(x->occursin(outstem,x), readdir(workdir))
```

Here `chrpairing_phase=22` denotes only bivalent formations in parental phasing, and `chrpairing=44` denotes bivalent and quadrivalent formations in ancestral inference. `isphysmap` specifies if the input map is a physical map, `isrefine` specifies if refinning the input map, `workdir` specifies the directory of input and output files, and `outstem` specifies the stem for output filenames.

