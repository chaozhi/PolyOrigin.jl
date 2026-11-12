# Run `polyorigin`

Run function [`polyorigin`](@ref) for haplotype reconstruction in diploid or tetraploid multi-generation pedigrees. 

```@example 1
using PolyOrigin 
workdir = abspath(pkgdir(PolyOrigin),"docs","run_polyorigin") 
dataid = "sim_4x"
genofile = dataid*"_geno.vcf.gz"   
pedfile = dataid*"_ped.csv"
truegenofile = abspath(workdir,dataid*"_truegeno.vcf.gz")
trueorigfile = abspath(workdir,dataid*"_trueorig.vcf.gz")
outstem = dataid*"_out"
phasedgeno, polyancestry, joinedancestry = polyorigin(genofile,pedfile;
    chrpairing_phase=22,    
    chrpairing=44,        
    isphysmap = false,
    isrefine = true,    
    isplot = false, 
    maxrun = 1, # hide
    refhapfile = truegenofile, 
    trueorigfile = trueorigfile,        
    workdir,
    outstem,    
)
filter(x->occursin(outstem,x), readdir(workdir))
```

where `chrpairing_phase=22` denotes only bivalent formations in parental phasing, and `chrpairing=44` denotes bivalent and quadrivalent formations in ancestral inference. `isphysmap` specifies if the input map is a physical map, `isrefine` specifies if refining the input map, `workdir` specifies the directory of input and output files, and `outstem` specifies the stem for output filenames.

The optional true files are used only for calibrating absolute phase (that is the ordering of homologies) for calculating accuracy in the following. 

```@example 1
acc = calaccuracy!(truegenofile,trueorigfile, joinedancestry)
df = PolyOrigin.acc2df(acc) # hide
show(df, row_labels=nothing,eltypes=false, summary=false) # hide
```
where the accuracies are for haplotype phasing and ancestral inference and for founders, intermediate parents, and leaf offspring. The permutation (or switch for diploid) phasing error is given by npermerr/npermoppty. 

```@example 1
trueorig = readtruegeno!(trueorigfile,joinedancestry)        
plotcondprob(joinedancestry;                 
    ismidparent = true, 
    individual = 4,    
    trueorig,
)
```
where the gray crosses ("x") denote true founder origins for the 4-th parent that is an intermediate parent. 


```@example 1
rm.(filter(x->!isdir(x) && occursin(outstem,x), readdir(workdir; join=true)))
nothing # hide
```
Clean up output files. 