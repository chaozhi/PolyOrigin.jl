# Run polyOrigin

Run function [`polyOrigin`](@ref) for haplotype reconstruction, which consists of two major steps: parental phasing and ancestral inference for each offspring.

Here is a simulated example of a tetraploid population, consists of four subpopulations: diallel crosses among 3 parents and one selfing of the first parent. The genofile contains genetic map with 100 SNPs and dosage data for 3 parents and 100 offspring.

```@example
using PolyOrigin
workdir = joinpath(pkgdir(PolyOrigin),"docs","run_polyOrigin")
genofile = "geno_disturbedmap.csv"
pedfile = "ped.csv"
outstem="outstem"
polyancestry = polyOrigin(genofile,pedfile;
    chrpairing=44,    
    refinemap=true,
    refineorder=true,
    workdir,
    outstem
)
filter(x->occursin(outstem,x), readdir(workdir))
```
Here `chrpairing=44` denotes bivalent and quadrivalent formations in ancestral inference. `refinemap` and `refineorder` specify to refine inter-marker distances and local marker ordering. `workdir` specifies the directory of input and output files, and `outstem` specifies the filename stem for output files.

The polyOrigin function produces the following outputfiles:
* `outstem.log`: log file saves messages that are printed on console.
* `outstem_maprefined.csv`: same as the input geno file, except that input marker map is refined.
* `outstem_parentphased.csv`: same as the input geno file, except that parental genotypes are phased.
* `outstem_parentphased_corrected.csv`: the phased parent genotypes are further corrected.
* `outstem_polyancestry.csv`: saves the returned variable polyancestry. See also [`polyOrigin`](@ref) for the description of dataframes that are saved in this output file.
* `outstem_genoprob.csv`: a concise version of the above file, including genetic map, phased parental genotypes, and posterior genotype probabilities.
* `outstem_postdoseprob.csv`: same as the input geno file, except that the parent genotypes are phased and the offspring genotypes are given by the posterior probabilities of dosage.

Here  `outstem_maprefined.csv`, `outstem_parentphased.csv`,  `outstem_parentphased_corrected.csv`, or `outstem_postdoseprob.csv` may be iteratively used as the input genofile.

# Visualize conditional probability

The `outstem_polyancestry.csv` can be read back by the function [`readPolyAncestry`](@ref), and ancestral conditional probability can be visualized by [`plotCondprob`](@ref) or [`animCondprob`](@ref).

```@setup setdir
using PolyOrigin
workdir = joinpath(pkgdir(PolyOrigin),"docs","run_polyOrigin")
```

```@example setdir
polyancestry = readPolyAncestry("outstem_polyancestry.csv",workdir=workdir)
truefile = "true.csv"
truegeno = readTruegeno!(truefile,polyancestry,workdir=workdir)
animCondprob(polyancestry,truegeno=truegeno)
```

# Evaluate estimated map

Compare input map and estimated map with true map

```@setup settrue
using PolyOrigin
workdir = joinpath(pkgdir(PolyOrigin),"docs","run_polyOrigin")
genofile = "geno_disturbedmap.csv"
pedfile = "ped.csv"
polyancestry = readPolyAncestry("outstem_polyancestry.csv",workdir=workdir)
truefile = "true.csv"
truegeno = readTruegeno!(truefile,polyancestry,workdir=workdir)
```

```@example settrue
polygeno=readPolyGeno(genofile,pedfile,workdir=workdir)
fig = plotMapComp(truegeno.truemap,polygeno.markermap,
    xlabel="True position (cM)",
    ylabel="Input position (cM)")
fig2 = plotMapComp(truegeno.truemap, polyancestry.markermap,
    xlabel="True position (cM)",
    ylabel="Estimated position (cM)")
using Plots
plot(fig,fig2)
```
where tau is the Kendall rank correlation between true map and comparing map.

```@setup deloutput
using PolyOrigin
workdir = joinpath(pkgdir(PolyOrigin),"docs","run_polyOrigin")
outstem="outstem"
outfiles = filter(x->occursin(outstem,x), readdir(workdir;join=true))
rm.(outfiles)
```
