# PolyOrigin.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chaozhi.github.io/PolyOrigin.jl/dev)
[![Build Status](https://github.com/chaozhi/PolyOrigin.jl/workflows/CI/badge.svg)](https://github.com/chaozhi/PolyOrigin.jl/actions)
[![Coverage](https://codecov.io/gh/chaozhi/PolyOrigin.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/chaozhi/PolyOrigin.jl)

A package for haplotype reconstruction in connected polyploid F1 populations

## Features

- Apply to connected F1 in tetraploid (TODO for higher ploidy levels)
- Apply to SNP array and GBS data
- Robust to dosage errors in SNP array data
- Robust to dosage uncertainties in low read depth GBS data
- Robust to errors in input genetic or physical map

## Installation
From the julia (>v1.5.0, 64-bit) REPL model, type `]` to enter the Pkg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Vignettes

[Haplotype reconstruction in a simulated tetraploid population from SNP array data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simarray/step3_tetraploid_simarray.md)

[Haplotype reconstruction in a simulated tetraploid population from GBS data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simgbs/step2_tetraploid_simgbs.md)

[Haplotype reconstruction in a real 3x3 half-diallel potato population](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_realpotato/tetraploid_realpotato.md)


## Usage

From the julia REPL model, run
```julia-repl
using PolyOrigin
```

## Help

From the julia REPL model, run
```julia-repl
?polyOrigin
```
where `?` enters Help REPL model.

`polyOrigin(genofile, pedfile, keyargs...)` performs parental phasing and ancestral inference from input files. Only ancestral inference is performed in the case of phased parents.

### Positional arguments

`genofile::AbstractString`: filename for genotypic data file.

`pedfile::AbstractString`:  filename for pedigree information.

### Keyword arguments

`delimchar::AbstractChar=','`:  text delimiter.

`missingstring::AbstractString="NA"`: string code for missing value.

`commentstring::AbstractString="#"`: rows that begin with commentstring will be ignored.

`isphysmap::Bool=false`: if true, input markermap is physical map, where
marker locations are in unit of base pair(bp).

`recomrate::Real=1`: recombination rate in unit of 1 cM/Mbp (centiMorgan per million base pair). Valid only if `isphysmap=true`.

`epsilon::Real=0.01`: genotyping error probability.

`seqerr::Real=0.001`: sequencing read error probability for GBS data.

`chrpairing_phase::Integer=22`: chromosome pairing in parental phasing, with 22 being only bivalent formations and 44 being bivalent and quadrivalent formations.

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome, with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=false`: if true, multicore computing over chromosomes.

`delmarker::Bool=true`: if true, delete markers during parental phasing.

`delsiglevel::Real=0.05`: significance level for deleting markers.

`maxstuck::Integer=5`: the max number of consecutive iterations that are rejected
in a phasing run.

`maxiter::Integer=30`: the max number of iterations in a phasing run.

`minrun::Integer=3`: if the number of phasing runs having the same parental phases
reaches minrun, phasing algorithm will stop before reaching the maxrun.

`maxrun::Integer=10`: the max number of phasing runs.

`byparent::Union{Nothing,Bool}=nothing`: if true, update parental phases
parent by parent; if false, update parental phases one subpopulation by subpopulation. The nothing denotes that it is true if a connected component is a single F1 cross, and false otherwise.

`byneighbor::Union{Nothing,Bool}=nothing`: if ture, udpate the combination of bivalent or multivalents in parents by their neighbors; if false, consider all the possible combinations. The nothing denotes that it is true if max ploidy>=6, and false otherwise.

`refhapfile::Union{Nothing,AbstractString} = nothing`: reference haplotype file
for setting absolute parental phases. It has the same format as the input genofile,
except that parental genotypes are phased and offspring genotypes are ignored if exist.

`correctthreshold::AbstractFloat=0.15`: a candidate marker is selected for
parental error correction if the fraction of offspring genotypic error >= correctthreshold.

`refinemap::Bool=false`: if true, refine marker map.

`refineorder::Bool=false`: if true, refine marker mordering, valid only if refinemap=true

`maxwinsize::Integer=50`: max size of sliding windown in map refinning.

`inittemperature::Real=4`: initial temperature of simulated annealing in map refinning.

`coolingrate::Real=0.5`: cooling rate of annealing temperature in map refinning.

`stripdis::Real=20`: a chromosome end in map refinement is removed if it has a distance gap > stripdis
(centiMorgan) and it contains less than 5% markers.

`maxepsilon::Real=0.5`: markers in map refinement are removed it they have error
rates > maxepsilon.

`skeletonsize::Integer=50`: the number of markers in the skeleton map that is used
to re-scale inter-map distances.

`missingstring::AbstractString="NA"`: string code for missing value.

`isplot::Bool=false`: if true, plot condprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`: log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

# Output files

Output file | Description
------------- |----------------
`outstem.log`  |  log file
`outstem_maprefined.csv` |  same as input genofile except that marker map being refined
`outstem_parentphased.csv` |  same as input genofile except that  parents being phased
`outstem_parentphased_corrected.csv` |  exported if there exist detected parental errors
`outstem_polyancestry.csv` | genoprob and others (e.g. valent configurations)
`outstem_genoprob.csv` |  posterior genotype probabilities for all offspring
`outstem_postdoseprob.csv` |  posterior dosage probabilities for all offspring
`outstem_plots` |  a folder contains plots of condprob for all offspring if isplot = true


## Example
```julia-repl
julia> polyOrigin("geno.csv","ped.csv")
```

## Citing PolyOrigin

 If you use PolyOrigin in your analyses and publish your results, please cite the article:

  *Zheng C, Amadeu RR, Munoz PR, and Endelman JB. 2020. Haplotype reconstruction in connected tetraploid F1 populations. doi: https://doi.org/10.1101/2020.12.18.423519.*
