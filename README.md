# PolyOrigin.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chaozhi.github.io/PolyOrigin.jl/dev)
[![Build Status](https://github.com/chaozhi/PolyOrigin.jl/workflows/CI/badge.svg)](https://github.com/chaozhi/PolyOrigin.jl/actions)
[![Coverage](https://codecov.io/gh/chaozhi/PolyOrigin.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/chaozhi/PolyOrigin.jl)

A package for haplotype reconstruction in connected polyploid F1 populations

See [PolyOriginCmd](https://github.com/chaozhi/PolyOriginCmd) for running PolyOrigin in command line interface

See [PolyOriginR](https://github.com/chaozhi/PolyOriginR) for running PolyOrigin in R


## Features

- Apply to connected F1 in tetraploid (TODO for higher ploidy levels. Limited tests for diploids.)
- Apply to SNP array and GBS data
- Robust to dosage errors in SNP array data
- Robust to dosage uncertainties in low read depth GBS data
- Robust to errors in input genetic or physical map

## Installation

From the julia (>=v1.6.7, 64-bit) REPL mode, run

```
julia>]add https://github.com/chaozhi/PolyOrigin.jl
```
where `]` enters Pkg REPL mode.

To update the package, run

```
julia>]up PolyOrigin
```

## Vignettes

[Haplotype reconstruction in a simulated tetraploid population from SNP array data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simarray/step3_tetraploid_simarray.md)

[Haplotype reconstruction in a simulated tetraploid population from GBS data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simgbs/step2_tetraploid_simgbs.md)

[Haplotype reconstruction in a real 3x3 half-diallel potato population](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_realpotato/tetraploid_realpotato.md)

## Usage

From the julia REPL mode, run
```
julia>using PolyOrigin
```

## Help

List the names (types, functions, etc) exported by the package
```
julia>names(PolyOrigin)
```
To get help on a name, type `?` and the name, e.g. polyOrigin
```
julia>?polyOrigin
```

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
```
julia> polyOrigin("geno.csv","ped.csv")
```

## Citing PolyOrigin

 If you use PolyOrigin in your analyses and publish your results, please cite the article:

  *Zheng C, Amadeu RR, Munoz PR, and Endelman JB. 2021. Haplotype reconstruction in connected tetraploid F1 populations. https://doi.org/10.1093/genetics/iyab106
