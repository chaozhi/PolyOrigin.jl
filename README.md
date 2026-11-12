# PolyOrigin.jl

A package for haplotype reconstruction in multigenerational diploid and polyploid pedigrees.


## Features

- Apply to multi-generation pedigrees, including connected F1 populations
- Apply to SNP array and GBS data
- Account for preferential pairing and quadrivalent formation in polyploid
- Account for allelic bias and over-dispersion in GBS data

## Documentation

- [Stable Documentation](https://chaozhi.github.io/PolyOrigin.jl/stable/)
- [Development Documentation](https://chaozhi.github.io/PolyOrigin.jl/dev/)

## Installation

From the julia (>=v1.12, 64-bit) REPL mode, run

```
julia>]add https://github.com/chaozhi/PolyOrigin.jl
```
where `]` enters Pkg REPL mode.

To update the package, run

```
julia>]up PolyOrigin
```

## Usage

From the julia REPL mode, run
```
julia> using PolyOrigin
julia> polyorigin(genofile,pedfile)
```

From command line interface (CLI), run
```
> polyorigin -g genofile -p pedfile
```

## Help

List the names (types, functions, etc) exported by the package
```
julia>names(PolyOrigin)
```
To get help on a name, type `?` and the name, e.g. polyorigin
```
julia>?polyorigin

julia>?polyorigin!
```

From command line interface (CLI), run
```
>polyorigin --help
```


# Output files

Output file | Description
------------- |----------------
`outstem.log`  |  log file
`outstem_maprefined.vcf.gz` |  same as input genofile except that parents are phased and marker map is refined
`outstem_imputed.vcf.gz`| same as the above genofile except that offspring are imputed
`outstem_phased.vcf.gz` |  same as the above genofile except that offspring are phased
`outstem_polyancestry.csv.gz` | saves polyancestry.
`outstem_joinedpolyancestry.csv.gz` | saves polyancestry with joined state space.
`outstem_plots.tar` |  a folder contains all figures if isplot = true
`outstem_chrpairing.csv` | summary of the inference results on chromosomal pairing. 


## Citing PolyOrigin2

 If you use PolyOrigin2 in your analyses and publish your results, please cite the article:

  * Zheng C, Resende MFR, Munoz PR, and Endelman JB. 2026. Haplotype reconstruction in multi-generational diploid and tetraploid populations. 