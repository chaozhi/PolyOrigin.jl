# PolyOrigin.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chaozhi.github.io/PolyOrigin.jl/stable)
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
From the julia (>v1.5.0, 64-bit) REPL model, type `]` to enter the pKg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Examples

Examples of using PolyOrigin.jl

[Haplotype reconstruction in a simulated tetraploid multi-parental population from SNP array data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simarray/step3_tetraploid_simarray.md)

[Haplotype reconstruction in a simulated tetraploid multi-parental population from GBS data](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_simgbs/step2_tetraploid_simgbs.md)

[Haplotype reconstruction in a real 3x3 half-diallel potato population](https://github.com/chaozhi/PolyOrigin_Examples/tree/master/tetraploid_realpotato/tetraploid_realpotato.md)


## Citing PolyOrigin

 If you use PolyOrigin in your analyses and publish your results, please cite the article:

  *Zheng C, Amadeu R, Munoz P, and Endelman J. 2020. Haplotype reconstruction in connected tetraploid F1 populations. Manuscript.*
