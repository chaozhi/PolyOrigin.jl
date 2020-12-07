```@meta
CurrentModule = PolyOrigin
```

# PolyOrigin.jl

**Haplotype reconstruction in connected polyploid F1 populations**

## Installation
From the julia (>v1.5.0) REPL model, type `]` to enter the pKg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Features

- Apply to connected F1 in tetraploid (TODO for higher ploidy levels)
- Apply to SNP array and GBS data
- Robust to dosage errors in SNP array data
- Robust to dosage uncertainties in low coverage GBS data
- Robust to errors in input genetic or physical map

## Citing PolyOrigin.jl

 If you use PolyOrigin.jl in your analyses and publish your results, please cite the article:

  *Zheng C, Amadeu R, Munoz P, and Endelman J. 2020. Haplotype reconstruction in connected tetraploid F1 populations. Manuscript.*
