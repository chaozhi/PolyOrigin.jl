```@meta
CurrentModule = PolyOrigin
```

# PolyOrigin.jl

**Haplotype reconstruction in connected polyploid F1 populations**

## Features

- Apply to connected F1 in tetraploid (TODO for higher ploidy levels)
- Apply to SNP array and GBS data
- Robust to dosage errors in SNP array data
- Robust to dosage uncertainties in low coverage GBS data
- Robust to errors in input genetic or physical map

## Installation
From the julia (>v1.5.0, 64-bit) REPL model, type `]` to enter the Pkg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Citing PolyOrigin.jl

If you use PolyOrigin in your analyses and publish your results, please cite the article:

  *Zheng C, Amadeu RR, Munoz PR, and Endelman JB. 2020. Haplotype reconstruction in connected tetraploid F1 populations. doi: https://doi.org/10.1101/2020.12.18.423519.*
