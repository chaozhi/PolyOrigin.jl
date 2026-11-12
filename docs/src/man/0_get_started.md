# Get started

### Installation

From the julia (>v1.12.0, 64-bit) REPL model, type `]` to enter the Pkg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

The installation will also build PolyOrigin, which might result in a warning message like: 

    To run polyorigin from anywhere in your terminal, add ~/.julia/bin to your PATH

### Usage

To use `PolyOrigin` in your project,

```julia
using PolyOrigin
```

To use `PolyOrigin` with parallel computation of n processes

```julia
using Distributed
nprocs() < n && addprocs(n-nprocs())
addprocs(n)
@everywhere using PolyOrigin
```
See `addprocs` for launching processes via the specified cluster manager or on remote
machines via SSH.

To use PolyOrigin CLI from a terminal, directly run the following command line
```
>polyorigin -g genofile -p pedfile 
```

if ~/.julia/bin is in your PATH, or replace polyorigin by ~/.julia/bin/polyorigin