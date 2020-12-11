# Get started

## Installation
From the julia (>v1.5.0, 64-bit) REPL model, type `]` to enter the Pkg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Usage
To use ```PolyOrigin.jl``` in your project,

```julia
using PolyOrigin
```

To use ```PolyOrigin.jl``` with parallel computation of n workers at chromosome level

```julia
using Distributed
addprocs(n)
@everywhere using PolyOrigin
```
and set ```isparallel=true``` in function [`polyOrigin`](@ref). See `addprocs`
for launching worker processes via the specified cluster manager or on remote
machines via SSH.
