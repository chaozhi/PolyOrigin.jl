# Get started

## Installation
From the julia REPL model, type `]` to enter the pKg REPL mode and run

```pkg
add https://github.com/chaozhi/PolyOrigin.jl
```

## Usage
To use ```PolyOrigin.jl``` in your project,

```julia
using PolyOrigin
```

To use ```PolyOrigin.jl``` with multi-core computation at chromosome level

```julia
using Distributed
addprocs(4) # set n = min(#CPU threads, #chromosomes)
@everywhere using PolyOrigin
```
and set ```isparallel=true``` in function [`polyOrigin`](@ref). See `addprocs`
for launching worker processes via the specified cluster manager or on remote
machines via SSH.
