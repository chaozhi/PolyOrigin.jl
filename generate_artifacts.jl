using Pkg
Pkg.add("ArtifactUtils")
using ArtifactUtils

cd(@__DIR__)

add_artifact!(
    "Artifacts.toml", 
    "polyorigin2_datasets", 
    "https://github.com/chaozhi/PolyOrigin.jl/releases/download/datasets-v2/polyorigin2_datasets.tar.gz",
    force=true
)