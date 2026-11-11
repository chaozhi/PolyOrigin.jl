using Pkg
Pkg.add("ArtifactUtils")
using ArtifactUtils

# shell to create datasets and then upload to github
# tar -czf polyorigin2_datasets.tar.gz -C polyorigin2_datasets .

cd(@__DIR__)

add_artifact!(
    "Artifacts.toml", 
    "polyorigin2_datasets", 
    "https://github.com/chaozhi/PolyOrigin.jl/releases/download/PolyOrigin2-datasets/polyorigin2_datasets.tar.gz",
    force=true
)