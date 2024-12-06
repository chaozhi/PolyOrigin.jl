# __precompile__()

"""
PolyOrigin is a package for haplotype reconstruction in connected polyploid F1
populations. See also: [`polyOrigin`](@ref), [`polyOrigin!`](@ref).
"""
module PolyOrigin

# TODO: speed up change dataprob matrix to vector of vector for each marker, for a given offspring

using LinearAlgebra, SparseArrays, Random, Statistics
using Pkg, DelimitedFiles, Printf,Serialization, Dates
using DataFrames, CSV, Combinatorics, StatsBase
using DataStructures, Distributions, Distributed
using LightGraphs, GraphRecipes, StatsPlots, Plots
using ProgressMeter

export

    exportall, 
    
    # public
    PolyGeno, PolyAncestry,
    readPolyGeno,readPolyAncestry,savePolyAncestry,
    polyOrigin,polyOrigin!,
    plotCondprob, animCondprob,plotMapComp,calvalentfreq,plotvalentfreq,
    readTruegeno!,calAccuracy!

include("basis.jl")
include("localbrent.jl")
include("posteriordecode.jl")
include("posteriordecode_fraction.jl")
include("polygeno.jl")
include("polyancestry.jl")
include("priorprocess.jl")
include("datalikelihood.jl")
include("polyphase.jl")
include("setabsphase.jl")
include("polymarkerdel.jl")
include("polymaprefine.jl")
include("polyreconstruct.jl")
include("polyorigin_fun.jl")
include("parentcorrect.jl")
include("infer_prefpair.jl")
include("polyplot.jl")
include("calaccuracy.jl")

end # module
