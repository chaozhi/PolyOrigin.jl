# __precompile__()

"""
    PolyOrigin 
a package for haplotype reconstruction in connected polyploid F1
populations. See also: [`polyOrigin`](@ref), [`polyOrigin!`](@ref).
"""
module PolyOrigin

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
    plotCondprob, animCondprob,plotMapComp,cal_valentsum,plot_valentsum, 
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
