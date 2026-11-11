# __precompile__()

"""
PolyOrigin is a package for haplotype reconstruction in connected polyploid F1
populations. See also: [`polyorigin`](@ref), [`polyorigin!`](@ref).
"""
module PolyOrigin

# TODO: speedup:  change dataprob matrix to vector of vector for each marker, for a given offspring
using LinearAlgebra, SparseArrays, Random, Statistics
using Pkg, DelimitedFiles, GZip, Printf,Serialization, Dates
using DataFrames, CSV, Combinatorics, StatsBase, SpecialFunctions
using DataStructures, Distributions, Distributed
using LightGraphs, GraphRecipes, StatsPlots, Plots
using JLD2, Tar, CodecZlib
using ProgressMeter, ThreadsX
using LRUCache

export

    exportall, 
    
    ErrorSetting, 
    PolyGeno, PolyPed, PolyAncestry, JoinedPolyAncestry, 

    polyorigin,polyorigin!,    
    formpolygeno,readpolygeno, savepolygeno, savegenodata, 
    readpolyancestry,savepolyancestry,        
    plotcondprob, animcondprob, plotmapcomp, 
    readtruegeno!,calaccuracy!

# deactivate plot display 
function __init__()    
    ENV["GKSwstype"]="nul"
    # @info string("MagicBase init: ENV[\"GKSwstype\"]=",ENV["GKSwstype"])
end    

include("base/basis.jl")
include("base/vcfbase.jl")
include("base/localbrent.jl")
include("base/posteriordecode.jl")
include("base/posteriordecode_fraction.jl")
include("base/viterbidecode.jl")
include("base/pedigreeui.jl")
include("base/check_infiles.jl")
# include("base/adjusted_boxplot.jl")
include("polyped.jl")
include("polyancestry.jl")
include("joinedpolyancestry.jl")
include("polygeno.jl")
include("parsegeno.jl")
include("parsevcf.jl")
include("priorprocess.jl")
include("likelihood.jl")
include("infer_bvpair.jl")
include("infer_pphase.jl")
include("infer_distance.jl")
include("infer_distance_skeleton.jl")
include("infer_error.jl")
include("infer_order.jl")
include("infer_prefpair.jl")
include("parentcorrect.jl")
include("polyphase.jl")
include("setabsphase.jl")
include("polymarkerdel.jl")
include("polyrefine.jl")
include("polyreconstruct.jl")
include("polyorigin_fun.jl")
include("polyplot.jl")
include("polyplot_valent.jl")
include("calaccuracy.jl")

end # module
