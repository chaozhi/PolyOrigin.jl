# __precompile__()

"""
PolyOrigin is a package for haplotype reconstruction in connected polyploid F1
populations. See also: [`polyOrigin`](@ref), [`polyOrigin!`](@ref).
"""
module PolyOrigin

# add Revise Documenter DocumenterTools
# add DataFrames CSV Combinatorics StatsBase
# add DataStructures Distributions Distributed
# add LightGraphs GraphRecipes Plots
using LinearAlgebra, SparseArrays, Random, Statistics
using Pkg, DelimitedFiles, Printf,Serialization, Dates
using DataFrames, CSV, Combinatorics, StatsBase
using DataStructures, Distributions, Distributed
using LightGraphs, GraphRecipes, StatsPlots,Plots

export
    # private function
    # readDesign,parsemarkermap!,splitindex,parseinputgeno,
    # getstatespace,parseprobcell,calpermhomolog,connected_parents,
    # calpermhomolog,calnphaseerr,gethomologdict,getstateorder,ordergenoprob!,
    # printconsole,calabsolutehap,parentgeno2df,condprob2df,calabsaccuracy0,
    # caldataprobset,getpriorprocess,randbvpair,updateepsilonls!, updatedistance!,
    # delsnpeps!,stripchrend!,getskeleton_seg,setmarkerincl!,setabsphase!,
    # calparentacc,calancestryacc,toindexancestry,ancestrycall,getsubPolyGeno,
    # readparenthaplo,calabsolutehap,gethomologdict,calmapkendall,
    # setabsphase0!,calabsaccuracy0,kindofgeno,chrreconstruct,getdict2group,
    # randinitbvpair,getfhaplo,getderiveddose,getbvpairprop,calmarglogl,
    # caldataprob,caldataprob_prob,
    # calchrdoseprob,getpriorstatespace,

    # public
    # polyPhase,polyPhase!,polyMapRefine!,polyReconstruct!,
    # setAbsPhase!,savegenodata,savegenoprob,readTruegeno!,
    # valentprob2df,statespace2df,stringjoin,readdlm2dict,parsestatespace,groupby,
    # getskeleton_seg,splitindex,gridpartition,
    # setAbsPhase!,readparenthaplo,setabsphase0!,

    # public
    PolyGeno, PolyAncestry,
    readPolyGeno,readPolyAncestry,savePolyAncestry,
    polyOrigin,polyOrigin!,
    plotCondprob, animCondprob,plotMapComp,calvalentfreq,plotvalentfreq,
    readTruegeno!,calAccuracy!

include("basis.jl")
include("localbrent.jl")
include("posteriordecode.jl")
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
include("polyplot.jl")
include("calaccuracy.jl")

# precompile(polyOrigin,(String,String,))

end # module
