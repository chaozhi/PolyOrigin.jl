"""
performs parental phasing and ancestral inference. 

# Arguments

- `genofile`: filename for genotypic data file.
- `pedfile`:  filename for pedigree info.

# Options

- `--delimchar`:  text delimiter.
- `--missingstring`: string code for missing value.
- `--commentstring`: rows that begin with commentstring will be ignored.
- `--recomrate`: recombination rate in unit of 1 cM/Mbp, valid only if `isphysmap`.

# Flags

- `--isphysmap`: input markermap is physical map, where marker locations are in unit of base pair(bp).
- `--isrefine`: refine genetic map
"""
Comonicon.@main function polyorigincli(genofile::AbstractString,pedfile::AbstractString;
    formatpriority::AbstractVector=["AD","GP", "GT"],
    delimchar::AbstractChar=',',
    missingstring::AbstractString="NA",
    commentstring::AbstractString="##",
    isphysmap::Bool=false,
    recomrate::Real=1, # 1 cM per Mbp        
    chrpairing_phase::Union{Nothing,Integer}=nothing,
    chrpairing::Integer=44,
    chrsubset=nothing,
    snpsubset=nothing,        
    bindis::Real=1e-3, # cM
    isparentdoseerror::Union{Nothing,Bool}=nothing, 
    phaseerror_midparent::Union{Nothing,Real} = 0.01, 
    iscorrectparent::Union{Nothing,Bool}=nothing,     
    isdelvuong::Bool=false,
    delsiglevel::Real=0.01,
    maxstuck::Integer=3,
    maxiter::Integer=50,
    minrun::Integer=3,
    maxrun::Union{Nothing,Integer}=minrun,        
    byparent::Union{Nothing,Bool}=nothing,            
    byneighbor::Union{Nothing,Bool}=nothing,
    refhapfile::Union{Nothing,AbstractString} = nothing,    
    trueorigfile::Union{Nothing,AbstractString} = nothing,    
    isrefine::Bool=false,
    isrefineorder::Bool=false,
    isrefinedis::Union{Nothing,Bool}=nothing,
    isrefineerror::Union{Nothing,Bool}=nothing,     
    ismarkerspecific::Union{Nothing,Bool}=nothing,
    doseerror::ErrorSetting=ErrorSetting(name="doseerror"),     
    allelicbias::ErrorSetting=ErrorSetting(name="allelicbias"), 
    overdispersion::ErrorSetting=ErrorSetting(name="overdispersion"), 
    baseerror::ErrorSetting=ErrorSetting(name="baseerror"), 
    tukeyfence::Union{Nothing,Real}=nothing, 
    isgreedy::Bool=false,    
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan    
    skeletonsize::Union{Nothing,Integer}=100,    
    isimputefirst::Bool = false, 
    isparallel::Union{Nothing,Bool}=nothing,     
    isplot::Bool=false,
    nplot_subpop::Integer=5, 
    samplesize::Integer=1000,
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString = ".vcf.gz", 
    logfile::AbstractString = string(outstem,".log"),
    workdir::AbstractString = pwd(),
    verbose::Bool=false
)    
    polyorigin(genofile,pedfile;
        formatpriority, delimchar, missingstring, commentstring,
        isphysmap, recomrate,
        chrpairing_phase,chrpairing,chrsubset,snpsubset, bindis,
        phaseerror_midparent,iscorrectparent,  isparentdoseerror,  
        isparallel,isdelvuong,delsiglevel,
        maxstuck,maxiter,minrun,maxrun, byparent, byneighbor,
        refhapfile,trueorigfile, isrefine,isrefineorder,isrefinedis, 
        maxwinsize,inittemperature,coolingrate,stripdis,skeletonsize,
        isrefineerror, doseerror,baseerror,allelicbias,overdispersion, 
        tukeyfence, isgreedy, ismarkerspecific,samplesize, 
        isimputefirst, isplot, nplot_subpop, outstem, outext, logfile,workdir,verbose
    )       

end