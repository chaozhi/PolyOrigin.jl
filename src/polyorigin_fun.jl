
"""
    polyOrigin(genofile, pedfile, delimchar=',', missingstring="NA",
        commentstring="#", keyargs...)

performs parental phasing and ancestral inference from input files
and return polyancestry::PolyAncestry. Only ancestral inference is performed
in the case of phased parents.

# Positional arguments

`genofile::AbstractString`: filename for genotypic data file.

`pedfile::AbstractString`:  filename for pedigree info.

see [`readPolyGeno`](@ref) for the requirements of `genofile` and `pedfile`.

# Keyword arguments

`delimchar::AbstractChar=','`:  text delimiter.

`missingstring::AbstractString="NA"`: string code for missing value.

`commentstring::AbstractString="#"`: rows that begin with commentstring will be ignored.

`isphysmap::Bool=false`: if true, input markermap is physical map, where
marker locations are in unit of base pair(bp).

`recomrate::Real=1`: recombination rate in unit of 1 cM/Mbp (centiMorgan per million base pair).
Valid only if `isphysmap=true`.

see keyargs in polyOrigin!(polygeno::PolyGeno, keyargs...)

# Output files

`outstem.log`: log file saves messages that were printed on console.

`outstem_maprefined.csv`: same as the input genofile, except that input marker map is refined.

`outstem_parentphased.csv`: same as the input genofile, except that parental genotypes are phased.

`outstem_parentphased_corrected.csv`: parented genotypes are further corrected.

`outstem_polyancestry.csv`: saves the returned polyancestry. See [`savePolyAncestry`](@ref).

`outstem_genoprob.csv`: a concise version of the above file, including genetic map,
phased parental genotypes, and posterior origin-genotype probabilities. See [`savegenoprob`](@ref).

`outstem_postdoseprob.csv`: same as the input genofile, except that  parent genotypes
are phased and offspring genotypes are given by the posterior dose probabilities.

# Example
```julia-repl
julia> polyOrigin(genofile,pedfile)
```
"""
function polyOrigin(genofile::AbstractString,pedfile::AbstractString;
    delimchar::AbstractChar=',',
    missingstring::AbstractString="NA",
    commentstring::AbstractString="#",
    isphysmap::Bool=false,
    recomrate::Real=1, # 1 cM per Mbp
    epsilon::Real=0.01,
    seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=false,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,
    maxiter::Integer=30,
    minrun::Integer=3,
    maxrun::Integer=10,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    correctthreshold::AbstractFloat=0.15,
    refinemap::Bool=false,
    refineorder::Bool=false,
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan
    maxepsilon::Real=0.5,
    skeletonsize::Integer=50,
    isplot::Bool=false,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    if logfile === nothing
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
        else
            io=logfile
        end
    end
    msg=string("Julia version: ",VERSION)
    printconsole(io,verbose,msg)
    printpkgst(io,verbose,"PolyOrigin")
    starttime = time()
    printconsole(io,verbose,string(repeat("=",33),"polyOrigin",repeat("=",33)))
    msg = string("PolyOrigin, polyOrigin, logfile=", logfile,
        ", ", Dates.now())
    printconsole(io,verbose,msg)
    polygeno = readPolyGeno(genofile,pedfile,
        delimchar=delimchar,
        missingstring=missingstring,
        commentstring=commentstring,
        isphysmap = isphysmap,
        recomrate = recomrate,
        workdir=workdir,
        verbose=false)
    msg = string("list of input files: \n",
            "genofile = ", genofile, "\n",
            "pedfile = ", pedfile)
    printconsole(io,false,msg)
    msg = string("file input/output options: \n",
        "delimchar = ", delimchar, "\n",
        "missingstring = ", missingstring, "\n",
        "commentstring = ", commentstring)
    printconsole(io,false,msg)
    polyancestry=polyOrigin!(polygeno,
        epsilon=epsilon,
        seqerr=seqerr,
        chrpairing_phase=chrpairing_phase,
        chrpairing=chrpairing,
        chrsubset=chrsubset,
        snpsubset=snpsubset,
        isparallel=isparallel,
        delmarker=delmarker,
        delsiglevel=delsiglevel,
        maxstuck=maxstuck,
        maxiter=maxiter,
        minrun=minrun,
        maxrun=maxrun,
        byparent=byparent,
        byneighbor=byneighbor,
        refhapfile=refhapfile,
        refinemap=refinemap,
        refineorder=refineorder,
        maxwinsize=maxwinsize,
        inittemperature=inittemperature,
        coolingrate=coolingrate,
        stripdis=stripdis,
        maxepsilon=maxepsilon,
        skeletonsize=skeletonsize,
        isplot=isplot,
        outstem = outstem,
        correctthreshold=correctthreshold,
        logfile=io,
        workdir= workdir,
        verbose=verbose)
    printconsole(io,verbose,string("End, ", Dates.now(),", total time used = ",
        round(time()-starttime), " seconds by polyOrigin"))
    printconsole(io,verbose,repeat("=",76))
    if typeof(logfile) <: AbstractString
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polyancestry
end

"""
    polyOrigin!(polygeno::PolyGeno, keyargs...)

performs parental phasing and ancestral inference from polygeno
and return polyancestry::PolyAncestry. Only ancestral inference is performed
in the case of phased parents. See [`readPolyGeno`](@ref) for creating
polygeno from inputfiles.

# Positional arguments

`polygeno::PolyGeno`:  a struct storing genotypic data and pedigree info.

# Keyword arguments

`epsilon::Real=0.01`: genotyping error probability.

`seqerr::Real=0.001`: sequencing read error probability for GBS data.

`chrpairing_phase::Integer=22`: chromosome pairing in parental phasing, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome,
with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=false`: if true, multicore computing over chromosomes.

`delmarker::Bool=true`: if true, delete markers during parental phasing.

`delsiglevel::Real=0.05`: significance level for deleting markers.

`maxstuck::Integer=5`: the max number of consecutive iterations that are rejected
in a phasing run.

`maxiter::Integer=30`: the max number of iterations in a phasing run.

`minrun::Integer=3`: if the number of phasing runs having the same parental phases
reaches minrun, phasing algorithm will stop before reaching the maxrun.

`maxrun::Integer=10`: the max number of phasing runs.

`byparent::Union{Nothing,Bool}=nothing`: if true, update parental phases
parent by parent; if false, update parental phases one subpopulation by subpopulation.
The nothing denotes that it is true if a connected component is a single F1 cross,
and false otherwise.

`byneighbor::Union{Nothing,Bool}=nothing`: if ture, udpate the combination of bivalent
or multivalents in parents by their neighbors; if false, consider all the possible combinations.
The nothing denotes that it is true if max ploidy>=6, and false otherwise.

`refhapfile::Union{Nothing,AbstractString} = nothing`: reference haplotype file
for setting absolute parental phases. It has the same format as the input genofile,
except that parental genotypes are phased and offspring genotypes are ignored if exist.

`correctthreshold::AbstractFloat=0.15`: a candidate marker is selected for
parental error correction if the fraction of offspring genotypic error >= correctthreshold.

`refinemap::Bool=false`: if true, refine marker map.

`refineorder::Bool=false`: if true, refine marker mordering, valid only if refinemap=true

`maxwinsize::Integer=50`: max size of sliding windown in map refinning.

`inittemperature::Real=4`: initial temperature of simulated annealing in map refinning.

`coolingrate::Real=0.5`: cooling rate of annealing temperature in map refinning.

`stripdis::Real=20`: a chromosome end in map refinement is removed if it has a distance gap > stripdis
(centiMorgan) and it contains less than 5% markers.

`maxepsilon::Real=0.5`: markers in map refinement are removed it they have error
rates > maxepsilon.

`skeletonsize::Integer=50`: the number of markers in the skeleton map that is used
to re-scale inter-map distances.

`missingstring::AbstractString="NA"`: string code for missing value.

`isplot::Bool=false`: if true, plot condprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

# Example
```julia-repl
julia> polygeno = readPolyGeno(genofile,pedfile)
julia> polyOrigin!(polygeno)
```
"""
function polyOrigin!(polygeno::PolyGeno;
    epsilon::Real=0.01,
    seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=false,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,maxiter::Integer=30,
    minrun::Integer=3,maxrun::Integer=10,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    correctthreshold::AbstractFloat=0.15,
    refinemap::Bool=false,
    refineorder::Bool=false,
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan
    maxepsilon::Real=0.5,
    skeletonsize::Integer=50,
    missingstring::AbstractString="NA",
    isplot::Bool=false,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    if logfile === nothing
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
            printpkgst(io,verbose,"PolyOrigin")
            starttime = time()
            printconsole(io,verbose,string(repeat("=",33),"polyOrigin",repeat("=",33)))
            printconsole(io,verbose,string("PolyOrigin, polyOrigin, logfile=", logfile, ", ", Dates.now()))
        else
            io=logfile
        end
    end
    gdesign= PolyOrigin.plotdesign(polygeno)
    verbose && display(gdesign)
    if isplot && (!isnothing(outstem))
        figdir = joinpath(workdir, outstem * "_plots")
        isdir(figdir) || mkdir(figdir)
        fn = joinpath(figdir,outstem * "_pedigree.png")
        savefig(gdesign,fn)
        msg = string("save design plot in ", fn)
        printconsole(io,verbose,msg)
    end
    if kindofgeno(polygeno.parentgeno) == "phasedgeno"
        phasedgeno = getsubPolyGeno!(polygeno,chrsubset=chrsubset,snpsubset=snpsubset)
    else
        printconsole(io,verbose,string(repeat("-",34),"phasing",repeat("-",35)))
        phasedgeno = polyPhase!(polygeno,
            epsilon=epsilon, seqerr=seqerr,
            chrpairing_phase=chrpairing_phase,
            chrsubset=chrsubset,
            snpsubset=snpsubset,
            isparallel=isparallel,
            delmarker=delmarker,
            delsiglevel=delsiglevel,
            maxstuck=maxstuck,maxiter=maxiter,
            minrun=minrun,maxrun=maxrun,
            byparent=byparent,byneighbor=byneighbor,
            refhapfile=refhapfile,
            missingstring=missingstring,
            outstem=outstem,
            logfile=io,
            workdir=workdir,
            verbose=verbose
        )
    end
    if refinemap
        printconsole(io,verbose,string(repeat("-",32),"maprefinning",repeat("-",32)))
        if correctthreshold<1.0
            polyReconstruct!(phasedgeno,
                epsilon=epsilon, seqerr=seqerr,
                chrpairing=chrpairing,
                chrsubset=nothing,
                snpsubset=nothing,
                isparallel=isparallel,
                correctthreshold=correctthreshold,
                missingstring=missingstring,
                outstem=nothing,
                logfile=io,
                workdir=workdir,
                verbose=verbose
            )
        end
        inputmap = deepcopy(phasedgeno.markermap)
        polyMapRefine!(phasedgeno,
            epsilon=epsilon,
            chrpairing=chrpairing,
            chrsubset=nothing,
            snpsubset=nothing,
            isparallel=isparallel,
            refineorder=refineorder,
            maxwinsize=maxwinsize,
            inittemperature=inittemperature,
            coolingrate=coolingrate,
            stripdis=stripdis,
            maxepsilon=maxepsilon,
            skeletonsize=skeletonsize,
            outstem=outstem,
            logfile=io,
            workdir=workdir,
            verbose=verbose
        )
        if isplot && (!isnothing(outstem))
            figdir = joinpath(workdir, outstem * "_plots")
            isdir(figdir) || mkdir(figdir)
            fn = joinpath(figdir,outstem * "_mapcomp.png")
            plotMapComp(
                inputmap,
                phasedgeno.markermap,
                xlabel = "Input map position (cM)",
                ylabel = "Estimated map position (cM)",
            )
            savefig(fn)
            msg = string("mapcomp file: ", joinpath(splitpath(fn)[end-1:end]...))
            printconsole(io,verbose,msg)
        end
    end
    printconsole(io,verbose,string(repeat("-",31),"reconstructing",repeat("-",31)))
    polyancestry = polyReconstruct!(phasedgeno,
        epsilon=epsilon, seqerr=seqerr,
        chrpairing=chrpairing,
        chrsubset=nothing,
        snpsubset=nothing,
        isparallel=isparallel,
        correctthreshold=correctthreshold,
        missingstring=missingstring,
        isplot=isplot,
        outstem=outstem,
        logfile=io,
        workdir=workdir,
        verbose=verbose
    )
    printconsole(io,verbose,repeat("-",76))
    if refhapfile != nothing && size(polyancestry.correction,1)>0
        try
            setAbsPhase!(refhapfile,phasedgeno,workdir=workdir,io=io,verbose=verbose)
        catch err
            msg = string("warning: could not set absolute parental phases,", err)
            verbose && @warn msg
            printconsole(io,false,msg)
        end
    end
    if typeof(logfile) <: AbstractString
        printconsole(io,verbose,string("End, ", Dates.now(),", total time used = ",
            round(time()-starttime), " seconds by polyOrigin"))
        printconsole(io,verbose,repeat("=",76))
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polyancestry
end
