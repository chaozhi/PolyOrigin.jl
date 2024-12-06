
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
    doseerr::Real=0.01,
    seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    chrpairing_refine::Integer=22,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    isparalleloffspring::Bool=true, 
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
    isinfererror::Bool=false,     
    refinemap::Bool=false,
    refineorder::Bool=false,
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan
    maxdoseerr::Real=0.5,
    skeletonsize::Integer=50,    
    isplot::Bool=false,
    nplot_subpop::Integer=10, 
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    if isnothing(logfile)
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
    printconsole(io,verbose,string(repeat("=",27),"polyOrigin",repeat("=",28)))
    msg = string("polyOrigin, logfile=", logfile,
        ", ", Dates.now())
    printconsole(io,verbose,msg)
    polygeno = readPolyGeno(genofile,pedfile;
        delimchar,missingstring,commentstring,
        isphysmap,recomrate,
        workdir,verbose)
    msg = string("list of file args/options: \n",
            "genofile = ", genofile, "\n",
            "pedfile = ", pedfile, "\n",
            "delimchar = ", delimchar, "\n",
            "missingstring = ", missingstring, "\n",
            "commentstring = ", commentstring)
    printconsole(io,verbose,msg)
    try  
        gdesign= PolyOrigin.plotdesign(polygeno)
        if isplot && !isnothing(outstem)
            figdir = joinpath(workdir, outstem * "_plots")
            isdir(figdir) || mkdir(figdir)        
            fn = joinpath(figdir,outstem * "_pedigree.png")
            savefig(gdesign,fn)
            msg = string("save design plot in ", fn)
            printconsole(io,verbose,msg)
        end
    catch err
        printconsole(io, false, "Failed to plot pedigree")
        @warn err
    end
    polyancestry=polyOrigin!(polygeno;
        doseerr,seqerr,chrpairing_phase,chrpairing_refine,chrpairing,
        chrsubset,snpsubset,isparallel,isparalleloffspring, delmarker,delsiglevel,
        maxstuck,maxiter,minrun,maxrun,byparent,byneighbor,
        refhapfile,refinemap,refineorder,
        maxwinsize,inittemperature,coolingrate,stripdis,
        maxdoseerr,skeletonsize,correctthreshold,
        isinfererror, 
        isplot,nplot_subpop, outstem,logfile=io,workdir,verbose)
    printconsole(io,verbose,string("End, ", Dates.now(),", total time used = ",
        round(time()-starttime), " seconds by polyOrigin"))
    printconsole(io,verbose,repeat("=",66))
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

`doseerr::Real=0.01`: genotyping error probability.

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

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`isparalleloffspring::Bool=true`: if true, multicore computing over offspring in ancestral inference.

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

`refinemap::Bool=false`: if true, refine genetic map. 

`refineorder::Bool=false`: if true, refine marker mordering.

`maxwinsize::Integer=50`: max size of sliding windown in map refinning.

`inittemperature::Real=4`: initial temperature of simulated annealing in map refinning.

`coolingrate::Real=0.5`: cooling rate of annealing temperature in map refinning.

`stripdis::Real=20`: a chromosome end in map refinement is removed if it has a distance gap > stripdis
(centiMorgan) and it contains less than 5% markers.

`maxdoseerr::Real=0.5`: markers in map refinement are removed it they have error
rates > maxdoseerr.

`skeletonsize::Integer=50`: the number of markers in the skeleton map that is used
to re-scale inter-map distances.

`isinfererror::Bool=false`: if true, infer dosage error rate per marker during hapotype reconstruction. 

`missingstring::AbstractString="NA"`: string code for missing value.

`isplot::Bool=false`: if true, plot condprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log"))`:
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
    doseerr::Real=0.01,
    seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    chrpairing_refine::Integer=22,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    isparalleloffspring::Bool=true,
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
    maxdoseerr::Real=0.5,
    skeletonsize::Integer=50,
    isinfererror::Bool=false, 
    isplot::Bool=false,
    nplot_subpop::Integer=10, 
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    if isnothing(logfile)
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
            printpkgst(io,verbose,"PolyOrigin")
            starttime = time()
            printconsole(io,verbose,string(repeat("=",27),"polyOrigin",repeat("=",28)))
            printconsole(io,verbose,string("polyOrigin, logfile=", logfile, ", ", Dates.now()))
        else
            io=logfile
        end
    end    
    isparallel,isparalleloffspring= reset_parallel(isparallel,isparalleloffspring, 
        isinfererror,length(polygeno.markermap); io, verbose)
    if kindofgeno(polygeno.parentgeno) != "phasedgeno" || refinemap
        phase_refine!(polygeno;
            doseerr,seqerr,chrpairing_phase,chrpairing_refine,chrpairing,
            chrsubset,snpsubset,isparallel,delmarker,delsiglevel,
            maxstuck,maxiter,minrun,maxrun,byparent,byneighbor,
            refhapfile,refinemap,refineorder,
            maxwinsize,inittemperature,coolingrate,stripdis,
            maxdoseerr,skeletonsize,correctthreshold,isplot,
            outstem,logfile=io,workdir,verbose)
    end
    printconsole(io,verbose,string(repeat("-",26),"reconstructing",repeat("-",26)))
    polyancestry = polyReconstruct!(polygeno;
        doseerr, seqerr, chrpairing,
        chrsubset=nothing, snpsubset=nothing, 
        isparallel, isparalleloffspring, 
        correctthreshold, byneighbor, isinfererror, 
        missingstring, isplot, nplot_subpop, 
        outstem, logfile=io, workdir, verbose
    )
    printconsole(io,verbose,repeat("-",66))
    if !isnothing(refhapfile) && size(polyancestry.correction,1)>0
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
        printconsole(io,verbose,repeat("=",66))
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polyancestry
end

function reset_parallel(isparallel,isparalleloffspring, isinfererror,nchr; io, verbose)
    isparallel = isparallel && nprocs()>1 && nchr > 1
    isparalleloffspring = isparalleloffspring && nprocs()>1   
    if isparalleloffspring && isparallel
        if isinfererror 
            isparalleloffspring = false # error estimation cannot be paralleled           
            msg = string("reset isparalleloffspring=",isparalleloffspring, " since isparallel=",isparallel,
                " and isinfererror=", isinfererror)                
            printconsole(io, false, "WARN: "*msg)       
            @warn msg
        else
            isparallel = false       
            msg = string("reset isparallel=",isparallel, " since isparalleloffspring=",isparalleloffspring,
                " and isinfererror=", isinfererror)      
            printconsole(io, verbose, msg)                   
        end
    end
    isparallel,isparalleloffspring
end

function phase_refine!(polygeno::PolyGeno;
    doseerr::Real=0.01,
    seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    chrpairing_refine::Integer=22,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
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
    maxdoseerr::Real=0.5,
    skeletonsize::Integer=50,
    missingstring::AbstractString="NA",
    isplot::Bool=false,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    if logfile === nothing
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
            printpkgst(io,verbose,"PolyOrigin")
            starttime = time()
            printconsole(io,verbose,string(repeat("=",27),"phase_refine!",repeat("=",27)))
            printconsole(io,verbose,string("phase_refine!, logfile=", logfile, ", ", Dates.now()))
        else
            io=logfile
        end
    end
    if kindofgeno(polygeno.parentgeno) == "phasedgeno"
        getsubPolyGeno!(polygeno,chrsubset=chrsubset,snpsubset=snpsubset)
    else
        printconsole(io,verbose,string(repeat("-",29),"phasing",repeat("-",30)))
        polyPhase!(polygeno;
            doseerr, seqerr, chrpairing_phase, chrsubset, snpsubset,
            isparallel, delmarker, delsiglevel, maxstuck, maxiter,
            minrun, maxrun, byparent, byneighbor, refhapfile,
            missingstring, outstem, logfile=io, workdir, verbose
        )
    end
    if refinemap
        printconsole(io,verbose,string(repeat("-",28),"refinning",repeat("-",29)))
        if correctthreshold<1.0
            polyReconstruct!(polygeno;
                doseerr, seqerr, chrpairing=chrpairing_refine, 
                isinfererror = true, 
                chrsubset=nothing, snpsubset=nothing,
                isparallel, isparalleloffspring = false, 
                correctthreshold, byneighbor,missingstring,
                outstem=nothing, isplot=false, logfile=io, workdir, verbose
            )
        end        
        polyMapRefine!(polygeno;
            doseerr, chrpairing, chrsubset=nothing, snpsubset=nothing,
            isparallel, byneighbor, refineorder, 
            maxwinsize, inittemperature,
            coolingrate, stripdis, maxdoseerr, skeletonsize,
            outstem, isplot, logfile=io, workdir, verbose
        )
    end
    if typeof(logfile) <: AbstractString
        printconsole(io,verbose,string("End, ", Dates.now(),", total time used = ",
            round(time()-starttime), " seconds by phase_refine!"))
        printconsole(io,verbose,repeat("=",66))
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polygeno
end
