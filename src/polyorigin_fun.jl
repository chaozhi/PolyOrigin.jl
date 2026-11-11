
"""
    polyorigin(genofile, pedfile, delimchar=',', missingstring="NA",
        commentstring="#", keyargs...)

performs parental phasing and ancestral inference from input files. 

# Positional arguments

`genofile::AbstractString`: filename for genotypic data file.

`pedfile::AbstractString`:  filename for pedigree info.

see [`formpolygeno`](@ref) for the requirements of `genofile` and `pedfile`.

# Keyword arguments

`delimchar::AbstractChar=','`:  text delimiter.

`missingstring::AbstractString="NA"`: string code for missing value.

`commentstring::AbstractString="##"`: rows that begin with commentstring will be ignored.

`isphysmap::Bool=false`: if true, input markermap is physical map, where
marker locations are in unit of base pair(bp).

`recomrate::Real=1`: recombination rate in unit of 1 cM/Mbp (centiMorgan per million base pair).
Valid only if `isphysmap=true`.

see keyargs in polyorigin!(polygeno::PolyGeno, keyargs...)

# Output files

`outstem.log`: log file saves messages that were printed on console.

`outstem_maprefined.vcf.gz`: same as the input genofile, except that parents are phased and marker map is refined.

`outstem_imputed.csv`: same as the input genofile, except that  parents are phased, marker map is refined, 
and offspring are imputed. 

`outstem_phased.vcf.gz`: same as the input genofile, except that parents and offspring are phased and marker map is refined.

`outstem_polyancestry.csv`: saves the returned polyancestry. See [`savepolyancestry`](@ref).

`outstem_joinedancestry.csv`: saves polyancestry with joined state space.  

`outstem_plot.tar`: saves all figures. 

`outstem_chrpairing.csv`: summary of the inference results on chromosomal pairing. 

# Example
```julia-repl
julia> polyorigin(genofile,pedfile)
```
"""
function polyorigin(genofile::AbstractString,pedfile::AbstractString;
    formatpriority::AbstractVector=["AD","GP", "GT"],
    delimchar::AbstractChar=',',
    missingstring::AbstractString="NA",
    commentstring::AbstractString="##",
    isphysmap::Bool=false,
    recomrate::Real=1, # 1 cM per Mbp        
    chrpairing_phase::Union{Nothing,Integer}=nothing,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    markersubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,        
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
    isrefine::Bool=true,
    isrefineorder::Bool=false,
    isrefinedis::Union{Nothing,Bool}=nothing,
    isinfererror::Union{Nothing,Bool}=nothing,     
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
    isimputefirst::Bool = true, 
    isparallel::Union{Nothing,Bool}=nothing,     
    isplot::Bool=true,
    nplot_subpop::Integer=1, 
    samplesize::Integer=1000,
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString = ".vcf.gz", 
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)    
    starttime = time()    
    io = set_logfile_begin(logfile, workdir, "polyorigin"; verbose,delim="≡")    
    check_infiles(genofile,pedfile; commentstring, workdir,io, verbose)
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        printconsole(io,verbose,msg)
        printpkgst(io,verbose,"PolyOrigin")
    end    
    polygeno = formpolygeno(genofile,pedfile;
        formatpriority, delimchar,
        missingstring = unique([missingstring, "missing"]),
        commentstring,
        isphysmap,recomrate, workdir)
    msg = string("list of options: \n",
            "genofile = ", genofile, "\n",
            "pedfile = ", pedfile, "\n",
            "isphysmap = ", isphysmap, "\n",
            "recomrate = ", recomrate, "\n",
            "delimchar = ", delimchar, "\n",
            "missingstring = ", missingstring, "\n",
            "commentstring = ", commentstring)
    printconsole(io,verbose,msg)    
    phasedpolygeno, polyancestry, joinedpolyancestry =polyorigin!(polygeno;
        chrpairing_phase,chrpairing,chrsubset,markersubset, bindis,
        phaseerror_midparent,iscorrectparent,  isparentdoseerror,  
        isparallel,isdelvuong,delsiglevel,
        maxstuck,maxiter,minrun,maxrun, byparent, byneighbor,
        refhapfile,trueorigfile, isrefine,isrefineorder,isrefinedis, 
        maxwinsize,inittemperature,coolingrate,stripdis,skeletonsize,
        isinfererror, doseerror,baseerror,allelicbias,overdispersion, 
        tukeyfence, isgreedy, ismarkerspecific,samplesize, 
        isimputefirst, isplot, nplot_subpop, outstem, outext, logfile=io,workdir,verbose
    )    
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, io, tused,"polyorigin"; verbose,delim="≡")
    phasedpolygeno, polyancestry, joinedpolyancestry
end

"""
    polyorigin!(polygeno::PolyGeno, keyargs...)

performs parental phasing and ancestral inference from polygeno
and return polyancestry::PolyAncestry. Only ancestral inference is performed
in the case of phased parents. See [`formpolygeno`](@ref) for creating
polygeno from inputfiles.

# Positional arguments

`polygeno::PolyGeno`:  a struct storing genotypic data and pedigree info.

# Keyword arguments

`chrpairing_phase::Union{Nothing,Integer}=nothing`: chromosome pairing in parental phasing, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome,
with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`markersubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`doseerror::ErrorSetting=ErrorSetting(name="doseerror")`: setting for the inference of dosserror. 
The doseerror.init is used before the inference. 

`baseerror::ErrorSetting=ErrorSetting(name="baseerror")`: setting for the inference of baseerror. 
The baseerror.init is used before the inference. 

`allelicbias::ErrorSetting=ErrorSetting(name="allelicbias")`: setting for the inference of allelicbias. 
The allelicbias.init is used before the inference. 

`overdispersion::ErrorSetting=ErrorSetting(name="overdispersion")`: setting for the inference of overdispersion.

`isparallel::Union{Nothing,Bool}=nothing`: if true, multicore computing over chromosomes.

`isdelvuong::Bool=false`: if true, delete markers during parental phasing.

`delsiglevel::Real=0.01`: significance level for deleting markers.

`maxstuck::Integer=3`: the max number of consecutive iterations that are rejected
in a phasing run.

`maxiter::Integer=50`: the max number of iterations in a phasing run.

`minrun::Integer=3`: if the number of phasing runs having the same parental phases
reaches minrun, phasing algorithm will stop before reaching the maxrun.

`maxrun::Union{Nothing,Integer}=minrun`: the max number of phasing runs.

`byneighbor::Union{Nothing,Bool}=nothing`: if ture, udpate the combination of bivalent
or multivalents in parents by their neighbors; if false, consider all the possible combinations.
The nothing denotes that it is true if max ploidy>=6, and false otherwise.

`refhapfile::Union{Nothing,AbstractString} = nothing`: reference haplotype file
for setting absolute parental phases. It has the same format as the input genofile,
except that parental genotypes are phased and offspring genotypes are ignored if exist.

'trueorigfile::Union{Nothing,AbstractString} = nothing': true ancestral origin file produced by PolyoriginSim. 

`isrefine::Bool=true`: if true, refine marker map.

`isrefineorder::Bool=false`: if true, refine marker ordering, valid only if isrefine=true

`isrefinedis::Union{Nothing,Bool}=nothing`: if true, refine marker distances, valid only if isrefine=true
If nothing, it is set to true only if population size >= 50. 

`maxwinsize::Integer=50`: max size of sliding windown in map refinning.

`inittemperature::Real=4`: initial temperature of simulated annealing in map refinning.

`coolingrate::Real=0.5`: cooling rate of annealing temperature in map refinning.

`stripdis::Real=20`: a chromosome end in map refinement is removed if it has a distance gap > stripdis
(centiMorgan) and it contains less than 5% markers.

`skeletonsize::Union{Nothing,Integer}=100`: the number of markers in the skeleton map that is used
to re-scale inter-map distances.

`missingstring::AbstractString="NA"`: string code for missing value.

`isplot::Bool=true`: if true, plot condprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

# Example
```julia-repl
julia> polygeno = formpolygeno(genofile,pedfile)
julia> polyorigin!(polygeno)
```
"""
function polyorigin!(polygeno::PolyGeno;    
    chrpairing_phase::Union{Nothing,Integer}=nothing,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    markersubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,        
    byparent::Union{Nothing,Bool}=nothing,            
    byneighbor::Union{Nothing,Bool}=nothing,
    bindis::Real=1e-3,  # centiMorgan        
    maxstuck::Integer=3,
    maxiter::Integer=50,    
    minrun::Integer=3,
    maxrun::Union{Nothing,Integer}=minrun,
    refhapfile::Union{Nothing,AbstractString} = nothing,    
    trueorigfile::Union{Nothing,AbstractString} = nothing,    
    iscorrectparent::Union{Nothing,Bool}= nothing, 
    isrefine::Bool=true,
    isrefineorder::Bool=false,
    isrefinedis::Union{Nothing,Bool}=nothing,
    isinfererror::Union{Nothing,Bool}=nothing,    
    ismarkerspecific::Union{Nothing,Bool}=nothing,
    isparentdoseerror::Union{Nothing,Bool}=nothing, 
    phaseerror_midparent::Union{Nothing,Real} = 0.01,     
    doseerror::ErrorSetting=ErrorSetting(name="doseerror"),     
    allelicbias::ErrorSetting=ErrorSetting(name="allelicbias"), 
    overdispersion::ErrorSetting=ErrorSetting(name="overdispersion"), 
    baseerror::ErrorSetting=ErrorSetting(name="baseerror"),     
    tukeyfence::Union{Nothing,Real}=nothing,     
    isgreedy::Bool=false,     
    isdelvuong::Bool=false,
    delsiglevel::Real=0.01,
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,    
    stripdis::Real=20, # centiMorgan    
    skeletonsize::Union{Nothing,Integer}=100,
    missingstring::AbstractString="NA",
    isimputefirst::Bool = true, 
    isparallel::Union{Nothing,Bool}=nothing,     
    isplot::Bool=true,
    nplot_subpop::Integer=5, 
    samplesize::Integer=1000,
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString = ".vcf.gz", 
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    starttime = time()
    io = set_logfile_begin(logfile, workdir, "polyorigin!"; verbose,delim="=")
    isa(logfile, AbstractString) && printpkgst(io,verbose,"PolyOrigin")             
    doseerror,baseerror, allelicbias,overdispersion = reset_error_bais(doseerror,baseerror, allelicbias,overdispersion; io, verbose)    
    # length(parentformatset) > 1 && @error string("TODO for mixed parent genoformat: ", parentformatset)
    # printconsole(io, verbose, string("parent genoformat: ", join(parentformatset,", ")))
    # offformatset = get_genoformat(polygeno; isparentformat=false)
    # length(offformatset) > 1 && @error string("TODO for mixed offspring genoformat: ", offformatset)
    # printconsole(io, verbose, string("offspring genoformat: ", join(offformatset,", ")))
    # priorspace = get_priorstatespace(polygeno.polyped,chrpairing_phase)    
    # isnothing(io) || CSV.write(io, polygeno.polyped.parentinfo; append=true,header=true)            
    # if maxploidy > 2 && isnothing(chrpairing_phase)
    #     popsum = get_popsum(priorspace)
    #     isnothing(io) || CSV.write(io, popsum; append=true,header=true)
    #     if verbose
    #         show(popsum, show_row_number=true, eltypes=false, summary=false)
    #         println()
    #     end            
    # end
    if verbose 
        parentinfo = copy(polygeno.polyped.parentinfo)
        rename!(parentinfo, "individual"=>"parent")    
        show(parentinfo, show_row_number=true, eltypes=false, summary=false)
        println()
    end
    if allelicbias.init != 1.0
        msg = string("allelicbias=", allelicbias, "; allelicbias.init ≠ 1.0 (no bias)")
        @warn msg
        printconsole(io, false, "WARN: "*msg)
    end
    if isplot && (!isnothing(outstem))
        ped2 = extract_compactped(polygeno.polyped.ped)
        gdesign =  plotped(ped2)        
        figdir = joinpath(workdir, outstem * "_plots")
        isdir(figdir) || mkdir(figdir)
        fn = joinpath(figdir,outstem * "_pedigree.png")
        savefig(gdesign,fn)
        msg = string("save design plot in ", fn)
        printconsole(io,verbose,msg)
    end        
    empty!(lru_genolike_GTunphased)
    empty!(lru_genolike_GTphased)    
    empty!(lru_genolike_AD)
    empty!(lru_site_derivedgeno)         
    phasedgeno = polyphase!(polygeno;
        isparallel = isnothing(isparallel) ? true : isparallel, 
        chrpairing_phase, chrsubset, markersubset,byparent, byneighbor, bindis, 
        phaseerror_midparent, isparentdoseerror, 
        iscorrectparent, isinfererror ,            
        doseerror, baseerror, allelicbias, overdispersion, 
        isdelvuong,  delsiglevel,isgreedy, ismarkerspecific, tukeyfence, 
        maxstuck, maxiter,minrun, maxrun, refhapfile,
        missingstring, isplot, outstem, outext, logfile=io, workdir, verbose
    )    
    if isrefine         
        polyrefine!(phasedgeno;
            isparallel = isnothing(isparallel) ? true : isparallel,             
            chrpairing, chrsubset=nothing, markersubset=nothing, bindis,             
            iscorrectparent, isinfererror,              
            phaseerror_midparent, doseerror,baseerror,allelicbias,overdispersion, tukeyfence,     
            byneighbor, isrefineorder, isrefinedis, isdelvuong,  delsiglevel,            
            maxwinsize, inittemperature,
            coolingrate, stripdis, skeletonsize, 
            isplot, outstem, outext, logfile=io, workdir, verbose
        )
        if !isnothing(refhapfile) 
            try                
                refparenthap = setabsphase!(refhapfile,phasedgeno; workdir)
                phaseerrdf = calphaseerr(refparenthap,phasedgeno)
                savephaseerr(phaseerrdf; io, verbose)
            catch err
                msg = string("warning: could not calculate/save phaseing errors,", err)
                verbose && @warn msg
                printconsole(io,false,msg)                
                # throw(err)
            end
        end
    end        
    polyancestry, joinedancestry = polyreconstruct!(phasedgeno;
        chrpairing, chrsubset=nothing, markersubset=nothing, bindis, 
        trueorigfile, 
        tukeyfence = isnothing(tukeyfence) ? 3 : tukeyfence, 
        iscorrectparent = isnothing(iscorrectparent) ? true : iscorrectparent, 
        doseerror=isnothing(doseerror.init) ? 0.01 : doseerror.init, 
        baseerror=baseerror.init,  samplesize, 
        isimputefirst, isparallel, isparallelind=nothing, missingstring, 
        isplot, nplot_subpop, outstem, outext, logfile=io, workdir, verbose
    )        
    empty!(lru_genolike_GTunphased)
    empty!(lru_genolike_GTphased)    
    empty!(lru_genolike_AD)
    empty!(lru_site_derivedgeno)
    if !isnothing(outstem)         
        if !isnothing(refhapfile) 
            try
                phaseerrdf,phaseerrfile = savephaseerr(refhapfile,polyancestry; workdir, outstem)
                savephaseerr(phaseerrdf; io, verbose)  
                msg = string("phaseerrfile: ",relpath(phaseerrfile, workdir))
                printconsole(io,verbose, msg)
            catch err
                msg = string("warning: could not calculate/save phasing errors,", err)
                verbose && @warn msg
                printconsole(io,false,msg)           
            end
        end
        if !isnothing(trueorigfile) 
            try
                ancestryerrfile = saveancestryerror(trueorigfile, joinedancestry; workdir, outstem)
                msg = string("ancestryerrfile: ",relpath(ancestryerrfile, workdir))
                printconsole(io,verbose, msg)
            catch err
                msg = string("warning: could not calculate/save ancestry errors,", err)
                verbose && @warn msg
                printconsole(io,false,msg)     
            end
        end
        if !isnothing(refhapfile)  && !isnothing(trueorigfile) 
            try
                acc = calaccuracy!(refhapfile,trueorigfile, joinedancestry;workdir,verbose=false)
                verbose && display(acc)
                msg = join([string(k=>v) for (k,v) in acc],"\n")
                printconsole(io,false, msg) 
            catch err
                # throw(err)
                msg = string("warning: could not calculate accuracies,", err)
                verbose && @warn msg
                printconsole(io,false,msg)       
            end
        end    
        if isplot 
            figdir = joinpath(workdir, outstem * "_plots")
            if isdir(figdir) 
                try
                    tarplotfile = abspath(workdir, outstem * "_plots.tar")
                    create_tar(figdir, getabsfile(workdir,tarplotfile))                            
                    msg = string("tarplotfile: ",relpath(tarplotfile, workdir))
                    printconsole(io,verbose,msg)                  
                finally 
                    # rm.(readdir(figdir;join=true);force=true)
                    # rm(figdir; force=true,recursive=true)
                end
            end
        end      
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, io, tused,"polyorigin!"; verbose,delim="=")
    phasedgeno, polyancestry, joinedancestry
end

function consistent_dirac_init(rate::Union{ErrorSetting,ErrorSetting})    
    if isa(rate.prior, Dirac) 
        if isnothing(rate.init)
            false
        else
            mean(rate.prior) ≈ rate.init 
        end
    else 
        true
    end
end

function reset_error_bais(doseerror::ErrorSetting, 
    baseerror::ErrorSetting, allelicbias::ErrorSetting,
    overdispersion::ErrorSetting; 
    io::Union{Nothing, IO}, verbose::Bool)
    errls = [doseerror, baseerror, allelicbias,overdispersion]
    errlabls = ["doseerror","baseerror","allelicbias","overdispersion"]
    for i in eachindex(errls,errlabls)
        if !consistent_dirac_init(errls[i])
            (; prior, init, softthresh)  = errls[i]     
            init = mean(prior)
            if i <= 2           
                errls[i] = ErrorSetting(prior, init, softthresh)
            else
                errls[i] = ErrorSetting(prior, init, softthresh)
            end
            msg = string("inconsistent between init=",init, " and prior=", prior, 
                "\n reset ", errlabls[i], "=", errls[i])
            verbose && @warn msg
            printconsole(io,false,msg)        
        end
    end
    Tuple(errls)
end

function reset_parallel(isparallel,isparallelind, nchr, maxploidy;io, verbose)
    isparallel = isparallel && nprocs()>1  
    isparallelind = isparallelind && nprocs()>1   
    if isparallelind && isparallel
        if maxploidy <= 6 
            isparallel = true
            isparallelind = false
        else
            isparallel = false
            isparallelind = true
        end
        msg = string("reset isparallel=",isparallel, " and isparallelind=",isparallelind)
        printconsole(io, verbose, msg)                   
    end
    isparallel,isparallelind
end
