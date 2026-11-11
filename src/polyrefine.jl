
"""
    polyrefine!(polygeno::PolyGeno, keyargs...)

performs marker map refinning for polygeno with phased parent genotypes.
    Modifies polygeno.markermap into a refined genetic map.

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosomes
to be considered, with nothing denoting all chromosomes.
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

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`isrefinedis::Union{Nothing,Bool}=nothing`: if true, refined intermarker distaince. 
If nothing, it is true only if population size >= 50.

`isrefineorder::Bool=false`: if true, refine marker mordering.

`maxwinsize::Integer=50`: max size of sliding window in map refinning.

`inittemperature::Real=4`: initial temperature of simulated annealing in map refinning.

`coolingrate::Real=0.5`: cooling rate of annealing temperature in map refinning.

`stripdis::Real=20`: a chromosome end in map refinement is removed if it has a distance gap > stripdis
(centiMorgan) and it contains less than 5% markers.

`skeletonsize::Union{Nothing,Integer}=100`: the number of markers in the skeleton map that is used
to re-scale inter-map distances.

`logfile::Union{AbstractString,IOStream}=string(outstem,".log")`: output filenames
or stream for writing log.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyrefine!(polygeno::PolyGeno;    
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    markersubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,    
    byneighbor::Union{Nothing,Bool}=nothing,    
    bindis::Real=1e-3, #cM
    isdelvuong::Bool=false,
    delsiglevel::Real=0.01,    
    isrefineorder::Bool=false,
    isrefinedis::Union{Nothing,Bool}=nothing,     
    isinfererror::Union{Nothing,Bool}=nothing,     
    tukeyfence::Union{Nothing,Real}=nothing, 
    iscorrectparent::Union{Nothing,Bool}=nothing,     
    phaseerror_midparent::Union{Nothing,Real} = 0.01, 
    doseerror::ErrorSetting=ErrorSetting(name="doseerror"),     
    allelicbias::ErrorSetting=ErrorSetting(name="allelicbias"), 
    overdispersion::ErrorSetting=ErrorSetting(name="overdispersion"), 
    baseerror::ErrorSetting=ErrorSetting(name="baseerror"),        
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan    
    skeletonsize::Union{Nothing,Integer}=100,
    # missingstring::AbstractString="NA",
    isplot::Bool=true, 
    isparallel::Bool=true,
    outstem::Union{Nothing,AbstractString}="outstem",    
    outext::AbstractString = ".vcf.gz", 
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true,
    verbosemore::Bool=false)
    starttime = time()    
    io = set_logfile_begin(logfile, workdir, "polyrefine!"; verbose)        
    if isnothing(byneighbor)
        byneighbor = max(polygeno.polyped.parentinfo[!,:ploidy]...) >=6
    end        
    doseerror,baseerror, allelicbias,overdispersion = reset_error_bais(doseerror,baseerror, allelicbias,overdispersion; io, verbose)
    msg = string("list of option values: \n",        
        "chrpairing = ",chrpairing, "\n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "markersubset = ", isnothing(markersubset) ? "all markers" : markersubset,"\n",        
        "byneighbor = ",byneighbor,"\n",         
        "bindis = ", bindis, "cM\n",        
        "isdelvuong = ",isdelvuong,"\n",
        "delsiglevel = ",delsiglevel,"\n",
        "isrefineorder = ", isrefineorder, "\n",
        "isrefinedis = ", isrefinedis, "\n",
        "isinfererror = ", isinfererror, "\n",
        "tukeyfence = ", tukeyfence, "\n",
        "iscorrectparent = ", iscorrectparent, "\n",
        "doseerror = ", doseerror, "\n",
        "baseerror = ", baseerror, "\n",
        "allelicbias = ", allelicbias, "\n",
        "overdispersion = ", overdispersion, "\n",        
        "maxwinsize = ",maxwinsize,"\n",
        "inittemperature = ",inittemperature, "\n",
        "coolingrate = ", coolingrate, "\n",
        "stripdis = ", stripdis, "cM\n",        
        "skeletonsize = ", skeletonsize, "\n",
        "isparallel = ",isparallel,"\n",
        "isplot = ",isplot,"\n",
        # "missingstring = ",missingstring,"\n",
        "outstem = ",outstem,"\n",
        "outext = ",outext,"\n",
        "logfile = ",io,"\n",
        "workdir = ",workdir,"\n",
        "verbose = ",verbose,"\n",
        "verbosemore = ",verbosemore)
    printconsole(io,verbose,msg)
    info_polygeno(polygeno; io, verbose)    
    isnothing(chrsubset) || printconsole(io, verbose, string("chrsubset=",chrsubset))    
    isnothing(markersubset) || printconsole(io, verbose, string("markersubset=",markersubset))    
    getsubpolygeno!(polygeno,chrsubset=chrsubset,markersubset=markersubset)
    inputmap = deepcopy(polygeno.markermap)
    if isparallel && length(inputmap) == 1
        isparallel = false
        printconsole(io, verbose, "reset isparallel = false")
    end
    if isnothing(isrefinedis)
        popsize = size(polygeno.polyped.offspringinfo,1)
        isrefinedis = popsize >= 20
        printconsole(io, verbose, string("reset isrefinedis = ", isrefinedis))
    end
    if isnothing(isinfererror)
        popsize = size(polygeno.polyped.offspringinfo,1)
        isinfererror = popsize >= 20
        printconsole(io, verbose, string("reset isinfererror = ", isinfererror))
    end
    if isnothing(tukeyfence)
        popsize = size(polygeno.polyped.offspringinfo,1)
        tukeyfence = popsize >= 20 ? 3 : Inf
        printconsole(io, verbose, string("reset tukeyfence = ", tukeyfence))
    end
    if isnothing(iscorrectparent) 
        iscorrectparent = any([isrefinedis,isinfererror, isrefineorder])
        printconsole(io, verbose, string("reset iscorrectparent = ", iscorrectparent))
    end
    b = [isrefinedis,isinfererror, isrefineorder, iscorrectparent]
    if !any(b) 
        msg = string("polyrefine is not performed since none of [isrefinedis,isinfererror, isrefineorder, iscorrectparent] is true")
        @warn msg
        printconsole(io, false, "Warn: "*msg) 
        return polygeno
    end
    # maprefinement and error estimation
    tused = @elapsed isdel=maprefine_allchr!(polygeno; chrpairing, bindis, 
        isdelvuong, delsiglevel, 
        iscorrectparent, isinfererror, 
        phaseerror_midparent,doseerror,baseerror,allelicbias,overdispersion, tukeyfence,     
        stripdis,skeletonsize,isrefineorder,isrefinedis,           
        maxwinsize,inittemperature,coolingrate,
        isparallel,byneighbor,io,verbose,verbosemore)
    msg = string("tused=", round(tused,digits=1), "s for all chromosomes")
    printconsole(io,verbose,msg)	        
    nmarker = sum(size.(polygeno.markermap,1))
    delmarkerdf = setparentphase!(polygeno,isdel)
    if size(delmarkerdf,1)>0
        msg = string("delete ", size(delmarkerdf,1), " out of ", nmarker, " markers: ",
            join(sort(string.(delmarkerdf[!,:marker])),", "))
    else
        msg = "no marker deleted"        
    end 
    printconsole(io,verbose,msg)   
    info_polygeno(polygeno; io, verbose)
    if !isnothing(outstem)
        outfile = string(outstem,"_maprefined"*outext)
        savegenodata(outfile,polygeno,workdir=workdir)
        msg = string("refined genofile: ", outfile)
        printconsole(io,verbose,msg)
        if isplot            
            figdir = getabsfile(workdir, outstem * "_plots")
            isdir(figdir) || mkdir(figdir)
            # map comparison
            fn = getabsfile(figdir,outstem * "_polyrefine_mapcomp.png")
            plotmapcomp(inputmap, polygeno.markermap;
                maplabels = ["Input map position (cM)",
                    "Estimated map position (cM)"]
            )
            savefig(fn)
            msg = string("mapcomp file: ", joinpath(splitpath(fn)[end-1:end]...))
            printconsole(io,verbose,msg)

            # plot errors
            markermap = reduce(vcat,polygeno.markermap)
            figerr = plotmarkererror(markermap;tukeyfence)
            if !isnothing(figerr)                
                errfile = getabsfile(figdir,outstem * "_polyrefine_markererror.png")
                try 
                    savefig(figerr,getabsfile(workdir,errfile))
                catch err
                    msg = string("Could not savefig. ",err)
                    @warn msg
                    printconsole(io, verbose, "Warning: "*msg)
                end
            end
        end
    end
    if isrefineorder
        tau = round.(calmapkendall(inputmap,polygeno.markermap),digits=4)
        msg = string("Kendall tau between inputmap and refinedmap = ", tau)
        printconsole(io,verbose,msg)
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, io, tused,"polyrefine!"; verbose)
    polygeno
end
    
function maprefine_allchr!(inputpolygeno::PolyGeno;
    chrpairing::Integer,    
    bindis::Real,
    isdelvuong::Bool=false,
    delsiglevel::Real=0.01,
    iscorrectparent::Bool, 
    isinfererror::Bool, 
    isrefinedis::Bool,          
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::ErrorSetting,
    baseerror::ErrorSetting, 
    allelicbias::ErrorSetting,
    overdispersion::ErrorSetting, 
    tukeyfence::Real, 
    stripdis::Real,skeletonsize::Union{Nothing,Integer},
    isrefineorder::Bool,maxwinsize::Integer,inittemperature::Real,coolingrate::Real,
    isparallel::Bool,byneighbor::Bool,
    io::Union{Nothing,IO},
    verbose::Bool,verbosemore::Bool)
    # exclude outlier
    outlier0=inputpolygeno.polyped.offspringinfo[!,:isoutlier]
    outlier = skipmissing(outlier0)
    isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
    noutlier == 0 && !isempty(outlier) && printconsole(io,verbose,string("no offspring excluded"))
    if noutlier > 0
        offoutid = inputpolygeno.polyped.offspringinfo[outlier0,:individual]
        msg = string("exclude ", noutlier, " outlier offspring: ",join(offoutid,", "))
        printconsole(io,verbose,msg)
        polygeno = deepcopy(inputpolygeno)
        deloutlier!(polygeno)
    else
        polygeno = inputpolygeno
    end
    maxiter=30
    nchr=length(polygeno.markermap)
    res=Vector{Vector}(undef,nchr)
    if isparallel && nprocs()>1
        polygenols=[getsubpolygeno(polygeno,chrsubset=[chr]) for chr=1:nchr]
        resmap = pmap(x->maprefine_chr!(x,1; chrpairing, bindis,
            isdelvuong, delsiglevel, 
            iscorrectparent, isinfererror, 
            phaseerror_midparent,doseerror,baseerror,allelicbias,overdispersion,tukeyfence,     
            byneighbor, stripdis, skeletonsize,isrefineorder, isrefinedis, 
            maxwinsize,maxiter,inittemperature,coolingrate,
            io = nothing, verbose,verbosemore),polygenols)
        @everywhere GC.gc()
        for chr=1:nchr
            res[chr], iobuffer = resmap[chr]
            if isnothing(io)
                close(iobuffer)
            else
                write(io,String(take!(iobuffer)))
                flush(io)
            end
        end
    else
        for chr=1:nchr
            res[chr] = first(maprefine_chr!(polygeno,chr; chrpairing, bindis,
                isdelvuong, delsiglevel, 
                iscorrectparent, isinfererror, 
                phaseerror_midparent,doseerror,baseerror,allelicbias,overdispersion,tukeyfence,     
                byneighbor, stripdis, skeletonsize,isrefineorder,isrefinedis, 
                maxwinsize,maxiter,inittemperature,coolingrate,io,verbose,verbosemore))
        end
    end
    if noutlier > 0
        polygeno = inputpolygeno
    end
    isdel=Vector{Vector}(undef,nchr)
    for chr=1:nchr
        chrparentgeno, snporder, isdel[chr], snpid, snppos, targetls, doseerrorls, baseerrorls, allelicbiasls, overdispersionls = res[chr]
        polygeno.markermap[chr] = polygeno.markermap[chr][snporder, :]
        polygeno.markermap[chr][!,:marker] = snpid
        polygeno.markermap[chr][!,:poscm] = snppos
        polygeno.parentgeno[chr] = Matrix{Any}([ismissing(i) ? i : string.(i) for i in chrparentgeno[snporder,:]])
        polygeno.offspringgeno[chr] = polygeno.offspringgeno[chr][snporder,:]
        in("doseerror",targetls) && (polygeno.markermap[chr][!,:doseerror] = doseerrorls[snporder])
        in("baseerror",targetls) && (polygeno.markermap[chr][!,:baseerror] = baseerrorls[snporder])
        in("allelicbias",targetls) && (polygeno.markermap[chr][!,:allelicbias] = allelicbiasls[snporder])
        in("overdispersion",targetls) && (polygeno.markermap[chr][!,:overdispersion] = overdispersionls[snporder])
    end
    isdel
end



function maprefine_chr!(polygeno::PolyGeno, chr::Integer;     
    chrpairing::Integer,        
    bindis::Real, 
    isdelvuong::Bool=false,
    delsiglevel::Real=0.01,
    iscorrectparent::Bool, 
    isinfererror::Bool,
    isrefinedis::Bool, 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::ErrorSetting,
    baseerror::ErrorSetting,
    allelicbias::ErrorSetting,
    overdispersion::ErrorSetting, 
    tukeyfence::Real, 
    byneighbor::Bool,
    stripdis::Real,
    skeletonsize::Union{Nothing,Integer},
    priorlength::Real=100, 
    isrefineorder::Bool,
    maxwinsize::Integer,
    maxiter::Integer,
    inittemperature::Real,
    coolingrate::Real,    
    io::Union{IO,Nothing},
    verbose::Bool,
    verbosemore::Bool)
    starttime = time()
    io===nothing && (io=IOBuffer(append=true))
    polyped = polygeno.polyped
    maxploidy_off = maximum(polyped.offspringinfo[!,:ploidy])
    priorspace = get_priorstatespace(polyped,chrpairing)
    idlen=max([length(i[1,:linkagegroup]) for i=polygeno.markermap]...)
    chrid=lpad(polygeno.markermap[chr][1,:linkagegroup],idlen)
    chrparentgeno= getchrparentgeno(polygeno.parentgeno[chr])
    chrmultiallele = polygeno.markermap[chr][!,:multiallele]
    derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)
    chroffgeno = get_chroffgeno(polygeno, chr)
    chroffformat = polygeno.markermap[chr][!,:offspringformat]        
    targetls = doseerror.isinfer ?  ["doseerror"] : String[] 
    if !all(occursin.(r"^GT_",chroffformat))         
        baseerror.isinfer && push!(targetls,"baseerror")
        allelicbias.isinfer && push!(targetls,"allelicbias")
        overdispersion.isinfer && push!(targetls,"overdispersion")
    end   
    isinfererror = isinfererror && !isempty(targetls)
    snporder=collect(1:size(chroffgeno,1))    
    ismissdoseerror = all(ismissing.(polygeno.markermap[chr][!,:doseerror]))
    ismissbaseerror = all(ismissing.(polygeno.markermap[chr][!,:baseerror]))
    ismissbias = all(ismissing.(polygeno.markermap[chr][!,:allelicbias]))
    ismissod = all(ismissing.(polygeno.markermap[chr][!,:overdispersion]))    
    
    doseerrinit = isnothing(doseerror.init) ? 0.01 : doseerror.init
    doseerrorls = reset_errorls(polygeno.markermap[chr][!,:doseerror],doseerrinit)
    baseerrorls = reset_errorls(polygeno.markermap[chr][!,:baseerror],baseerror.init)    
    allelicbiasls = reset_errorls(polygeno.markermap[chr][!,:allelicbias],allelicbias.init)    
    overdispersionls = reset_errorls(polygeno.markermap[chr][!,:overdispersion],overdispersion.init)        
    if ismissdoseerror || !in("doseerror",targetls)
        doseerrorprior = doseerror.prior    
    else
        avg = max(0.01,mean(doseerrorls))
        doseerrorprior = Beta(1,round(1/avg - 1,digits=1))         
        msg = string("chr=",chrid,", reset doseerror.prior = ", doseerrorprior)
        printconsole(io,verbose,msg)        
    end
    if ismissbaseerror || !in("baseerror",targetls) 
        baseerrorprior = baseerror.prior    
    else        
        avg = max(0.005,mean(baseerrorls))
        baseerrorprior = Beta(1,round(1/avg - 1,sigdigits=2))         
        msg = string("chr=",chrid,", reset baseerror.prior = ", baseerrorprior)
        printconsole(io,verbose,msg)        
    end    
    if ismissbias || !in("allelicbias",targetls)  
        biasprior = allelicbias.prior    
    else
        priorlogstd = round(std(log.(rand(allelicbias.prior,10000))),digits = 2)                 
        biasstd = round(std(log.(allelicbiasls)),digits = 2)
        if priorlogstd < biasstd
            biasprior = LogNormal(0, max(biasstd,0.1))
            msg = string("chr=",chrid,", reset allelicbias.prior = ", biasprior)
            printconsole(io,verbose,msg)       
        else
            biasprior = allelicbias.prior                 
            msg = string("chr=",chrid,", keep allelicbias.prior = ", biasprior)
            printconsole(io,verbose,msg)       
        end
    end
    if ismissod || !in("overdispersion",targetls) 
        overdispersionprior = overdispersion.prior    
    else        
        avg = max(0.05,round(mean(overdispersionls),digits=2))
        overdispersionprior = Exponential(avg) 
        msg = string("chr=",chrid,", reset overdispersion.prior = ", overdispersionprior)
        printconsole(io,verbose,msg)        
    end
    errorprior = (doseerror=doseerrorprior, baseerror = baseerrorprior, 
            allelicbias = biasprior, overdispersion = overdispersionprior)
    # errorprior = (doseerror=doseerror.prior, baseerror = baseerror.prior, 
    #     allelicbias = allelicbias.prior, overdispersion = overdispersion.prior)
    errorsoftthresh = (doseerror = doseerror.softthresh, baseerror = baseerror.softthresh, 
        allelicbias = allelicbias.softthresh, overdispersion = overdispersion.softthresh)
    errorhardthresh = (doseerror = doseerror.hardthresh, baseerror = baseerror.hardthresh, 
        allelicbias = allelicbias.hardthresh, overdispersion = overdispersion.hardthresh)    

    maxvalent = max(digits(chrpairing)...)
    priorprocess = get_priorprocess(polygeno,chr,maxvalent; bindis)
    bvpair = init_bvpair(priorspace,polyped)    
    nowlogl =refinebvpair!(bvpair, chrparentgeno,chroffgeno,chroffformat, chrmultiallele, 
        priorspace,priorprocess, polyped; 
        phaseerror_midparent,byneighbor, snporder,
        doseerrorls, baseerrorls, allelicbiasls, overdispersionls, 
    )    
    msg_err = get_msg_err(targetls, priorprocess; snporder, 
        doseerror = doseerrorls, baseerror = baseerrorls, 
        allelicbias = allelicbiasls, overdispersion = overdispersionls)
    chrlen = round(get_chrlen(priorprocess),digits=1)                   
    msg = string("chr=",chrid,", it=0", msg_err, ", len=",chrlen, "cM", 
        ", logl=",round(nowlogl,digits=1))
    printconsole(io,verbose,msg)
    
    maxnstuck = isrefineorder ? 3 : 5
    temperature = isrefineorder ? inittemperature : 0.0
    nstuck = 0.0
    nstuck_del = 0    
    nstuck_correct = 0
    ndel_err = 0
    winsize=maxwinsize    
    for it=1:maxiter
        msg =string("chr=",chrid,", it=",it)
        startt = time()
        # update ordering        
        if isrefineorder            
            snporder,priorprocess,logbook =updateorder!(chrparentgeno, chroffgeno, chroffformat, chrmultiallele, 
                derivedgeno, bvpair, priorspace, priorprocess;
                phaseerror_midparent, snporder, 
                doseerrorls, baseerrorls, allelicbiasls,overdispersionls, isreversechr = iseven(it), temperature, 
            )            
            newwinsize = length(logbook)+1
            if temperature ≈ 0.0
                nstuck +=  newwinsize == winsize ? 1.0 : 0.5
            else                
                nstuck = 0
            end
            winsize=newwinsize
            nowlogl = logbook[end][2]
        else
            logbook = []
            nstuck += 1.0 # 5 iterations when refine only distance, since maxnstuck = 5 in this case
        end
        # update doseerror
        if isinfererror            
            loglls_err, ndel_err = infer_errorls!(nothing, nothing, chrparentgeno, derivedgeno, chroffgeno, 
                chroffformat, chrmultiallele, bvpair, priorspace, priorprocess, polyped, targetls;     
                phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, errorprior, 
                errorsoftthresh,tukeyfence, errorhardthresh, bindis, snporder, temperature
            )                                            
            msg_err = get_msg_err(targetls, priorprocess; snporder, 
                doseerror = doseerrorls, baseerror = baseerrorls, 
                allelicbias = allelicbiasls, overdispersion = overdispersionls)
            msg *= string(msg_err, ", #del_err=",ndel_err)
            logl_err = last(loglls_err)
            deltlogl = round(logl_err-nowlogl,digits=2)
            push!(logbook,["#update_error",join(round.(loglls_err,digits=2),"|"),deltlogl,1])
            nowlogl = logl_err                        
            empty!(lru_genolike_GTunphased)
            empty!(lru_genolike_GTphased)    
            empty!(lru_genolike_AD)    
        end 
        if isrefinedis 
            # update inter-marker distances
            logldis = updatedistance!(snporder,bvpair,priorspace,priorprocess,
                chrparentgeno, chroffgeno, chroffformat, chrmultiallele, derivedgeno; 
                phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, 
                bindis, priorlength, temperature,
            )                    
            push!(logbook,["#update_distance",logldis,round(logldis-nowlogl,digits=2),1])    
            nowlogl = logldis    
            ndel3 = stripchrend!(priorprocess,stripdis=stripdis+5*temperature)                                        
            ndel3 > 0 && (msg *= string(", #del_strip=",ndel3))            
        end        
        if isdelvuong && (temperature ≈ 0.0)
            deltt = polymarkerdel!(bvpair, priorspace, priorprocess,
                chrparentgeno, chroffgeno, chroffformat, chrmultiallele; 
                phaseerror_midparent,doseerror=doseerrorls, baseerror=baseerrorls, 
                allelicbias=allelicbiasls, overdispersion = overdispersionls,
                delsiglevel, snporder,bindis, 
            )               
            ndel= length(deltt)                  
            msg *= string(", #del_test=",ndel)
            if ndel == 0  
                nstuck_del += 1
            else
                nstuck_del = 0                
            end
            (nstuck_del >= 2) && (isdelvuong = false)                                   
        end            
        if iscorrectparent && (temperature ≈ 0.0)            
            changelogl, changedf = parentcorrect(chrparentgeno, chroffgeno, chroffformat, chrmultiallele, 
                priorspace, priorprocess, polyped; 
                phaseerror_midparent, doseerror=doseerrorls, baseerror=baseerrorls, 
                allelicbias=allelicbiasls, overdispersion = overdispersionls, 
                snporder,  bvpair, 
                callthreshold = 0.7, minerrdiff = 3, 
            )                                 
            correctapt = changelogl > 0 
            correctapt && correct_chrphalo!(chrparentgeno, changedf)             
            ncorrect = size(changedf,1)
            msg *= string(", #correct=", ncorrect)
            ncorrect > 0 && (msg *= correctapt ? "(accepted)" : "(rejected)")
            if ncorrect == 0 || !correctapt
                nstuck_correct += 1                
            else
                nstuck_correct = 0
                derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)            
            end 
            (nstuck_correct >= 2) && (iscorrectparent = false)                       
        end  
        # update bvpair
        loglbv = refinebvpair!(bvpair, chrparentgeno,chroffgeno,chroffformat, chrmultiallele, 
            priorspace,priorprocess, polyped; 
            phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls,byneighbor, snporder,
            isrand = temperature > 0, 
        )        
        push!(logbook,["#update_valent",loglbv,round(loglbv-nowlogl,digits=2),1])    
        nowlogl = loglbv


        colid = Symbol.([string("#refine_",chrid),"logl","delt_logl","accept"])
        dflogbook = DataFrame(permutedims(reduce(hcat,logbook)),colid)        
        msgmore = isrefineorder ? string(", win=",winsize,", T=",round(temperature,digits=4)) : ""
        chrlen = round(get_chrlen(priorprocess),digits=1)                   
        msg *= string(", len=",chrlen, "cM", 
            ", logl=",round(nowlogl,digits=1),            
             msgmore, 
             ", tused=", round(time()-startt, digits=1), "s"
            )
        printconsole(io,verbose,msg)
        if verbosemore
            @info dflogbook
            CSV.write(io, dflogbook,header=true,append=true)
        end        
        isnochange = ndel_err == 0 && !isdelvuong && !iscorrectparent
        if it==maxiter || (temperature ≈ 0.0 && nstuck >= maxnstuck && isnochange && isodd(it))                                    
            pri1 = first(values(priorprocess))            
            deltdls = pri1.markerdeltd[pri1.markerincl][1:end-1]
            nseg = sum(deltdls .> 0)
            if isrefinedis && (isnothing(skeletonsize) || nseg >skeletonsize)                                
                skeletonprior = infer_skeletonmap(chrparentgeno, chroffgeno, chroffformat, chrmultiallele, 
                    derivedgeno, bvpair, priorspace, priorprocess, polyped,targetls,
                    phaseerror_midparent, isinfererror, 
                    doseerrorls, baseerrorls, allelicbiasls, overdispersionls, errorprior, 
                    errorsoftthresh,tukeyfence,errorhardthresh;
                    startit = it+1, skeletonsize, bindis, snporder, chrid, io, verbose, 
                )                
                rescale_map!(priorprocess,skeletonprior,bindis)    
                ndel3 = stripchrend!(priorprocess,stripdis=stripdis+5*temperature)                                        
                if ndel3 > 0 
                    msg = string("chr=",chrid, ", #del_strip=",ndel3)
                    printconsole(io,verbose,msg)
                end
            else
                msg = string("chr=",chrid, ", #seg=",nseg, ", no distance rescaling by skeleton markers")
                printconsole(io,verbose,msg)
            end            
            break
        end
        if temperature > 0.5 
            temperature *= coolingrate
        elseif temperature > 0.1
            temperature *= max(0.5,coolingrate^2)
        else
            temperature = 0.0
        end
    end
    chrlen = round(get_chrlen(priorprocess),digits=1)                
    nmarker = sum(first(values(priorprocess)).markerincl)
    msg = string("chr=",chrid,        
        ", done, elapsed=",round(time()-starttime), "seconds",
        ", len=",chrlen, "cM", ", #marker=",nmarker)
    printconsole(io,verbose,msg)
    markerid=polygeno.markermap[chr][!,:marker]
    checksnporder(snporder,priorprocess,markerid) || @error("inconsistent snporder and markerid!")
    checktranprobseq(priorprocess,bindis) || @error("inconsistent tranprobseq and deltd!")
    pri1=first(values(priorprocess))
    ismarkerdel = .!(pri1.markerincl)    
    snpid = pri1.markerid
    dls = copy(pri1.markerdeltd)
    dls[isnothing.(dls)] .= 0.0
    snppos = round.(vcat([0],100*accumulate(+,dls[1:end-1])),digits=2)
    # println("largeerr=",sort(doseerrorls[snporder[pri1.markerincl]],rev=true)[1:10])
    [chrparentgeno, snporder, ismarkerdel, snpid, snppos, targetls, doseerrorls, baseerrorls, allelicbiasls,overdispersionls], io
end

function get_chrlen(priorprocess::AbstractDict)
    pri1=first(values(priorprocess))
    100*sum(pri1.markerdeltd[pri1.markerincl][1:end-1])
end

function refinebvpair!(bvpair::NamedTuple, 
    chrparentgeno::AbstractMatrix, chroffgeno::AbstractMatrix, 
    chroffformat::AbstractVector, chrmultiallele::AbstractVector, 
    priorspace::AbstractDict, priorprocess::AbstractDict, polyped::PolyPed;     
    snporder::Union{Nothing,AbstractVector}=nothing,
    phaseerror_midparent::Union{Nothing,Real},
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector, 
    overdispersionls::AbstractVector, 
    byneighbor::Bool=false,
    isrand::Bool=false,
    incldelmarker::Bool=false)    
    phaploset=[[[i] for i=j] for j=eachcol(chrparentgeno)]
    phaploindex=[ones(Int, length(j)) for j=phaploset]    
    siblogl = init_siblogl(bvpair)
    rand_update_bvpair!(siblogl, bvpair,phaploindex,phaploset, priorspace, priorprocess,
        chroffgeno,chroffformat, chrmultiallele, polyped; 
        phaseerror_midparent,doseerror=doseerrorls, baseerror=baseerrorls, 
        allelicbias=allelicbiasls, overdispersion = overdispersionls,
        byneighbor,  isrand, snporder, incldelmarker, 
    )    
    # println("siblogl.parent=",siblogl.parent,", siblogl.offspring=",siblogl.offspring)
    # println("siblogl.parent=",sum(skipmissing(siblogl.parent)),", siblogl.offspring=",sum(skipmissing(siblogl.offspring)))
    # println("bvpair.parent=",bvpair.parent,", bvpair.offspring=",bvpair.offspring)
    sum_margloglls(siblogl)
end


function checksnporder(snporder::AbstractVector,priorprocess::AbstractDict,
    markerid::AbstractVector)
    pri1=first(values(priorprocess))
    markerid[snporder] == pri1.markerid
end

function checktranprobseq(priorprocess::AbstractDict,bindis::Real)
    pri1=first(values(priorprocess))
    incl = pri1.markerincl
    deltd = pri1.markerdeltd[incl][1:end-1]
    transeq=[getgametetran(i,pri1.nvalent,bindis) for i=deltd]
    transeq2 = pri1.tranprobseq[incl][1:end-1]

    issametran = map((x,y)->any(isnothing.([x,y])) || isapprox(x,y,atol=1e-4), transeq, transeq2)
    all(issametran) || @warn string("inconsistent markerdeltd and tranprobseq")

    b = [(isnothing(transeq[i]) || isnothing(transeq2[i])) ? true : transeq[i] ≈ transeq2[i] 
        for i in eachindex(transeq, transeq2) ]
    for i in findall(.!b)
        println("i=",i, ", d=",deltd[i], ",didff=",sum(abs.(transeq[i] .- transeq2[i])))
    end
    all(b)
end




function deloutlier!(polygeno::PolyGeno)
    outlier = polygeno.polyped.offspringinfo[!,:isoutlier]
    any(ismissing.(outlier)) && return polygeno
    if sum(outlier)>0
        incl = .! outlier
        polygeno.offspringgeno=[i[:,incl] for i=polygeno.offspringgeno]
        (;ped, designinfo, parentinfo, offspringinfo) = polygeno.polyped
        polygeno.polyped = PolyPed(ped, designinfo, parentinfo, offspringinfo[incl,:]) 
    end
    polygeno
end

function calmapkendall(inputmap::Vector{DataFrame},refinedmap::Vector{DataFrame})
    [begin
        markers = inputmap[chr][!,:marker]
        inputorder = Dict(markers .=> 1:length(markers))
        markers = refinedmap[chr][!,:marker]
        estorder = collect(skipmissing([get(inputorder,i,missing) for i=markers]))
        corkendall(estorder,collect(1:length(estorder)))
    end for chr=1:length(inputmap)]
end