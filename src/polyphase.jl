function connected_polygenos(polygeno::PolyGeno)
    polyped = polygeno.polyped
    ccparent = connected_parents(polyped)
    length(ccparent) ==1 && return [polygeno]
    ccpop = [findall(sum(Matrix(polyped.designinfo[:,i .+ 1]),dims=2)[:,1] .> 0) for i= ccparent]
    ccpopid=[polyped.designinfo[!,:family][i] for i=ccpop]
    offpop=polyped.offspringinfo[!,:family]
    ccoff=[findall([i in j for i=offpop]) for j=ccpopid]
    [begin
        ppgeno = [i[:,ccparent[c]] for i=polygeno.parentgeno]
        offgeno = [i[:,ccoff[c]] for i=polygeno.offspringgeno]
        ppinfo = polyped.parentinfo[ccparent[c],:]
        offinfo =  polyped.offspringinfo[ccoff[c],:]
        cols = ccparent[c] .+ 1
        pushfirst!(cols,1)
        desinfo = polyped.designinfo[ccpop[c],cols]
        ped = getsubped(polyped.ped, offinfo[:,:individual])        
        newpolyped = PolyPed(ped, desinfo, ppinfo,offinfo)
        PolyGeno(newpolyped, polygeno.markermap, ppgeno,offgeno,polygeno.misc)
    end for c = eachindex(ccparent)]
end

function connected_polypeds(polyped::PolyPed)    
    ccparent = connected_parents(polyped)
    length(ccparent) ==1 && return [polygeno]
    ccpop = [findall(sum(Matrix(polyped.designinfo[:,i .+ 1]),dims=2)[:,1] .> 0) for i= ccparent]
    ccpopid=[polyped.designinfo[!,:family][i] for i=ccpop]
    offpop=polyped.offspringinfo[!,:family]
    ccoff=[findall([i in j for i=offpop]) for j=ccpopid]
    [begin        
        ppinfo = polyped.parentinfo[ccparent[c],:]
        offinfo =  polyped.offspringinfo[ccoff[c],:]
        cols = ccparent[c] .+ 1
        pushfirst!(cols,1)
        desinfo = polyped.designinfo[ccpop[c],cols]
        ped = getsubped(polyped.ped, offinfo[:,:individual])        
        PolyPed(ped, desinfo, ppinfo,offinfo)        
    end for c = eachindex(ccparent)]
end


function connected_parents(polyped::PolyPed)
    parent2off = get_parent2progeny(polyped)
    pls = polyped.parentinfo[!,:individual]
    np = length(pls)
    adj = zeros(Int, np,np)
    for i=1:np
        offls = parent2off[pls[i]]
        for j in 1:np
            p = pls[j]
            if in(p, offls) || !isempty(intersect(offls, parent2off[p]))
                adj[i,j] = 1
            end
        end    
    end
    adj=sign.(adj + adj')
    g = SimpleGraph(adj)
    connected_components(g)
end

# same as get_parent_2_immediateoff(polyped)
function get_parent2progeny(polyped::PolyPed)    
    priorspace = get_priorstatespace(polyped,22)
    pls = polyped.parentinfo[!,:individual]
    parent2off = Dict([i=>[] for i in pls])
    for (key, val) in priorspace
        offls = vcat(val["offspring"], val["midparent"])
        pp = val["parent"]
        for p in pp
            union!(parent2off[p], offls)
        end    
    end
    parent2off
end


"""
    polyphase!(polygeno::PolyGeno, keyargs...)

performs parental phasing and return phasedgeno::PolyGeno
with  phasedgeno.parentgeno and polygeno.parentgeno being phased .

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`doseerror::ErrorSetting=ErrorSetting(name="doseerror")`: setting for the inference of dosserror. 
The doseerror.init is used before the inference. 

`baseerror::ErrorSetting=ErrorSetting(name="baseerror")`: setting for the inference of baseerror. 
The baseerror.init is used before the inference. 

`allelicbias::ErrorSetting=ErrorSetting(name="allelicbias")`: setting for the inference of allelicbias. 
The allelicbias.init is used before the inference. 

`overdispersion::ErrorSetting=ErrorSetting(name="overdispersion")`: setting for the inference of overdispersion.

`chrpairing_phase::Union{Nothing,Integer}=nothing`: chromosome pairing in parental phasing, with 22 being only
bivalent formations and 44 being bi- and quadri-valent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosomes
to be considered, with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`isdelvuong::Bool=false`: if true, delete markers during parental phasing.

"delsiglevel::Real=0.01": significance level for deleting markers.

`maxstuck::Integer=3`: the max number of consecutive iterations that are rejected
in a phasing run.

`maxiter::Integer=50`: the max number of iterations in a phasing run."

`minrun::Integer=3`: if the min number of phasing runs that are at the same local maximimum or
have the same parental phases reaches minrun, phasing algorithm will stop before reaching the maxrun.

`maxrun::Union{Nothing,Integer}=minrun`: the max number of phasing runs.
 
`byneighbor::Union{Nothing,Bool}=nothing`: if ture, udpate the combination of bivalent
or multivalents in parents by their neighbors; if false, consider all the possible combinations.
The nothing denotes thtat it is true if max ploidy>=6, and false otherwise.

`refhapfile::Union{Nothing,AbstractString} = nothing`: reference haplotype file
for setting absolute parental phases. It has the same format as the input genofile,
except that parental genotypes are phased and offspring genotypes are ignored if they exist.

`missingstring::AbstractString="NA"`: string code for missing value.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyphase!(polygeno::PolyGeno;    
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    chrpairing_phase::Union{Nothing,Integer}=nothing,    
    byparent::Union{Nothing,Bool}=nothing,            
    byneighbor::Union{Nothing,Bool}=nothing,     
    bindis::Real=1e-3,     
    iscorrectparent::Union{Nothing,Bool}=nothing,     
    isrefineerror::Union{Nothing,Bool}=nothing,    
    isparentdoseerror::Union{Nothing,Bool}=nothing, 
    phaseerror_midparent::Union{Nothing,Real} = 0.01, 
    tukeyfence::Union{Nothing,Real}=nothing, 
    doseerror::ErrorSetting=ErrorSetting(name="doseerror"),     
    allelicbias::ErrorSetting=ErrorSetting(name="allelicbias"), 
    overdispersion::ErrorSetting=ErrorSetting(name="overdispersion"), 
    baseerror::ErrorSetting=ErrorSetting(name="baseerror"),         
    isgreedy::Bool=false, 
    ismarkerspecific::Union{Nothing,Bool}=nothing, 
    isdelvuong::Bool=false, 
    delsiglevel::Real=0.01,     
    maxstuck::Integer=3,
    maxiter::Integer=50,
    minrun::Integer=3,
    maxrun::Union{Nothing,Integer}=minrun,
    runblocksize::Union{Nothing,Integer}=minrun, 
    refhapfile::Union{Nothing,AbstractString} = nothing,
    missingstring::AbstractString="NA",
    isparallel::Bool=true,
    isplot::Bool=true, 
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString = ".vcf.gz", 
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    tempdirectory::AbstractString = tempdir(),
    verbose::Bool=true)
    starttime = time()    
    io = set_logfile_begin(logfile, workdir, "polyphase!"; verbose)    
    minrun = max(1, minrun)    
    if isparallel && nworkers()  == 1
        isparallel = false
        printconsole(io, verbose, string("reset isparallel = true because of only one worker. ")) 
    end
    msg = string("list of option values: \n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",        
        "chrpairing_phase = ",chrpairing_phase, "\n",        
        "byparent = ", byparent,"\n",                            
        "byneighbor = ", byneighbor,"\n",                        
        "bindis = ", bindis,"cM\n",                                
        "isparentdoseerror = ", isparentdoseerror, "\n",
        "iscorrectparent = ", iscorrectparent, "\n", 
        "isrefineerror = ", isrefineerror, "\n",                
        "tukeyfence = ",tukeyfence,"\n",                
        "doseerror = ", doseerror, "\n",
        "baseerror = ",baseerror,"\n",        
        "allelicbias = ",allelicbias,"\n",        
        "overdispersion = ",overdispersion,"\n",                
        "isgreedy = ",isgreedy,"\n",        
        "ismarkerspecific = ",ismarkerspecific,"\n",        
        "isdelvuong = ", isdelvuong, "\n",
        "delsiglevel = ", delsiglevel, "\n",
        "maxstuck = ",maxstuck, "\n",
        "maxiter = ",maxiter,"\n",
        "minrun = ",minrun, "\n",
        "maxrun = ",maxrun,"\n",
        "runblocksize = ",runblocksize,"\n",
        "refhapfile = ",refhapfile,"\n",
        "missingstring = ",missingstring,"\n",
        "isparallel = ", isparallel, isparallel && nworkers() > 0  ?  string("(nworkers=", nworkers(), ")") : "", "\n",
        "isplot = ",isplot,"\n",
        "outstem = ",isnothing(outstem) ? "no output files" : outstem,"\n",        
        "outext = ",outext,"\n",
        "logfile = ",io,"\n",
        "workdir = ",workdir,"\n",
        "tempdirectory = ",tempdirectory,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
    info_polygeno(polygeno; io, verbose)    
    isnothing(chrsubset) || printconsole(io, verbose, string("chrsubset=",chrsubset))    
    isnothing(snpsubset) || printconsole(io, verbose, string("snpsubset=",snpsubset))    
    getsubpolygeno!(polygeno,chrsubset=chrsubset,snpsubset=snpsubset)    
    if isnothing(isrefineerror)
        isrefineerror = true
        printconsole(io, verbose, string("reset isrefineerror = ", isrefineerror))
    end
    if isnothing(isparentdoseerror)
        maxploidy = maximum(polygeno.polyped.parentinfo[!,:ploidy])
        if maxploidy >= 6 
            isparentdoseerror = false
        else
            isparentdoseerror = true
            # if occursin("AD", join(reduce(vcat,[unique(i[!,:parentformat]) for i in polygeno.markermap])))
            #     isparentdoseerror = false
            # else
            #     isparentdoseerror = true
            # end
        end
        printconsole(io, verbose, string("reset isparentdoseerror = ", isparentdoseerror))
    end
    if isnothing(iscorrectparent) 
        iscorrectparent = true
        printconsole(io, verbose, string("reset iscorrectparent = ", iscorrectparent))
    end
    if isnothing(byparent)            
        avgploidy = mean(polygeno.polyped.parentinfo[!,:ploidy])
        byparent = avgploidy > 4 # byparent =false would be too slow for high ploidy 
        printconsole(io,verbose,string("reset byparent=",byparent))   
    end
    if isnothing(byneighbor)            
        avgploidy = mean(polygeno.polyped.parentinfo[!,:ploidy])
        byneighbor = avgploidy > 4 
        printconsole(io,verbose,string("reset byneighbor=",byneighbor))         
    end    
    parentformatset = get_genoformat(polygeno; isparentformat=true)
    if parentformatset == ["GT_phased"] 
        minrun = maxrun = 1
        printconsole(io,verbose,string("reset minrun=",minrun,", maxrun=",maxrun))         
    end
    if isnothing(maxrun)
        isgbsad = any([occursin("AD", join(reduce(vcat,[unique(i[!,col]) for i in polygeno.markermap]))) 
            for col in ["parentformat", "offspringformat"]])
        maxrun = isgbsad ? minrun : 2*minrun
        printconsole(io,verbose,string("reset maxrun=",maxrun))         
    end
    if isnothing(runblocksize)
        nchr = length(polygeno.markermap)
        runblocksize = nworkers() >= nchr ? minrun : 1
        printconsole(io,verbose,string("reset runblocksize=",runblocksize))         
    end
    if isnothing(chrpairing_phase)            
        avgploidy = mean(polygeno.polyped.parentinfo[!,:ploidy])
        if avgploidy > 4 
            chrpairing_phase = 22
        elseif avgploidy > 2
            chrpairing_phase = 44
        else
            chrpairing_phase = 22
        end
        if avgploidy > 2
            printconsole(io,verbose,string("reset chrpairing_phase=",chrpairing_phase))   
        end
    end
    if isnothing(tukeyfence)        
        popsize = size(polygeno.polyped.offspringinfo,1)
        tukeyfence = popsize >= 50 ? 3 : Inf
        printconsole(io, verbose, string("reset tukeyfence = ",tukeyfence))
    end    

    describe_phase_msg(io,verbose)
    reserrorls  = []
    phasedict = Dict()
    polygenols = connected_polygenos(polygeno)        
    if isplot && (!isnothing(outstem)) && length(polygenols) > 1
        figdir = joinpath(workdir, outstem * "_plots")
        isdir(figdir) || mkdir(figdir)
        for sub = eachindex(polygenols)
            subpolygeno = polygenols[sub]
            ped2 = extract_compactped(subpolygeno.polyped.ped)
            gdesign =  plotped(ped2)            
            fn = joinpath(figdir,string(outstem,"_ped", sub, ".png"))
            savefig(gdesign,fn)            
        end
        msg = string("save design plot in ", figdir)
        printconsole(io,verbose,msg)
    end

    for sub = eachindex(polygenols)
        subpolygeno = polygenols[sub]
        pp = subpolygeno.polyped.parentinfo[!,:individual]
        pp = [i<length(pp) ? string(pp[i],", ") : pp[i] for i=eachindex(pp)]
        msg = string("\nphasing ", sub, " out of ", length(polygenols), 
            " components: ", length(pp), " connected parents = ", string(pp...))
        printconsole(io,verbose,msg)        
        if !byparent
            print_parent_partition(polygeno.polyped; byparent, io,verbose)                
        end        
        phases,resdict_err =parentalphasing_global(subpolygeno; 
            refhapfile, chrpairing_phase, byparent, byneighbor, iscorrectparent,  bindis,  
            isrefineerror, isparentdoseerror, phaseerror_midparent, fixedphaseerror=0, 
            doseerror, baseerror, allelicbias, overdispersion, 
            tukeyfence, isgreedy, ismarkerspecific, isdelvuong, delsiglevel,  
            maxstuck,maxiter,minrun,maxrun,runblocksize, 
            isparallel, workdir, tempdirectory, io,verbose)
        merge!(phasedict, phases)
        push!(reserrorls, resdict_err)
    end    
    isrephase = length(polygenols) > 1        
    isdells = setphase_polygeno!(polygeno, phasedict,reserrorls; isrephase);    
    if isrephase
        pp = polygeno.polyped.parentinfo[!,:individual]
        msg = string("\nphasing all ", length(polygenols), 
            " components: ", length(pp), " parents = ", join(pp,", "))
        printconsole(io,verbose,msg)        
        if !byparent
            print_parent_partition(polygeno.polyped; byparent, io,verbose)                
        end        
        phases,resdict_err =parentalphasing_global(polygeno; 
            refhapfile, chrpairing_phase, byparent, byneighbor, iscorrectparent, bindis, 
            isrefineerror, isparentdoseerror, phaseerror_midparent, fixedphaseerror=0.01, 
            doseerror, baseerror, allelicbias,  overdispersion, 
            tukeyfence, isgreedy,ismarkerspecific, isdelvuong, delsiglevel, 
            maxstuck,maxiter,minrun=1,maxrun=1,runblocksize=1, 
            isparallel, workdir, tempdirectory, io,verbose
        )
        isdells = setphase_polygeno!(polygeno, phases,[resdict_err]; isrephase=false)
    end

    # del markers
    println()
    nmarker = sum(size.(polygeno.markermap,1))
    delmarkerdf = sum(sum.(isdells)) > 0  ? setparentphase!(polygeno,isdells) : DataFrame()    
    ndel = size(delmarkerdf,1)
    if size(delmarkerdf,1)>=1
        msg = string("delete ", ndel," out of ", nmarker, " markers: ",
            join(sort(string.(delmarkerdf[!,:marker])),", "))
        printconsole(io,verbose,msg)
    end          
    # set absphase
    if !isnothing(refhapfile )
        try            
            refparenthap = setabsphase!(refhapfile,polygeno; workdir)
            phaseerrdf = calphaseerr(refparenthap,polygeno)
            savephaseerr(phaseerrdf; io, verbose)
        catch err
            msg = "could not set absolute parental phases,"
            msg = string("warning: ", msg, " ", err)
            verbose && @warn msg
            printconsole(io,false,msg)
        end
    end
    # save res
    if !isnothing(outstem)
        outfile = string(outstem,"_phased"*outext)
        savegenodata(outfile,polygeno,
            missingstring=missingstring,workdir=workdir)
        msg = string("parentphased file: ", outfile)
        printconsole(io,verbose,msg)
        # plot errors
        # figdir = getabsfile(workdir, outstem * "_plots")
        # isdir(figdir) || mkdir(figdir)
        # markermap = reduce(vcat,polygeno.markermap)
        # figerr = plotmarkererror(markermap;tukeyfence=isnothing(tukeyfence) ? 3 : tukeyfence)
        # if !isnothing(figerr)                
        #     errfile = joinpath(figdir,outstem * "_polyphase_markererror.png")
        #     try 
        #         savefig(figerr,getabsfile(workdir,errfile))
        #     catch err
        #         msg = string("Could not savefig. ",err)
        #         @warn msg
        #         printconsole(io, verbose, "Warning: "*msg)
        #     end
        # end
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, io, tused,"polyphase!"; verbose)
    polygeno
end

function setphase_polygeno!(polygeno::PolyGeno,phasedict::AbstractDict,reserrorls::AbstractVector; isrephase::Bool)
    isdells = []    
    hh = [phasedict[i] for i=polygeno.polyped.parentinfo[!,:individual]]
    for chr in eachindex(polygeno.markermap)    
        inputpgeno = deepcopy(polygeno.parentgeno[chr])     
        chrparentgeno = reduce(hcat,[i[chr] for i in hh])        
        polygeno.parentgeno[chr] = Matrix{Any}([ismissing(j) ? j : string.(j) for j in chrparentgeno])
        if isrephase
            b = ismissing.(polygeno.parentgeno[chr])
            if any(b)
                gls = inputpgeno[b]
                formatls = polygeno.markermap[chr][!,:parentformat]
                formatset = unique(formatls)
                if issubset(formatset, ["GT_unphased","GT_phased","GT"])
                    polygeno.parentgeno[chr][b] .= [if ismissing(i) || isa(i, AbstractVector) 
                        i
                    elseif occursin("N",i) 
                        missing         
                    else
                        ls = occursin("/",i) ? split(i, "/") : split(i, "|")
                        allequal(ls) ? parse.(Int8, ls) : i 
                    end for i in gls]
                else
                    error(string("TODO for format=",formatset))
                end    
            end
            isallphased = all(isa.(polygeno.parentgeno[chr], AbstractVector))
            polygeno.markermap[chr][!,:parentformat] .= isallphased ? "GT_phased" : "GT"
        else
            isdel = [any(ismissing.(i)) for i in eachrow(polygeno.parentgeno[chr])]
            push!(isdells,isdel)
            polygeno.markermap[chr][!,:parentformat] .= "GT_phased" 
        end
    end        
    # save error inference
    if length(reserrorls) == 0
        resdict_err = first(reserrorls)
    else
        targetlsls = collect.(keys.(reserrorls))
        allequal(targetlsls) || @error string("unequal targetls=",targetlsls)    
        resdict_err = Dict()
        targetls = first(targetlsls)
        for target in targetls
            push!(resdict_err, target => mean([i[target] for i in reserrorls]))
        end
    end
    for chr in eachindex(polygeno.markermap)    
        for target in targetls
            polygeno.markermap[chr][!,target] .= resdict_err[target][chr] 
        end
    end
    isdells
end

function print_parent_partition(polyped; byparent=false, io=nothing,verbose=true)       
    goodparents, pindexls = get_pindexls(polyped; byparent)
    if byparent    
        printconsole(io,verbose,string("byparent=",byparent, ", partition = ",join(join.(pindexls,"-"),"|")))
    else   
        nsample = length(goodparents) == size(polyped.parentinfo,1) ? 1 : 5
        pindexlsls = unique([last(get_pindexls(polyped; byparent)) for _ in 1:nsample])
        msg = string("goodparents = ",goodparents)
        for i in eachindex(pindexlsls)
            if nsample == 1
                msg *= string("\npartition = ",join(join.(pindexlsls[i],"-"),"|"))
            else
                msg *= string("\nrandpartition",i," = ",join(join.(pindexlsls[i],"-"),"|"))
            end
        end
        printconsole(io,verbose,msg)
    end
end

function describe_phase_msg(io::Union{Nothing, IO},verbose::Bool)
    msg = "keywords in messages: \n"
    msg *= "\t#diff     ≡ #differences between proposal and current phasing\n"
    msg *= "\tΔlogp     ≡ increase of log posterior\n"
    msg *= "\t#stuck    ≡ measure of parental phasing being stuck\n"    
    msg *= "\talwaysapt ≡ if true, always accept proposal phasing\n"    
    msg *= "\tbyblock   ≡ if true, update parental phasing by blocks\n"    
    msg *= "\t#del_err  ≡ #deleted markers with too large error rates\n"
    msg *= "\t#del_test ≡ #deleted markers via Vuong's test\n"
    msg *= "\ttused     ≡ time (second) used by each sub-step in the iteration"
    printconsole(io, verbose,msg)
end

function setparentphase!(polygeno::PolyGeno,isdel::AbstractVector)
    markermap = polygeno.markermap
    nchr = length(markermap)
    # delmarkerdf=vcat([markermap[ch][isdel[ch],[:marker,:linkagegroup,:poscm]] for ch=1:nchr]...)
    delmarkerdf=vcat([markermap[ch][isdel[ch],:] for ch=1:nchr]...)
    isnondel=[.!i for i=isdel]
    polygeno.offspringgeno = [Matrix{Any}(polygeno.offspringgeno[i][isnondel[i],:])  for i=1:length(isnondel)]
    polygeno.parentgeno = [Matrix{Any}(polygeno.parentgeno[i][isnondel[i],:])  for i=1:length(isnondel)]
    polygeno.markermap= [markermap[i][isnondel[i],:]  for i=1:length(isnondel)]
    pushmisc!(polygeno, "deletion"=>delmarkerdf)    
    delmarkerdf
end

function parentalphasing_global(polygeno::PolyGeno;        
    refhapfile::Union{Nothing,AbstractString},
    chrpairing_phase::Union{Nothing,Integer},    
    byparent::Bool, 
    byneighbor::Bool,          
    bindis::Real,
    iscorrectparent::Bool, 
    isrefineerror::Bool,    
    isparentdoseerror::Bool,  
    phaseerror_midparent::Union{Nothing,Real},
    fixedphaseerror::Real,
    doseerror::ErrorSetting,
    baseerror::ErrorSetting,
    allelicbias::ErrorSetting,    
    overdispersion::ErrorSetting, 
    tukeyfence::Union{Nothing,Real}, 
    isgreedy::Bool,
    ismarkerspecific::Union{Nothing,Bool}, 
    isdelvuong::Bool,
    delsiglevel::Real,     
    maxstuck::Integer,
    maxiter::Integer,    
    minrun::Integer,
    maxrun::Integer,
    runblocksize::Integer,     
    isparallel::Bool,    
    workdir::AbstractString,
    tempdirectory::AbstractString,     
    io::Union{Nothing,IO}=nothing, 
    verbose::Bool)        
	# jldopen permissiond denied if tempdir() is in network drive    	
	jltempdir = mktempdir(tempdirectory; prefix="jl_polyorigin_phasing_", cleanup=true)    
	tempid = tempname(jltempdir,cleanup=false)	     
    chridls = [i[1,:linkagegroup] for i in polygeno.markermap]      
    nchr =  length(chridls)
    avgsize = 	mean(length(i) for i in polygeno.offspringgeno)
    tsleep = (isparallel && nworkers() > 1) ? avgsize*(5e-6) : 0 
    # popsize = size(polygeno.polyped.offspringinfo,1)
    
    polygenofilemtx = [begin 
        chr_tsleep = round(Int,tsleep*(chr-1))
        chrid = replace(string(chridls[chr]),"_"=>"-")
        filekey = string(tempid, "_chr", chr, "_", chrid, "_run",run,"_tsleep",chr_tsleep)
        getabsfile(tempdirectory,string(filekey,".csv.gz"))
    end for chr in 1:nchr, run in 1:maxrun]    
    tused = @elapsed for i in 1:nchr
        subpolygeno = getsubpolygeno(polygeno,chrsubset=[i]) 
        savepolygeno(polygenofilemtx[i,1], subpolygeno)
        for j in 2:maxrun        
            cp(polygenofilemtx[i,1],polygenofilemtx[i,j]; force=true)
        end
    end	        
    msg = string("prepare for multiple phasing runs, tused=",round(tused,digits=1), " seconds")
	printconsole(io,verbose,msg)
    
    if isnothing(refhapfile)
        refparenthap = nothing
    else
        refparenthap= readparenthaplo(refhapfile,polygeno.polyped.parentinfo,polygeno.markermap;  workdir)        
    end
    tempid = tempname(jltempdir,cleanup=false)	        
    outfile = getabsfile(tempdirectory, tempid*"_polyorigin_phasing.jld2")        
    phasedict = Dict()      
    tused = @elapsed if isparallel       
        chrls = collect(1:nchr)      
        runstart = 1
        runend = min(minrun, maxrun)
        nwork = nworkers()
        while true        
            while true
                ntask = length(chrls) * (runend - runstart + 1)
                if ntask >= nwork                     
                    break
                else                    
                    if (ntask + length(chrls) > nwork)                                                
                        break
                    else
                        runend = min(maxrun, runend+1)
                        runend == maxrun && break
                    end
                end
            end     
            (runend + 1) == maxrun && (runend += 1)                   
            msg = string("parallel phasing: chrls=",chrls, ", run=",runstart:runend)
            printconsole(io, verbose, msg)
            polygenofilels = vec(permutedims(polygenofilemtx[chrls,runstart:runend])) 
            # isrefineerror2 = runstart <= min(minrun,maxrun) ? isrefineerror : false  # no error estimation would result in smaller logp!!          
            parentalphasing_local(polygenofilels; chrpairing_phase, byparent, byneighbor, bindis,  
                iscorrectparent, tukeyfence,                                 
                isrefineerror,                 
                isparentdoseerror, phaseerror_midparent,fixedphaseerror, 
                doseerror, baseerror, allelicbias, overdispersion,                                  
                ismarkerspecific,isdelvuong, delsiglevel, isgreedy, 
                maxstuck, maxiter, isparallel, outfile, io, verbose)		             
            update_phasedict_chrls!(phasedict, chrls, chridls, outfile,polygenofilemtx,runend; 
                isrefineerror, 
                refparenthap, io, minrun, maxrun,verbose)      
            (isempty(chrls) || runend == maxrun) && break
            runstart = runend + 1
            runend = min(maxrun, runstart + runblocksize - 1)
        end
    else
        for chr in eachindex(chridls)
            chrls = [chr]            
            for run in 1:maxrun
                polygenofilels = [polygenofilemtx[chr,run]]                
                # isrefineerror2 = run <= min(minrun,maxrun) ? isrefineerror : false                
                parentalphasing_local(polygenofilels; chrpairing_phase, byparent, byneighbor,  bindis, 
                    iscorrectparent, tukeyfence,                                 
                    isrefineerror,  
                    isparentdoseerror, phaseerror_midparent,fixedphaseerror, 
                    doseerror, baseerror, allelicbias, overdispersion,  
                    ismarkerspecific, isdelvuong, delsiglevel, isgreedy, 
                    maxstuck, maxiter, isparallel=false, outfile, io, verbose)	
                if run >= min(minrun, maxrun)
                    update_phasedict_chrls!(phasedict, chrls, chridls, outfile,polygenofilemtx,run; 
                        isrefineerror, 
                        refparenthap, io, minrun, maxrun,verbose)                      
                    isempty(chrls) && break                   
                end
            end
        end  
    end
    # phasedict[chr]: phaplo, bvpair, logl, logprior_error,logl_parentdata, ndel, targetls, doseerror, baseerror, allelicbias, overdispersion,maxlogp                
    # error inference
    targetindex = 7
    targetls_allchr = [phasedict[chr][targetindex] for chr in 1:nchr]    
    allequal(targetls_allchr) || @error string("unequal targetls = ",targetls_allchr)
    targetls = first(targetls_allchr)
    doseerror_allchr = [phasedict[chr][targetindex+1] for chr in 1:nchr]    
    baseerror_allchr = [phasedict[chr][targetindex+2] for chr in 1:nchr]    
    allelicbias_allchr = [phasedict[chr][targetindex+3] for chr in 1:nchr]    
    overdispersion_allchr = [phasedict[chr][targetindex+4] for chr in 1:nchr]    
    maxlogp_allchr = [phasedict[chr][targetindex+5] for chr in 1:nchr]
    msg = string("tused=", round(tused, digits=1), 
        "s for all phasing runs, logp sum =",round(sum(maxlogp_allchr),digits=2), 
        ", logp perchr = ", round.(maxlogp_allchr, digits=2))
    printconsole(io, verbose, msg)
    resdict_err = Dict()
    in("doseerror",targetls) && push!(resdict_err, "doseerror" => doseerror_allchr)
    in("baseerror",targetls) && push!(resdict_err, "baseerror" => baseerror_allchr)
    in("allelicbias",targetls) && push!(resdict_err, "allelicbias" => allelicbias_allchr)    
    in("overdispersion",targetls) && push!(resdict_err, "overdispersion" => overdispersion_allchr)    
    # parent phasing
    parenthaplo= [collect(eachcol(first(phasedict[chr]))) for chr in 1:nchr]    
    parents = polygeno.polyped.parentinfo[!,:individual]    
    resdict_parent = Dict([parents[i] => [j[i] for j=parenthaplo] for i=eachindex(parents)])
    resdict_parent, resdict_err
end

function update_phasedict_chrls!(phasedict, chrls, chridls, outfile,polygenofilemtx,runend; 
    refparenthap::Union{Nothing,AbstractVector},
    isrefineerror::Bool, 
    io::Union{Nothing,IO}, minrun::Integer, maxrun::Integer,
    verbose::Bool)       
    for chr in chrls
        # resls: phasing results for all runs
        # resls[run]: phaplo, bvpair, logl, logprior_error, logl_parentdata,  ndel, targetls, doseerrorls, baseerrorls, allelicbiasls, overdispersionls         
        # phaplo: nmaker x nparent, phased parent genotypes, set to missing at deleted markers.                
        # logl_parentdata: list of logprior for each parent; set to nothing for fixed parent                
        resdict = load(outfile,string("chr",chr))
        resls = last.(sort(collect(resdict), by = x->parse(Int,replace(x[1],"run"=>""))))
        # doseerrorls, baseerrorls, allelicbiasls, overdispersionls                 
        phaplols = first.(resls)        
        loglls = [val[3] for val in resls] # logl inluding logprior_error
        loglerror_ls = [val[4] for val in resls]       
        
        logpls = loglls .+  loglerror_ls
        nparent = length(resls[1][5]) # 5: logl_parentdata is a list of logprior for each parent
        fixloglparentls = [begin 
            fixkey = (chr, p, "fixedparent")
            haskey(phasedict, fixkey) ? phasedict[fixkey] : nothing
        end for p in 1:nparent]
        loglparent_ls = [val[5] for val in resls]       
        for ls in loglparent_ls 
            fixls = findall(isnothing.(ls))
            # println("chr=",chr, ",fixls=",fixls, ", fixloglparentls[fixls] = ", fixloglparentls[fixls] )
            if !isempty(fixls)            
                ls[fixls] .= fixloglparentls[fixls] 
            end
            ls
        end         
        if !any(isnothing.(reduce(vcat,loglparent_ls)))
            logpls .+= sum.(loglparent_ls)
        end        
        ndells = [val[6] for val in resls]                  
        maxlogp, best = findmax(logpls)        
        deltlogpls = maxlogp .- logpls        
        if isnothing(refparenthap)
            refmsg = ""
        else
            truechrhaplo = refparenthap[chr]
            truedis = [sum.(PolyOrigin.cal_haplo_nmismatch(truechrhaplo, est)) for est=phaplols]            
            truedis_phase = first.(truedis)
            truedis_perm = last.(truedis)
            refmsg = string("true_nphaseerr(refhap ~ runs) = ", truedis_phase,", true_npermerr = ", truedis_perm)
        end
        bestphaplo = phaplols[best]
        nmarker,nparent = size(bestphaplo)
        thresdiff = round(Int,0.001*nmarker)
        dis = [PolyOrigin.cal_haplo_nmismatch(phaplols[best],est) for est=phaplols]
        disls_phase = first.(dis)        
        disls_perm =  last.(dis)
        fixedparentls = findall(sum([i .<= thresdiff for i in disls_perm]) .>= max(2,minrun)) 
        if 1 <= length(fixedparentls) < nparent || isrefineerror           
            # doseerrorls, baseerrorls, allelicbiasls, overdispersionls         
            besterrls = resls[best][end-3:end]
            remainfilels = polygenofilemtx[chr,runend+1:end]            
            update_remainfilels(remainfilels, bestphaplo, besterrls, fixedparentls)
            for p in fixedparentls
                fixkey = (chr, p, "fixedparent")
                if !haskey(phasedict, fixkey)
                    phasedict[fixkey]  = loglparent_ls[best][p]
                end
            end
        end
        dis_phase = sum.(disls_phase)
        dis_perm = sum.(disls_perm)
        # ncount = sum(dis_phase .== 0)
        ncount = sum(dis_perm[i] == 0 || deltlogpls[i] < 1 for i = eachindex(dis_perm,deltlogpls))        
        todparentls = setdiff(1:nparent,fixedparentls)
        summsg = string(", parents = ", isempty(fixedparentls) ? "all" : (isempty(todparentls) ? "none" : join(todparentls,"|")),             
            ", nphaseerr(run",best," ~ runs) = ", dis_phase,
            ", npermerr = ", dis_perm,                    
            ", ndel = ", ndells, 
            "\n\tlogl = ",round.(loglls,digits=2),                        
            ", logp = ",round.(logpls,digits=2),                        
            # "\n\tloglerror_ls = ",round.(loglerror_ls,digits=2),                        
            # "loglparent_ls = ",round.(sum.(loglparent_ls),digits=2),                        
        )
        # @info string("chr=", chridls[chr]) loglparent_ls fixloglparentls
        if ncount >= minrun || isempty(todparentls) || length(logpls) == maxrun 
            msg = string("chr=", chridls[chr], ", done", summsg)                 
            isempty(refmsg) || (msg *= "\n\t"*refmsg)
            printconsole(io,verbose,msg)                        
            resbest = (resls[best]...,maxlogp)
            push!(phasedict, chr=>resbest)
        else
            msg = string("chr=", chridls[chr], ", to_be_continued", summsg)
            isempty(refmsg) || (msg *= "\n\t"*refmsg)
            printconsole(io,verbose,msg)
        end                        
    end
    setdiff!(chrls,keys(phasedict))
    phasedict, chrls
end

function update_remainfilels(remainfilels::AbstractVector, bestphaplo::AbstractMatrix,
    besterrls, fixedparentls::AbstractVector)
    isempty(remainfilels) && return nothing
    firstgenofile =remainfilels[1]
    polygeno = readpolygeno(firstgenofile)    
    if !isempty(fixedparentls) || true
        ploidyls = polygeno.polyped.parentinfo[!,:ploidy] 
        for p in fixedparentls
            ploidy = ploidyls[p]
            polygeno.parentgeno[1][:,p] .= [ismissing(i) ? join(repeat(["N"],ploidy),"|") : string.(i) for i in bestphaplo[:,p]] # only one chromosome in each remainfilels    
        end
        polygeno.markermap[1][!, :parentformat] .= [if occursin(r"^GT",i) 
            "GT"
        elseif in(i, ["AD","GP","GL"])
            string(i,"&GT")    
        elseif occursin("&",i) 
            ls = sort(split(i,"&")) # sort to put GT in the end
            length(ls) == 2 || @error string("TODO for parentformat=",i)
            if occursin(r"^GT",ls[1]) 
                if occursin(r"^GT",ls[2]) 
                    "GT"
                else
                    @error string("TODO for parentformat=",i)
                    string(ls[i],"&GT")    
                end
            elseif occursin(r"^GT",ls[2]) 
                string(ls[1],"&GT")    
            else
                @error string("TODO for parentformat=",i)
                "mixed"
            end
        else
            @error string("TODO for parentformat=",i)
            "mixed"
        end for i in polygeno.markermap[1][!, :parentformat]]
    end
    # besterrls: doseerrorls, baseerrorls, allelicbiasls, overdispersionls         
    colls = [:doseerror, :baseerror, :allelicbias, :overdispersion]    
    for i in eachindex(colls)
        polygeno.markermap[1][!, colls[i]] .= besterrls[i] 
    end
    isdells = [all(ismissing.(i)) for i in eachrow(bestphaplo)]      
    polygeno.markermap[1][!,:filter] .= map((isdel,y)->if isdel 
        if ismissing(y) 
            "POLYORIGINDEL"
        else
            if occursin(r"^POLYORIGINDEL",y) 
                y
            else
                string("POLYORIGINDEL;",y)
            end
        end
    else
        y
    end, isdells, polygeno.markermap[1][!,:filter])    
    savepolygeno(firstgenofile,polygeno)
    for dst in remainfilels[2:end]
        cp(firstgenofile,dst; force=true) 
    end
    return nothing
end

function cal_haplo_nmismatch(refhap::AbstractMatrix, esthap::AbstractMatrix)
    size(refhap) == size(esthap) || @error "inconsistent sizes"
    nind = size(refhap,2)
    nmismatch_phase = zeros(Int, nind)
    nmismatch_perm = zeros(Int, nind)
    b = [any(ismissing.(reduce(vcat,i))) for i in eachrow(refhap)] 
    b = b .|| [any(ismissing.(reduce(vcat,i))) for i in eachrow(esthap)] 
    isnonmiss = .!b
    for p in 1:nind
        ref_p = view(refhap,isnonmiss, p)
        est_p = view(esthap,isnonmiss, p)
        ploidy = length(ref_p[1])
        ploidy == length(est_p[1]) || @error "inconsistent ploidy"
        permls = collect(Combinatorics.permutations(1:ploidy))
        dls = [sum(ref_p .!= [j[perm] for j in est_p]) for perm in permls]        
        d,pos = findmin(dls)
        nmismatch_phase[p] = d
        if d == 0
            npermerr = 0
        else
            perm = permls[pos]
            est_p2 = [j[perm] for j in est_p]
            (; npermerr) = calphaseerr(ref_p,est_p2)           
        end
         nmismatch_perm[p] = npermerr
    end
    [nmismatch_phase,nmismatch_perm]
end

function parentalphasing_local(polygenofilels::AbstractVector;
    chrpairing_phase::Union{Nothing,Integer},    
    byparent::Bool,
    byneighbor::Bool,    
    bindis::Real, 
    iscorrectparent::Bool, 
    isrefineerror::Bool,   
    isparentdoseerror::Bool,
    phaseerror_midparent::Union{Nothing,Real},
    fixedphaseerror::Real,
    doseerror::ErrorSetting,
    baseerror::ErrorSetting,
    allelicbias::ErrorSetting,    
    overdispersion::ErrorSetting, 
    tukeyfence::Real, 
    isgreedy::Bool,
    ismarkerspecific::Union{Nothing,Bool},
    isdelvuong::Bool,
    delsiglevel::Real,        
    maxstuck::Integer,
    maxiter::Integer,    
    isparallel::Bool,    
    outfile::AbstractString, 
    io::Union{Nothing, IO}, 
    verbose::Bool)    
	startt = time()	    
    chrrunls = [begin 
        filesegs = split(first(split_allext(basename(genofile))),"_") 
        run = replace(filesegs[end-1],"run"=>"") # filesegs[end] for tsleep
        chrid = filesegs[end-2]
        chr = replace(filesegs[end-3], "chr"=>"")
        (parse(Int, chr), chrid, run)
    end for genofile in polygenofilels]
    jldopen(outfile,"a") do file        
        jltempdir = dirname(first(polygenofilels))
        tempid = tempname(jltempdir,cleanup=true)
        tempfilels = [string(tempid, "_datalike_chr",chr, "_run",run,".jld2") for (chr, chrid, run) in chrrunls]
        if isparallel && nprocs()>1                                 
            logfilels = [string(tempid, "_chr",chr, "_run",run,".log") for (chr, chrid, run) in chrrunls]
            try
                resls = pmap((x,y,z)->parentalphasing_chr_run(x; chrpairing_phase, byparent, byneighbor, bindis,  
                    iscorrectparent, isrefineerror, isparentdoseerror, phaseerror_midparent, fixedphaseerror,
                    doseerror, baseerror, allelicbias,  overdispersion, chrliketempfile=z, 
                    tukeyfence, isgreedy,ismarkerspecific, isdelvuong, delsiglevel, 
                    maxstuck, maxiter, logio=y, verbose), polygenofilels, logfilels,tempfilels)				                    
                for i in eachindex(resls, logfilels,chrrunls)	
                    if isfile(logfilels[i])			
                        chr, chrid, run = chrrunls[i]
                        file[string("chr", chr, "/run",run)] = resls[i]
                    end
                end                    
            catch err 				
                printconsole(io,false, string("ERROR: ",err))			
                rethrow(err)
            finally
                if !isnothing(io)
                    for f in logfilels
                        isfile(f) && write(io,read(f))		                        
                    end
                end
                rm.(logfilels,force=true)           
            end
        else			
            # runres: phaplo, bvpair,logp, logphis	                                    
            for i in eachindex(polygenofilels,chrrunls,tempfilels)				
                runres = parentalphasing_chr_run(polygenofilels[i]; chrpairing_phase, byparent, byneighbor,  bindis,
                    iscorrectparent, isrefineerror, isparentdoseerror, fixedphaseerror,
                    phaseerror_midparent,doseerror, baseerror, allelicbias,overdispersion,   
                    tukeyfence, isgreedy, ismarkerspecific, isdelvuong, delsiglevel, chrliketempfile=tempfilels[i],
                    maxstuck, maxiter, logio=io, verbose)		
                chr, chrid, run = chrrunls[i]
                file[string("chr", chr, "/run",run)] = runres
            end            
        end					        
        if length(chrrunls) == 1
            chrrunls2 = string(chrrunls[1][2], "_run",chrrunls[1][3])
        else
            sort!(chrrunls, by = x->x[1])        
            iils = splitindex(first.(chrrunls))
            chrrunls2 = [string(chrrunls[ii[1]][2],"_runs",chrrunls[ii[1]][3],"-",chrrunls[ii[end]][3]) for ii in iils]
        end
        msg = string("tused =", round(time()-startt,digits=1), "s for phasing runs = ", chrrunls2)
        printconsole(io,verbose,msg)		
        rm.(tempfilels,force=true)     	
    end
end

function parentalphasing_chr_run(polygenofile::AbstractString; 
    chrpairing_phase::Union{Nothing,Integer},            
    byparent::Bool,
    byneighbor::Bool,    
    bindis::Real,
    iscorrectparent::Bool, 
    isrefineerror::Bool,    
    isparentdoseerror::Bool,
    phaseerror_midparent::Union{Nothing,Real},
    fixedphaseerror,
    doseerror::ErrorSetting,
    baseerror::ErrorSetting,
    allelicbias::ErrorSetting,    
    overdispersion::ErrorSetting, 
    tukeyfence::Real, 
    isgreedy::Bool,
    ismarkerspecific::Union{Nothing,Bool},
    isdelvuong::Bool,
    delsiglevel::Real,                
    maxstuck::Integer,
    maxiter::Integer,      
    chrliketempfile::AbstractString,
    logio::Union{Nothing,AbstractString,IO},    		
    verbose::Bool)    
    io = isa(logio,AbstractString) ? open(logio,"w") : logio	
	try	
        filesegs = split(first(split_allext(basename(polygenofile))),"_") 
        tsleep = parse(Int, replace(filesegs[end],"tsleep"=>""))
        run = parse(Int, replace(filesegs[end-1],"run"=>""))        		
		sleep(tsleep) # avoid reading data in parallel		        
		tused = @elapsed polygeno = readpolygeno(polygenofile)
		chr = 1      # polygeno contains data fdor a single chromosome
	    chrid = polygeno.markermap[chr][1,:linkagegroup]				        
		msg_run = string(", run=",run)
		nmarker = size(polygeno.markermap[chr],1)
        chrparentformat = unique(polygeno.markermap[chr][!,:parentformat])
        chroffformat = unique(polygeno.markermap[chr][!,:offspringformat])
		printconsole(io,verbose,string("chr=", chrid, msg_run, ", #markers=", nmarker, 
            ", parentformat=", String.(chrparentformat), ", offspringformat=", String.(chroffformat), 
            ", readpolygeno, tused=", round(tused,digits=1),"s"))		                        
        # fill_parent_missing!(polygeno,chr)
        res = parentalphasing_chr_run(polygeno, chr; chrpairing_phase, byparent, byneighbor,  bindis, 
            iscorrectparent, isrefineerror, isparentdoseerror, fixedphaseerror,
            phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, 
            tukeyfence, isgreedy, ismarkerspecific, isdelvuong, delsiglevel, 
            maxstuck,maxiter,run, chrliketempfile, io,verbose)                
        res
	finally
		isa(logio,AbstractString) && close(io)		
	end    
end    

function parentalphasing_chr_run(polygeno::PolyGeno, chr::Integer;     
    chrpairing_phase::Union{Nothing,Integer},              
    byparent::Bool,
    byneighbor::Bool,        
    bindis::Real,
    iscorrectparent::Bool, 
    isrefineerror::Bool,    
    isparentdoseerror::Bool,  
    phaseerror_midparent::Union{Nothing,Real},
    fixedphaseerror::Real, 
    doseerror::ErrorSetting,
    baseerror::ErrorSetting,
    allelicbias::ErrorSetting,    
    overdispersion::ErrorSetting, 
    tukeyfence::Union{Nothing,Real}, 
    isgreedy::Bool,
    ismarkerspecific::Union{Nothing,Bool},
    isdelvuong::Bool,
    delsiglevel::Real,            
    run::Integer, 
    maxstuck::Integer,
    maxiter::Integer,        
    chrliketempfile::AbstractString,
    io::Union{Nothing,IO},
    verbose::Bool)            
    chrid = polygeno.markermap[chr][1,:linkagegroup]
    chrparentgeno = polygeno.parentgeno[chr] # keep phased genotypes in format of vector{String}    
    chrparentformat = polygeno.markermap[chr][!,:parentformat]         
    chroffgeno = get_chroffgeno(polygeno,chr)      
    chrmultiallele = polygeno.markermap[chr][!,:multiallele]
    chroffformat = polygeno.markermap[chr][!,:offspringformat]                 
    # setting for error inference
    targetls = doseerror.isinfer ?  ["doseerror"] : String[] 
    hasoffAD = occursin("AD",join(chroffformat))    
    if hasoffAD
        baseerror.isinfer && push!(targetls,"baseerror")
        allelicbias.isinfer && push!(targetls,"allelicbias")
        overdispersion.isinfer && push!(targetls,"overdispersion")
    end        
    isrefineerror = isrefineerror && !isempty(targetls)    
    errorprior = (doseerror=doseerror.prior, baseerror = baseerror.prior, 
        allelicbias = allelicbias.prior, overdispersion = overdispersion.prior)
    errorsoftthresh = (doseerror = doseerror.softthresh, baseerror = baseerror.softthresh, 
        allelicbias = allelicbias.softthresh, overdispersion = overdispersion.softthresh)
    errorhardthresh = (doseerror = doseerror.hardthresh, baseerror = baseerror.hardthresh, 
        allelicbias = allelicbias.hardthresh, overdispersion = overdispersion.hardthresh)            
    # doseerror.init
    hasparentAD = occursin("AD",join(chrparentformat))    
    doseerrinit = isnothing(doseerror.init) ? (hasparentAD ? 0.001 : 0.01) : doseerror.init
    estdoseerror = [ismissing(i) ? doseerrinit : i for i in polygeno.markermap[chr][!,:doseerror]]    
    estbaseerror = [ismissing(i) ? baseerror.init : i for i in polygeno.markermap[chr][!,:baseerror]]
    estallelicbias = [ismissing(i) ? 1.0*allelicbias.init : i for i in polygeno.markermap[chr][!,:allelicbias]]
    estoverdispersion = [ismissing(i) ? overdispersion.init : i for i in polygeno.markermap[chr][!,:overdispersion]]            
    if isnothing(ismarkerspecific)
        fixtargetls = intersect(["doseerror","overdispersion","baseerror"],targetls)
        vartargetls = intersect(["allelicbias"],targetls)
        estdoseerror = mean(estdoseerror)
        estbaseerror = mean(estbaseerror)        
        estoverdispersion = mean(estoverdispersion)            
    else
        if ismarkerspecific
            fixtargetls = []
            vartargetls = intersect(["doseerror","allelicbias", "overdispersion","baseerror"],targetls)            
        else            
            fixtargetls = intersect(["doseerror","allelicbias", "overdispersion","baseerror"],targetls)            
            vartargetls = []
            estdoseerror = mean(estdoseerror)
            estbaseerror = mean(estbaseerror)
            estallelicbias = mean(estallelicbias)
            estoverdispersion = mean(estoverdispersion)            
        end
    end        
    
    # calculate prior and inital
    polyped = polygeno.polyped    
    priorspace = get_priorstatespace(polyped,chrpairing_phase)    
    maxvalent = isnothing(chrpairing_phase) ? 4 : max(digits(chrpairing_phase)...)    
    priorprocess = get_priorprocess(polygeno,chr,maxvalent; bindis)   
    isdells = [ismissing(i) ? false : occursin(r"^POLYORIGINDEL",i) for i in polygeno.markermap[chr][!,:filter]]        
    isnondells = .!isdells
    if any(isdells)        
        msg = string("chr=", chrid, ", run=", run, ", delete ", sum(isdells),"/", length(isdells), " markers before phasing")
        printconsole(io,verbose,msg)        
        setpriorprocess!(priorprocess,findall(isdells),bindis)
    end        
    
    if isrefineerror && in("doseerror",targetls)
        if isa(estdoseerror,AbstractVector)  
            estdoseerror .=  [i > 1e-4 ? rand(Uniform(1e-4, i)) : i for i in estdoseerror]
        else 
            if estdoseerror > 1e-4 
                estdoseerror = rand(Uniform(1e-4, estdoseerror)) 
            end
        end
    end
    # if isrefineerror && in("baseerror",targetls)
    #     if isa(estbaseerror,AbstractVector) 
    #         estbaseerror .=   [i > 1e-5 ? rand(Uniform(1e-5,i)) : i for i in estbaseerror] 
    #     else
    #         if estbaseerror > 1e-5
    #             estbaseerror = rand(Uniform(1e-5, estbaseerror)) 
    #         end
    #     end
    # end       
    phaploset, phaploweight = init_phaploprior(chrparentgeno, chrparentformat, chrmultiallele,polyped; 
        phaseerror = fixedphaseerror,
        doseerror = isparentdoseerror ? estdoseerror : 0.0,  
        baseerror = estbaseerror, 
        allelicbias = estallelicbias, overdispersion  = estoverdispersion)    
    phaploindex= init_phaploindex(phaploweight)    
    # phaploposterior:  posterior distribution of h^m after account for its parents' haplotypes, same as phaploweight if m is a founder. 
    phaploposterior  = deepcopy(phaploweight)
    if fixedphaseerror > 0 
        msg = string("TODO: set fixed parents when phaseerror > 0")
        @warn msg
        printconsole(io, false,msg)
    end
    fixedparentls = findall([unique(length.(view(i,isnondells))) == [1] for i in phaploset])        
    nparent = size(chrparentgeno,2)
    for p in 1:nparent
        gls = chrparentgeno[isnondells ,p]
        b = [isa(gls[t],AbstractVector) && eltype(gls[t]) <: AbstractString && occursin("GT", chrparentformat[t]) for t in eachindex(gls, chrparentformat[isnondells])] 
        if all(b)
            if !in(p, fixedparentls) 
                msg = string("phased parent = ",p, " is not in the fixedparentls")
                @warn msg
                printconsole(io, false,msg)
            end
        end 
    end
    popidls = collect(skipmissing([issubset(priorspace[i]["parentindex"],fixedparentls) ? missing : i for i in keys(priorspace)]))
    bvpair = init_bvpair(priorspace,polyped)
    siblogl = init_siblogl(bvpair)    
    rand_update_bvpair!(siblogl, bvpair, phaploindex,phaploset, priorspace, priorprocess,
        chroffgeno,chroffformat, chrmultiallele, polyped; 
        phaseerror_midparent, byneighbor, isrand = true, incldelmarker=false, 
        doseerror = estdoseerror, baseerror = estbaseerror, 
        allelicbias = estallelicbias, overdispersion = estoverdispersion,           
    )  

    msg_err = get_msg_err(targetls, priorprocess; 
        doseerror = estdoseerror, baseerror = estbaseerror, 
        allelicbias = estallelicbias, overdispersion = estoverdispersion)
    msg = string("chr=", chrid, ", run=", run, ", targetls=", targetls, 
        ", #fixedparents=", length(fixedparentls), msg_err)        
    printconsole(io,verbose,msg)                 
    
    correct_minerrdiff = min(3,max(1, round(Int,size(chroffgeno,2)/5)))      

    ndel_err = -1
    nstuck = 0  
    nstuck_correct = 0   
    it_1ststuck = -1
    ncorrect = iscorrectparent ? -1 : 0    
    ndiff = ntot = 2*length(chrparentgeno) 
    ndiffls = [2*ndiff,ndiff]
    israndls = falses(size(chrparentgeno,2))     
    alwaysaccept = isgreedy ? false : true
    upbyblock  = false    
    # upbyhalf = false
    istucking = false    
    tusedls = []            
    for it=1:maxiter        
        # starttime = time()                      
        # update phaloindex, phaploposterior, and bvpair        
        empty!(tusedls)      
        msg = string("chr=",chrid,", run=", run,", it=", it)
        if issetequal(1:nparent, fixedparentls)
            ndiff = 0            
            ntot = 0
            deltlogl = 0.0
            isaccept = true
        else
            tused = @elapsed ndiff,ntot, deltlogl, isaccept = update_phaplo_bvpair!(siblogl, 
                phaploindex, phaploset, phaploweight, phaploposterior, 
                priorspace, priorprocess, bvpair, 
                chrparentgeno, chrparentformat, chrliketempfile, 
                chroffgeno, chroffformat, chrmultiallele, polyped; 
                israndls, isreversechr = iseven(it), alwaysaccept, upbyblock, upbyhalf=upbyblock,  
                phaseerror_midparent, byneighbor, byparent,  fixedparentls, 
                doseerror = estdoseerror, baseerror = estbaseerror, 
                allelicbias = estallelicbias, overdispersion = estoverdispersion,           
            )   
            push!(tusedls,tused)
        end        
        if alwaysaccept
            deltlogl < 0 && (nstuck += 1)
        else
            nstuck = (deltlogl > 0 && ndiff > 0) ? 0 : (nstuck +1 + (ndiff == 0))
        end     
        if !issetequal(1:nparent, fixedparentls)
            msg *= string(", #diff=",join([ndiff,ntot],"/"))  
            msg *= string(", Δlogp=", round(deltlogl,digits=1), isaccept ? "" : "(rej)")
            msg *= string(", #stuck=", nstuck, ", alwaysaccept=", alwaysaccept, ", byblock=", upbyblock)        
        end
        if !istucking
            istucking = !alwaysaccept && (nstuck > 0 || (it >= maxiter - 1) || ndiff == 0)
            istucking && (it_1ststuck = it)              
        end        
        if isrefineerror && !isempty(targetls) && (!upbyblock || (upbyblock && (nstuck == 0 || it - it_1ststuck <= 0) || ndel_err>0))            
            vartargetls2 = istucking ? vartargetls : setdiff(vartargetls,["doseerror","baseerror"])            
            t_var= @elapsed if isempty(vartargetls2)                
                ndel_err = -1
            else
                errorbound = (it == 1 || (it == it_1ststuck) && in("doseerror", vartargetls2)) ? log(100) : log(5)
                ndel_err = last(infer_errorls!(chrparentgeno, chrparentformat, 
                    phaploindex,phaploset, chroffgeno, chroffformat, chrmultiallele, 
                    bvpair, priorspace, priorprocess, polyped, vartargetls2; 
                    phaseerror_midparent, isdellargeerr = istucking, 
                    doseerror = estdoseerror, baseerror = estbaseerror, 
                    allelicbias = estallelicbias, overdispersion = estoverdispersion,           
                    errorprior, errorsoftthresh,tukeyfence, errorhardthresh, bindis, 
                    errorbound,  errormaxit = 20,
                ))  
                if it >= 2 
                    errorprior = update_errorprior(errorprior; vartargetls = vartargetls2, estdoseerror, estbaseerror,estallelicbias,estoverdispersion)                         
                end                
                # @info tostring_errorprior(errorprior,targetls)                                                   
            end            
            fixtargetls2 = istucking ? fixtargetls : setdiff(fixtargetls,["doseerror","baseerror"])
            t_fix = @elapsed if !isempty(fixtargetls2)
                errorbound = (it == 1 || (it == it_1ststuck) && in("doseerror", fixtargetls2)) ? log(100) : log(5)
                estdoseerror, estbaseerror, estallelicbias,estoverdispersion = infer_avgerror(fixtargetls2, 
                    phaploindex, phaploset,  chrparentgeno, chrparentformat, 
                    chroffgeno, chroffformat, chrmultiallele, bvpair, priorspace, priorprocess,polyped;     
                    phaseerror_midparent, errorprior,  
                    doseerror = estdoseerror, baseerror = estbaseerror, 
                    errorbound,  errormaxit = 20,
                    allelicbias = estallelicbias, overdispersion = estoverdispersion,           
                )                    
            end
            t_up = @elapsed if !(targetls == ["doseerror"] && !isparentdoseerror)                
                updoseerror = isparentdoseerror ? estdoseerror : 0.0                
                up_parent_phaploprior!(phaploindex, phaploset, phaploweight, phaploposterior, chrparentgeno, chrparentformat, 
                    chrmultiallele,polyped; phaseerror = fixedphaseerror,
                    doseerror = updoseerror, 
                    baseerror = estbaseerror, 
                    allelicbias = estallelicbias, overdispersion = estoverdispersion,           
                )                                     
            end
            msg_err = get_msg_err(targetls, priorprocess; 
                doseerror = estdoseerror, baseerror = estbaseerror, 
                allelicbias = estallelicbias, overdispersion = estoverdispersion)
            if ndel_err >= 0 
                msg_err *= string(istucking ? ", #del_err=" : ", #large_err=", ndel_err) 
            end              
            tused = t_var + t_fix + t_up           
            push!(tusedls,tused)          
            empty!(lru_genolike_GTunphased)
            empty!(lru_genolike_GTphased)    
            empty!(lru_genolike_AD)          
        else
            msg_err = ""                            
        end
        if iscorrectparent && istucking
            # modify phaploindex, phaploset, and phaloweight
            # corrections might result in phases with zero weight; In cal_logl_parentdata, zero weight is replaced with 1e-8                      
            tused = @elapsed changelogp, correctdf = parentcorrect!(phaploindex, phaploset, phaploweight, phaploposterior, 
                chroffgeno, chroffformat, chrmultiallele, priorspace, priorprocess, polyped; 
                bvpair,                
                phaseerror_midparent, callthreshold = 0.7, minerrdiff = correct_minerrdiff, 
                doseerror = estdoseerror, baseerror = estbaseerror, 
                allelicbias = estallelicbias, overdispersion = estoverdispersion,      
            )                              
            correctapt = changelogp > 0
            push!(tusedls,tused)
            ncorrect = size(correctdf,1)
            msg_correct = string(", #correct=", ncorrect)
            ncorrect > 0 && (msg_correct *= correctapt ? "(accepted)" : "(rejected)")
            if ncorrect == 0 || !correctapt 
                nstuck_correct += 1
            else
                nstuck_correct = 0
            end
            nstuck_correct >= 2 && (iscorrectparent = false)            
        else
            correctapt = false
            msg_correct = ""
        end

        if isdelvuong && istucking 
            # priorprocess is modified
            tused = @elapsed deltt = polymarkerdel!(bvpair, phaploindex,phaploset, priorspace, priorprocess,
                chroffgeno, chroffformat, chrmultiallele; delsiglevel, phaseerror_midparent, 
                doseerror = estdoseerror, baseerror = estbaseerror, 
                allelicbias = estallelicbias, overdispersion = estoverdispersion, bindis,         
            )             
            push!(tusedls,tused)         
            ndel= length(deltt)               
            msg_del = string(", #del_test=", ndel)
            ndel == 0 && (isdelvuong = false)
        else
            ndel = -1
            msg_del = ""         
        end                                  
        if !all(isempty.([msg_err, msg_correct, msg_del]))            
            tusedls[1] += @elapsed rand_update_bvpair!(siblogl, bvpair, phaploindex,phaploset, priorspace, priorprocess,
                chroffgeno,chroffformat, chrmultiallele, polyped;                 
                phaseerror_midparent, byneighbor, isrand = true, incldelmarker=false, 
                doseerror = estdoseerror, baseerror = estbaseerror, 
                allelicbias = estallelicbias, overdispersion = estoverdispersion,  popidls,      
            )            
        end
        msg *= string(msg_err,  msg_correct, msg_del)    
        msg *= string(", tused=", join(round.(tusedls,digits=1),"|"),"s")  
        printconsole(io,verbose,msg)     

        push!(ndiffls, ndiff)
        if !alwaysaccept && upbyblock && ndel <= 0  && ndel_err <= 0         
            if (nstuck >= maxstuck && !iscorrectparent) || it == maxiter || (ndiff == 0 && ndiffls[end-1] == 0) 
                break
            end            
        end        
        # update alwaysaccept and upbyblock        
        if alwaysaccept            
            if ndiff <= max(1, 0.01*ntot) || (ndiff >= ndiffls[end-1] && it >= 4)                
                alwaysaccept = false        
                nstuck = 0                
            end
        end
        # if !upbyhalf && !alwaysaccept  
        #     if nstuck >  0 || ndiff == 0
        #         upbyhalf = true
        #         nstuck = 0                
        #     end
        # end
        if !upbyblock && !alwaysaccept  
            if nstuck >= maxstuck || ndiff == 0
                upbyblock = true                
                nstuck = 0            
            end       
        end    
    end    
    if !isempty(targetls) && isrefineerror
        msg = string("chr=",chrid,", run=", run, ", ")
        msg *= tostring_errorprior(errorprior,targetls)
        printconsole(io,verbose,msg)     
    end
    # calculate logp used for selectin best phasinfg run    
    rand_update_bvpair!(siblogl, bvpair,phaploindex,phaploset, priorspace, priorprocess,
        chroffgeno,chroffformat, chrmultiallele, polyped; 
        phaseerror_midparent, byneighbor, isrand = false, incldelmarker=true,  # isrand affects bvpair and thus logl for deleted marerks
        doseerror = estdoseerror, baseerror = estbaseerror, 
        allelicbias = estallelicbias, overdispersion = estoverdispersion,           
    )  
    logl = sum_margloglls(siblogl)    
    # uniform prior for  bvpair, so its prioir is not included. 
    # logprior_error: incldue deleted markers    
    logprior_error = cal_logprior_error(targetls, estdoseerror, estbaseerror, estallelicbias, estoverdispersion, errorprior)    
    # logl_parentdata: list of logprior for each parent
    logl_parentdata = cal_logl_parentdata(phaploindex, phaploweight)
    if !isempty(fixedparentls)
        logl_parentdata = [in(i,fixedparentls) ? nothing : logl_parentdata[i] for i in eachindex(logl_parentdata)]
    end 

    # set phaplo to missing at deleted markers
    phaplo = get_phaplo(phaploindex,phaploset)    
    pri1=first(values(priorprocess))
    excl=.!(pri1.markerincl)
    ndel = sum(excl)
    if ndel > 0
        phaplo = Matrix{Any}(phaplo)
        phaplo[excl,:] .= missing        
    end        

    phaplo, bvpair,logl, logprior_error, logl_parentdata, ndel, targetls, estdoseerror, estbaseerror, estallelicbias,estoverdispersion      
end

function tostring_errorprior(errorprior::NamedTuple, targetls::AbstractVector)
    isempty(targetls) && return ""
    join([string(i,".prior=", errorprior[Symbol(i)]) for i in targetls],", ")    
end

function update_errorprior(errorprior::NamedTuple; 
    vartargetls::AbstractVector,
    estdoseerror::Union{Real,AbstractVector}, 
    estbaseerror::Union{Real,AbstractVector},  
    estallelicbias::Union{Real,AbstractVector}, 
    estoverdispersion::Union{Real,AbstractVector})
    errnamels = string.(keys(errorprior))
    priorls = [begin 
        prior = errorprior[Symbol(errname)]
        if errname == "doseerror"
            errls = estdoseerror
        elseif errname == "baseerror"
            errls = estbaseerror
        elseif errname == "allelicbias"
            errls = estallelicbias
        elseif errname == "overdispersion"
            errls = estoverdispersion
        else
            @error string("unknow errname=",errname)
        end
        if isa(errls, AbstractVector) && in(errname, vartargetls)
            errmean = mean(errls)
            errvar = var(errls)
            if isa(prior,LogNormal)
                ls = log.(errls)
                prior = LogNormal(round(mean(ls),sigdigits=4),round(std(ls),sigdigits=4))                
            elseif isa(prior, Gamma) || isa(prior, Exponential)                                
                prior = Gamma(round(errmean^2/errvar,sigdigits=6),round(errvar/errmean,sigdigits=6))
            elseif isa(prior, Beta)                    
                prefact = (errmean*(1-errmean)/errvar - 1)
                prior = Beta(round(prefact*errmean,sigdigits=6),round(prefact*(1-errmean),sigdigits=6))
            elseif isa(prior, Uniform)    
                errmin,errmax = extrema(errls)
                prior = Uniform(errmin,errmax)            
            else
                if in(errname, ["doseerror", "baseerror"])                    
                    prefact = (errmean*(1-errmean)/errvar - 1)
                    prior = Beta(round(prefact*errmean,sigdigits=6),round(prefact*(1-errmean),sigdigits=6))
                elseif errname == "allelicbias"
                    ls = log.(errls)
                    prior = LogNormal(round(mean(ls),sigdigits=4),round(std(ls),sigdigits=4))                                    
                elseif errname == "overdispersion"                    
                    prior = Gamma(round(errmean^2/errvar,sigdigits=6),round(errvar/errmean,sigdigits=6))
                else
                    @error string("unknow errname=",errname)
                end
            end
        end
        prior 
    end for errname in errnamels]
    (; (keys(errorprior) .=> priorls)...)
end

function cal_israndls(polyped::PolyPed; minprogeny_perploidy=10)
    parent2off = get_parent2progeny(polyped)
    parentls = polyped.parentinfo[!,:individual]
    ploidyls = polyped.parentinfo[!,:ploidy]
    nls = [length(parent2off[p]) for p in parentls]
    nls .>= (minprogeny_perploidy .* ploidyls)
end

function cal_ndiff_phaploindex(oldphaploindex::AbstractVector, newphaploindex::AbstractVector, 
    phaploset::AbstractVector, markerincl,fixedparentls::AbstractVector)
    phaplo = get_phaplo(oldphaploindex, phaploset)
    newphaplo = get_phaplo(newphaploindex, phaploset)
    nparent = size(phaplo,2)
    if isempty(fixedparentls)
        ndiff = sum(cal_parentgeno_nphaseerr(phaplo[markerincl,:], newphaplo[markerincl,:])) 
    else
        nonfixedls = setdiff(1:nparent, fixedparentls)
        ndiff = sum(cal_parentgeno_nphaseerr(phaplo[markerincl,nonfixedls], newphaplo[markerincl,nonfixedls])) 
    end
    ntot = sum(markerincl) * (nparent -length(fixedparentls))
    ndiff, ntot
end

function cal_parentgeno_nphaseerr(refparentgeno::AbstractMatrix, estparentgeno::AbstractMatrix)
    refhapls = [permutedims(reduce(hcat,i)) for i in eachcol(refparentgeno)]
    esthapls = [permutedims(reduce(hcat,i)) for i in eachcol(estparentgeno)]
    nphaseerr=last(calpermhomolog(refhapls, esthapls)) # each vector element refers to a parent
    nphaseerr
end


function up_phaploindex!(phaploindex::AbstractVector, phaplo::AbstractMatrix, 
    phaploset::AbstractVector, 
    phaploweight::AbstractVector,
    phaploposterior::AbstractVector)
    for j in eachindex(phaploset)
        phaploset_j = phaploset[j]
        phaploweight_j = phaploweight[j]
        phaploposterior_j = phaploposterior[j]
        phaplo_j = phaplo[:,j]
        for i in eachindex(phaploset_j,phaploweight_j)
            hh,ww,pp = phaploset_j[i], phaploweight_j[i], phaploposterior_j[i]
            h = phaplo_j[i]
            hindex = findfirst(==(h),hh)
            if isnothing(hindex)
                push!(hh,h)
                ww .*= 0.99
                push!(ww, 0.01)
                pp .*= 0.99
                push!(pp, 0.01)
                phaploindex[j][i] = length(hh)
            else
                phaploindex[j][i] = hindex
            end
        end                        
    end
    phaploindex
end    

function up_parent_phaploprior!(phaploindex::AbstractVector,
    phaploset::AbstractVector, phaploweight::AbstractVector,    
    phaploposterior::AbstractVector, 
    chrparentgeno::AbstractMatrix, 
    chrparentformat::AbstractVector, 
    chrmultiallele::AbstractVector, 
    polyped::PolyPed; 
    phaseerror::Real, 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector})
    phaplo = get_phaplo(phaploindex, phaploset)
    phaploset2, phaploweight2 = init_phaploprior(chrparentgeno, chrparentformat, chrmultiallele,polyped; 
        phaseerror, doseerror, baseerror, allelicbias, overdispersion, onlyfounder=false)      
    nparent = length(phaploset)
    nmarker = length(first(phaploset))
    for p in 1:nparent
        phaploset2_p = phaploset2[p]
        phaploweight2_p = phaploweight2[p]
        
        phaploset_p = phaploset[p]        
        phaploweight_p = phaploweight[p]        
        phaploposterior_p = phaploposterior[p]
        for t in 1:nmarker
            current = phaploset_p[t]
            new = phaploset2_p[t]
            current == new && continue            
            if issubset(current,new)
                # replace current with new                
                phaploset_p[t] = copy.(phaploset2_p[t])                
                phaploweight_p[t] = phaploweight2_p[t]           
                phaploposterior_p[t] = phaploweight2_p[t]     
            else
                # merge current and new
                hls = union(current, new)
                wls = zeros(length(hls))
                for h in hls
                    if in(h,current)
                        i = findfirst(==(h),current)
                        wls[i] += phaploweight_p[t][i]
                    end
                    if in(h,new)
                        i = findfirst(==(h),new)
                        wls[i] += phaploweight2_p[t][i]
                    end
                end 
                wls ./= sum(wls)
                phaploweight_p[t] = wls
                phaploposterior_p[t] = copy(wls)
                phaploset_p[t] = copy.(hls)                
            end            
        end
    end                                                           
    phaploset .= phaploset2
    phaploweight .= phaploweight2
    phaploposterior .= deepcopy.(phaploweight)
    up_phaploindex!(phaploindex, phaplo, phaploset, phaploweight,phaploposterior)  
    nothing
end


# TODO: for other formats GT, AD, ..
function get_chroffgeno(polygeno::PolyGeno, chr::Integer)
    # [ismissing(i) ? missing : unphased_genocode_2_index(i) for i in polygeno.offspringgeno[chr]]
    res = copy(polygeno.offspringgeno[chr])
    chrformat = polygeno.markermap[chr][!,:offspringformat]
    b = chrformat .== "GT_unphased"
    if any(b)
        subres = view(res, b, :)
        subres .= [ismissing(i) ? missing : unphased_genocode_2_index(i) for i in subres]
    end
    res
end


function index2haplo(phaploset::AbstractVector,haploindex::AbstractVector)
    [reduce(hcat,map((x,y)->ismissing(y) ? missings(length(x[1])) : x[y],
        phaploset[i],haploindex[i]))' for i=eachindex(phaploset)]
end


