
"""
    polyMapRefine!(polygeno::PolyGeno, keyargs...)

performs marker map refinning for polygeno with phased parent genotypes.
    Modifies polygeno.markermap into a refined genetic map.

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`doseerr::Real=0.01`: genotypic error probability for offspring phasing.

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosomes
to be considered, with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

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

`logfile::Union{AbstractString,IOStream}=string(outstem,".log")`: output filenames
or stream for writing log.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyMapRefine!(polygeno::PolyGeno;
    doseerr::Real=0.01,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    byneighbor::Union{Nothing,Bool}=nothing,    
    refineorder::Bool=false,
    maxwinsize::Integer=50,
    inittemperature::Real=4,
    coolingrate::Real=0.5,
    stripdis::Real=20, # centiMorgan
    maxdoseerr::Real=0.5,
    skeletonsize::Integer=50,
    # missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    isplot::Bool=false, 
    logfile::Union{Nothing,AbstractString,IO}="outstem.log",
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    starttime = time()
    if isnothing(logfile)
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
            printpkgst(io,verbose,"PolyOrigin")
        else
            io=logfile
        end
    end
    printconsole(io,verbose,string("polyMapRefine!, logfile=", logfile, ", ", Dates.now()))
    if !isnothing(chrsubset)
        verbose && @info string("chrsubset=",chrsubset)
    end
    if !isnothing(snpsubset)
        verbose && @info string("snpsubset=",snpsubset)
    end
    if isnothing(byneighbor)
        byneighbor = max(polygeno.parentinfo[!,:ploidy]...) >=6
    end
    msg = string("list of option values: \n",
        "doseerr = ", doseerr, "\n",
        "chrpairing = ",chrpairing, "\n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",
        "isparallel = ",isparallel,"\n",
        "byneighbor = ",byneighbor,"\n",        
        "maxwinsize = ",maxwinsize,"\n",
        "inittemperature = ",inittemperature, "\n",
        "coolingrate = ", coolingrate, "\n",
        "stripdis = ", stripdis, "\n",
        "maxdoseerr = ", maxdoseerr, "\n",
        "skeletonsize = ", skeletonsize, "\n",
        # "missingstring = ",missingstring,"\n",
        "outstem = ",outstem,"\n",
        "logfile = ",io,"\n",
        "workdir = ",workdir,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
    getsubPolyGeno!(polygeno,chrsubset=chrsubset,snpsubset=snpsubset)
    inputmap = deepcopy(polygeno.markermap)
    isdel=maprefine_allchr!(polygeno,doseerr,chrpairing,  
        stripdis,maxdoseerr, skeletonsize,refineorder,maxwinsize,
        inittemperature,coolingrate, isparallel,byneighbor,io,verbose)
    delmarker = setparentphase!(polygeno,isdel)
    if size(delmarker,1)>0
        msg = string("delete ", size(delmarker,1)," markers: ",
            join(sort(string.(delmarker[!,:marker])),", "))
        printconsole(io,verbose,msg)
    end
    if !isnothing(outstem)
        outfile = string(outstem,"_maprefined.csv")
        savegenodata(outfile,polygeno,workdir=workdir)
        msg = string("maprefined file: ", outfile)
        printconsole(io,verbose,msg)
    end
    tau = round.(calmapkendall(inputmap,polygeno.markermap),digits=4)
    msg = string("Kendall tau between inputmap and refinedmap = ", tau)
    printconsole(io,verbose,msg)
    # plot 
    if !isnothing(outstem) && isplot 
        try 
            plotMapComp(inputmap, polygeno.markermap;
                xlabel = "Input map position(cM)", 
                ylabel = "Refined map position(cM)",     
            )
            figdir = joinpath(workdir, outstem * "_plots")
            isdir(figdir) || mkdir(figdir)       
            outfile = string(outstem,"_mapcompare.png")        
            savefig(getabsfile(figdir, outfile))
            msg = string("refinedmap plot: ", outfile)
            printconsole(io,verbose,msg)
        catch err
            @warn string("Failed to visualize the comparision between input and refined maps")
            @error err
        end
    end
    printconsole(io,verbose,string("End, ", Dates.now(),", time used = ",
        round(time()-starttime), " seconds by polyMapRefine!"))
    if typeof(logfile) <: AbstractString
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polygeno
end


function maprefine_allchr!(inputpolygeno::PolyGeno,
    doseerr::Real,chrpairing::Integer, 
    stripdis::Real,maxdoseerr::Real,skeletonsize::Integer,
    refineorder::Bool,maxwinsize::Integer,inittemperature::Real,coolingrate::Real,
    isparallel::Bool,byneighbor::Bool,
    io::Union{Nothing,IO},verbose::Bool)
    # exclude outlier
    outlier0=inputpolygeno.offspringinfo[!,:isoutlier]
    outlier = skipmissing(outlier0)
    isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
    noutlier == 0 && printconsole(io,verbose,string("no offspring excluded"))
    if noutlier > 0
        offoutid = inputpolygeno.offspringinfo[outlier0,:individual]
        msg = string("exclude ", noutlier, " outlier offspring: ",join(offoutid,", "))
        printconsole(io,verbose,msg)
        polygeno = deepcopy(inputpolygeno)
        delOutlier!(polygeno)
    else
        polygeno = inputpolygeno
    end
    maxiter=30
    nchr=length(polygeno.markermap)
    res=Vector{Vector}(undef,nchr)
    if isparallel && nprocs()>1
        polygenols=[getsubPolyGeno(polygeno,chrsubset=[chr]) for chr=1:nchr]
        resmap = pmap(x->maprefine_chr(x,1,doseerr,chrpairing,
            byneighbor, stripdis, maxdoseerr, skeletonsize,refineorder,
            maxwinsize,maxiter,inittemperature,coolingrate,nothing,verbose),polygenols)
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
            res[chr] = first(maprefine_chr(polygeno,chr,doseerr,chrpairing,
                byneighbor, stripdis, maxdoseerr, skeletonsize,refineorder,
                maxwinsize,maxiter,inittemperature,coolingrate,io,verbose))
        end
    end
    if noutlier > 0
        polygeno = inputpolygeno
    end
    isdel=Vector{Vector}(undef,nchr)
    for chr=1:nchr
        snporder, isdel[chr], snpid, snppos= res[chr]
        polygeno.markermap[chr][!,:marker] = snpid
        polygeno.markermap[chr][!,:position] = snppos
        polygeno.parentgeno[chr] = polygeno.parentgeno[chr][snporder,:]
        polygeno.offspringgeno[chr] = polygeno.offspringgeno[chr][snporder,:]
    end
    isdel
end

function maprefine_chr(polygeno::PolyGeno, chr::Integer,
    doseerr::Real,chrpairing::Integer,
    byneighbor::Bool,    
    stripdis::Real,maxdoseerr::Real,skeletonsize::Integer,
    refineorder::Bool,maxwinsize::Integer,maxiter::Integer,
    inittemperature::Real,coolingrate::Real,
    io::Union{IO,Nothing},verbose::Bool)
    starttime = time()
    isnothing(io) && (io=IOBuffer(append=true))
    priorspace = getpriorstatespace(polygeno,chrpairing)
    idlen=max([length(i[1,:chromosome]) for i=polygeno.markermap]...)
    chrid=lpad(polygeno.markermap[chr][1,:chromosome],idlen)
    fhaplo= getfhaplo(polygeno.parentgeno[chr])
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    chrdose = polygeno.offspringgeno[chr]
    nsnp = size(chrdose,1)
    snporder=collect(1:nsnp)
    doseerrls0=[doseerr for i=1:nsnp]
    doseerrls=Vector{Union{Missing,Float64}}(doseerrls0)
    maxvalent = max(digits(chrpairing)...)
    priorprocess = getpriorprocess(polygeno,chr,maxvalent)
    bvpair = randinitbvpair(priorspace,polygeno)
    bvpair,logl = refinebvpair!(bvpair, chr,doseerrls,chrdose,priorspace,priorprocess,
        polygeno;snporder,byneighbor)
    temperature = inittemperature
    nstuck=0.0
    winsize=maxwinsize
    for it=1:maxiter
        startt = time()
        # update ordering
        if refineorder
            snporder,priorprocess,logbook = updateorder!(snporder,bvpair,
                priorspace,priorprocess, chrdose, deriveddose, doseerrls, polygeno;
                reversechr=false, maxwinsize,temperature)
            newwinsize = length(logbook)+1
            if newwinsize==winsize
                nstuck += 1
            else
                winsize=newwinsize
                if temperature ≈ 0.0
                    nstuck += 0.5  # max 6 iterations at T= 0 
                else
                    nstuck = 0
                end                
            end
        else
            logbook = []
            nstuck += 1 # 3 iterations when refine only distance, since stop when nstuck>=3
        end
        # update inter-marker distances        
        logldis = updatedistance!(snporder,bvpair,priorspace,priorprocess,
            chrdose, deriveddose,doseerrls, polygeno)
        # logbook[end][2]
        push!(logbook,["#update_distance",logldis,round(logldis-logl,digits=2),1])
    
        # update bvpair
        bvpair,loglbv = refinebvpair!(bvpair, chr,doseerrls,chrdose,priorspace,priorprocess,
            polygeno;snporder,byneighbor)        
        push!(logbook,["#update_valent",loglbv,round(loglbv-logldis,digits=2),1.0])
        # update doseerr
        doseerrls,logleps = updatedoseerrls!(doseerrls, snporder,bvpair,chrdose,deriveddose,
            priorspace, priorprocess,polygeno)
        delsnpeps!(priorprocess,doseerrls,snporder,maxdoseerr=maxdoseerr)
        push!(logbook,["#update_doseerr",logleps,round(logleps-loglbv,digits=2),1.0])
        
        #marker deletion; priorprocess is modified
        # polymarkerdel!(dataprobset,bvpair,priorspace, priorprocess,polygeno,
        #     snporder=snporder,chrindex=chr)
        # logllist = calmarglogl(dataprobset, priorspace,priorprocess,polygeno,
        #     [[i] for i=bvpair],snporder=snporder)
        # logldel=sum(vcat(logllist...))
        # push!(logbook,["#update_deletion",logldel,round(logldel-loglbv,digits=2),1.0])
        # # print
        # colid = Symbol.([string("#refine_",chrid),"logl","delt_logl","accept"])
        # dflogbook = DataFrame(permutedims(reduce(hcat,logbook)),colid)
        pri1=first(values(priorprocess))
        chrlen = round(100*sum(pri1.markerdeltd),digits=1)
        logl =round(loglbv,digits=1)
        ndel = sum(first(values(priorprocess)).markerincl[1:end-1].== 0)
        msgmore = refineorder ? string(", win=",winsize,", T=",round(temperature,digits=4)) : ""
        msg =string("chr=",chrid,
            # nprocs()>1 ? string(", p=",myid()) : "",
            ", it=",it,
            ", logl=",logl,
            ", len=",chrlen,
            ", ndel=",ndel,
            ", <err>=",round(mean(skipmissing(doseerrls)),digits=4),
             msgmore,
            #  ", stuck=",nstuck,
             ", tused=",round(time()-startt,digits=1),"s",
            )
        printconsole(io,verbose,msg)
        # write(io, string("epsls: ", join(round.(doseerrls,digits=3),","), "\n"))
        # verbose && println(dflogbook)
        # CSV.write(io, dflogbook,header=true,append=true)
        if  ((temperature ≈ 0.0) && (winsize==2)) || nstuck>=3 || it==maxiter
            temperature = 0.0
            stripchrend!(priorprocess,stripdis=stripdis)
            # polymarkerdel!(bvpair,priorspace, priorprocess,polygeno,
            #     caldistance=true,snporder=snporder,chrindex=chr)
            for itdis=1:2
                startt = time()
                it+=1
                logldis = updatedistance!(snporder,bvpair,priorspace,priorprocess,
                    chrdose, deriveddose,doseerrls, polygeno)
                stripchrend!(priorprocess,stripdis=stripdis)
                pri1=first(values(priorprocess))
                chrlen = round(100*sum(pri1.markerdeltd),digits=1)
                logl =round(sum(logldis),digits=1)
                ndel = sum(first(values(priorprocess)).markerincl[1:end-1].== 0)
                msg =string("chr=",chrid,
                        # nprocs()>1 ? string(", p=",myid()) : "",
                        ", it=",it,", logl=",logl,", len=",chrlen,
                        ", ndel=",ndel,
                        ", tused=",round(time()-startt,digits=1),"s",
                        )
                printconsole(io,verbose,msg)
            end
            if sum(first(values(priorprocess)).markerincl)>skeletonsize
                # doseerrls and dataprobset will be modified
                # cd("C:\\Chaozhi\\Workspace\\JuliaWorkspace\\Workspace_Polyploid\\PolyOrigin\\examples\\dihaploid")
                # serialize("temp.test",[doseerrls,dataprobset,snporder,
                #     bvpair,priorspace,priorprocess,polygeno,
                #     chr,chrid,skeletonsize,verbose,io])
                skeletonprior =  getskeletonmap!(doseerrls,snporder,
                    bvpair,priorspace,priorprocess,polygeno,
                    chr,chrid,skeletonsize,verbose,io)
                rescalemap!(priorprocess,skeletonprior)
            end
            break
        end
        if temperature > 0.5
            temperature *= coolingrate
        elseif 0.05< temperature <= 0.5
            temperature *= coolingrate^2
        elseif temperature<=0.05
            temperature =0.0 
        end
    end
    pri1=first(values(priorprocess))
    chrlen = round(100*sum(pri1.markerdeltd),digits=1)
    nmarker = sum(first(values(priorprocess)).markerincl)
    msg = string("chr=",chrid,
        # nprocs()>1 ? string(", p=",myid()) : "",
        ", done, elapsed=",round(time()-starttime), "seconds",
        ", len=",chrlen, "cM", ", #marker=",nmarker)
    printconsole(io,verbose,msg)
    markerid=polygeno.markermap[chr][!,:marker]
    checksnporder(snporder,priorprocess,markerid) || @error("inconsistent snporder and markerid!")
    checktranprobseq(priorprocess) || @error("inconsistent tranprobseq and deltd!")
    pri1=first(values(priorprocess))
    excl = .!(pri1.markerincl)
    issnpdel = excl[snporder]
    snpid = pri1.markerid
    snppos = round.(vcat([0],100*accumulate(+,pri1.markerdeltd[1:end-1])),digits=2)
    [snporder, issnpdel, snpid, snppos], io
end

function refinebvpair!(bvpair::AbstractVector, chr::Integer,
    doseerrls::AbstractVector,chrdose::AbstractMatrix,
    priorspace::AbstractDict, priorprocess::AbstractDict,
    polygeno::PolyGeno;
    snporder::Union{Nothing,AbstractVector}=nothing,
    byneighbor::Bool=false)
    fhaploset=[[[i] for i=j] for j=eachcol(polygeno.parentgeno[chr])]
    fhaploindex=[ones(Int, length(j)) for j=fhaploset]
    siblogl = repeat([-Inf],length(bvpair))
    popidls = keys(priorspace)
    if byneighbor
        updatebvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,
            priorspace,priorprocess,polygeno,popidls;snporder)
    else
        randbvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,
            priorspace,priorprocess,polygeno,popidls; isrand=false,snporder)
    end
    bvpair, sum(siblogl)
end

# function randbvpair(snporder::AbstractVector,
#     priorspace::AbstractDict,priorprocess::AbstractDict,
#     chrdose::AbstractMatrix, deriveddose::AbstractDict,
#     doseerrls::AbstractVector, polygeno::PolyGeno)
#     bvpairprop = getbvpairprop(priorspace,polygeno)
#     logllist = calmarglogl(doseerrls,deriveddose,chrdose,priorspace,
#         priorprocess,polygeno,bvpairprop; snporder=snporder)
#     bvpair=[rand(findall(i .== max(i...))) for i = logllist]
#     logl = sum(map((x,y)->x[y], logllist,bvpair))
#     bvpair,logl
# end

function gettrandict(d::Real,priorprocess::AbstractDict)
    Dict([strkey=>getgametetran(d,pri.nvalent) for (strkey, pri) in priorprocess])
end

function getbvkeyls(bvpair::AbstractVector,priorspace::AbstractDict,
    polygeno::PolyGeno)
    bvkeyls = Vector{String}(undef,length(bvpair))
    popidls = polygeno.designinfo[!,:population]
    for popid=popidls
        keyls = priorspace[popid]["valentkey"]
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        bvkeyls[offls] =keyls[bvpair[offls]]
    end
    bvkeyls
end

function callogldis(d::Real,priorprocess::AbstractDict,
    bvkeyls::AbstractVector,dataprobls::AbstractVector,
    fwprobls::AbstractVector,bwprobls::AbstractVector)
    tranprobdict=gettrandict(d,priorprocess)
    tranprobls=[[get(tranprobdict,j,missing)' for j=split(i,"|")] for i=bvkeyls]
    noff = length(bvkeyls)
    sum(log(sum(kronvec(tranprobls[off]...,fwprobls[off]) .* dataprobls[off] .* bwprobls[off])) for off=1:noff)
end

function calinterdis(tnow::Integer,tnext::Integer,priorprocess::AbstractDict,
    bvkeyls::AbstractVector,dataprobls::AbstractVector,
    fwdict::AbstractDict,logbw::AbstractVector)
    # Haldane map function, here rfreq defined one chromsome pairing
    # so that the map function is the same as the diploid case
    # priorlength=min(0.1,max(0.01,20/length(markerdeltd)))
    pri1=first(values(priorprocess))
    markerdeltd = pri1.markerdeltd[1:end-1]
    # priorlength= 3*quantile(markerdeltd,0.95)
    accuracygoal, precisiongoal, itmax = 2, 2, 100
    lowbound,upbound = log(10^(-7)), log(10.0)
    fwprobls = [i[tnow,:] for i=fwdict["fwprob"]]
    bwprobls = [begin
        bw0 = logbw[off][tnext,:]
        exp.(bw0 .- max(bw0...))
    end for off=1:length(logbw)]
    # the term -exp(x) is to set x->0 for non-informative data.
    function f(x::Real)
        callogldis(exp(x),priorprocess,bvkeyls,dataprobls,fwprobls,bwprobls)-exp(x)
    end
    # b=sum([max(logbw[off][tnext,:]...) for off=1:noff])
    # c=sum([i[1] for i=fwdict["fwlogl"]])
    xold = log(max(exp(lowbound),markerdeltd[tnow]))
    res=brentMax(f,lowbound,upbound,xold,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)
    tnowdis = exp(res[1])
    tnowdis<=10^(-6.) && (tnowdis=0.0)
    logl = res[2]
    tnowdis, logl
end

function updatedistance!(snporder::AbstractVector,
    bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno;
    reversechr::Bool=false)
    if reversechr
        reverse!(snporder)
        reverse!(doseerrls)
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    bvkeyls = getbvkeyls(bvpair,priorspace,polygeno)
    pri1=first(values(priorprocess))
    tseq = findall(pri1.markerincl)
    logbw = callogbackward(chrdose, deriveddose, bvpair,priorspace,priorprocess,
        doseerrls, polygeno; snporder=snporder)
    snp = snporder[tseq[1]]
    dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
        chrdose,priorspace,polygeno)
    fwdict = calinitforward(tseq[1],dataprobls,bvpair,priorspace,priorprocess,
        polygeno,snporder=snporder)
    kkmax = length(tseq)-1
    for kk=1:kkmax
        snp = snporder[tseq[kk+1]]
        dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
            chrdose,priorspace,polygeno)
        tnowdis,logl= calinterdis(tseq[kk],tseq[kk+1],priorprocess,bvkeyls,
            dataprobls,fwdict,logbw)
        setdistanceat!(priorprocess,tseq[kk],tnowdis)
        calnextforward!(fwdict,tseq[kk],tseq[kk+1],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
    end
    logl = sum(calindlogl(fwdict,logbw, tseq[end]))
    if reversechr
        reverse!(snporder)
        reverse!(doseerrls)
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    logl
end

function updatedoseerrls!(doseerrls::AbstractVector,snporder::AbstractVector,
    bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,chr::Integer, polygeno::PolyGeno)
    chrdose=polygeno.offspringgeno[chr]
    fhaplo= getfhaplo(polygeno.parentgeno[chr])
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    updatedoseerrls!(doseerrls, snporder,bvpair,chrdose,deriveddose,
        priorspace, priorprocess,polygeno)
end

function updatedoseerrls!(doseerrls::AbstractVector,snporder::AbstractVector,
    bvpair::AbstractVector,chrdose::AbstractMatrix,deriveddose::AbstractDict,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno)
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)
    logbw = callogbackward(chrdose, deriveddose, bvpair,priorspace,priorprocess,
        doseerrls, polygeno; snporder=snporder)
    kk=1
    snp = snporder[tseq[kk]]
    doseerrls[snp]=first(callogldoseerrls(0,tseq[kk],doseerrls[snp],snporder,nothing,
        logbw,bvpair,deriveddose,chrdose,priorspace,priorprocess,polygeno))
    dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
        chrdose,priorspace,polygeno)
    fwdict = calinitforward(tseq[kk],dataprobls,bvpair,priorspace,priorprocess,
        polygeno,snporder=snporder)
    for kk=2:length(tseq)
        snp = snporder[tseq[kk]]
        # println("kk,t,snp,doseerr=", [kk,tseq[kk],snp,doseerrls[snp]])
        doseerrls[snp]=first(callogldoseerrls(tseq[kk-1],tseq[kk],doseerrls[snp],snporder,fwdict,
            logbw,bvpair,deriveddose,chrdose,priorspace,priorprocess,polygeno))
        dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
            chrdose,priorspace,polygeno)
        calnextforward!(fwdict,tseq[kk-1],tseq[kk],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
    end
    logl = sum(calindlogl(fwdict,logbw, tseq[end]))
    doseerrls,logl
end

function callogldoseerrls(tpre::Integer,tnow::Integer,doseerr::Real,
    snporder::AbstractVector,fwdict::Union{Nothing,AbstractDict},
    logbw::AbstractVector,
    bvpair::AbstractVector,deriveddose::AbstractDict,chrdose::AbstractMatrix,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno)
    accuracygoal, precisiongoal, itmax = 3, 3, 50
    lowbound,upbound = log(10^(-6.)), log(1.0-10^(-6.))
    # lowbound,upbound =10^(-3.), 1-10^(-3.)
    if tpre==0
        bvkeyls = getbvkeyls(bvpair,priorspace,polygeno)
        fwprobls = getstartprobls(priorprocess,bvkeyls)
    else
        fwprobls = [i[tpre,:] for i=fwdict["fwprob"]]
    end
    bwprobls = [begin
        bw0 = logbw[off][tnow,:]
        exp.(bw0 .- max(bw0...))
    end for off=1:length(logbw)]
    snp = snporder[tnow]
    # +x-exp(x)
    function f(x::Real)
        callogldoseerr(snp,exp(x),fwprobls,bwprobls,bvpair,deriveddose,chrdose,
            priorspace,polygeno)
    end
    xold = log(max(exp(lowbound),doseerr))
    res=brentMax(f,lowbound,upbound,xold,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)
    neweps = exp(res[1])
    logl = res[2]
    neweps, logl
end

function callogldoseerr(snp::Integer,doseerr::Real,
    fwprobls::AbstractVector,bwprobls::AbstractVector,
    bvpair::AbstractVector,deriveddose::AbstractDict,chrdose::AbstractMatrix,
    priorspace::AbstractDict,polygeno::PolyGeno)
    dataprobls = calsitedataprobls(snp,doseerr,bvpair,deriveddose,chrdose,
        priorspace,polygeno)
    noff = length(bvpair)
    sum([log(sum(fwprobls[off] .* dataprobls[off] .* bwprobls[off])) for off=1:noff])
end


function getsegforward(isreverse::Bool,tseq::AbstractVector, kkmin::Integer,
    kkmax::Integer,fwdict::AbstractDict,
    snporder::AbstractVector, bvpair::AbstractVector,
    priorspace::AbstractDict, priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,
    doseerrls::AbstractVector,polygeno::PolyGeno)
    bvkeyls = getbvkeyls(bvpair,priorspace,polygeno)
    if kkmin==1
        startprobls = getstartprobls(priorprocess,bvkeyls)
        leftlogl  = zeros(length(bvkeyls))
    else
        startprobls = [i[tseq[kkmin-1],:] for i=fwdict["fwprob"]]
        trandict = Dict([strkey=>(pri.tranprobseq[tseq[kkmin-1]])
            for (strkey, pri) in priorprocess])
        startprobls=[begin
                strkey = bvkeyls[i]
                tranprob=[get(trandict,i,missing)' for i=split(strkey,"|")]
                kronvec(tranprob..., startprobls[i])
            end for i=1:length(bvkeyls)]
        leftlogl = [i[tseq[kkmin-1]] for i=fwdict["fwlogl"]]
    end
    tseqseg = tseq[kkmin:kkmax]
    tranprobdict = Dict([strkey=> pri.tranprobseq[tseqseg[1:end-1]]  for (strkey, pri) in priorprocess])
    isreverse && (tranprobdict =Dict([strkey => reverse(val) for (strkey, val) in tranprobdict]))
    tranprobdict2 = Dict([begin
        tran = [get(tranprobdict,j,missing) for j=split(i,"|")]
        tran2 = map((x,y)->kron(x,y),tran...)
        i=>tran2
    end for i=unique(bvkeyls)])
    noff = length(bvkeyls)
    fwprob=Vector{Matrix{Float64}}(undef,noff)
    fwlogl=Vector{Vector{Float64}}(undef,noff)
    popidls = polygeno.designinfo[!,:population]
    tseqseg2 = isreverse ? reverse(tseqseg) : tseqseg
    snps = snporder[tseqseg2]
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        ploidy = polygeno.offspringinfo[first(offls),:ploidy]
        condstates = priorspace[popid]["condstate"]
        for off=offls
            bv = bvpair[off]
            deriveddose2 = deriveddose[popid][snps,condstates[bv]]
            dataprobseq = caldataprob(chrdose[snps,off],ploidy,deriveddose2,doseerrls[snps])
            startprob = Vector{Float64}(startprobls[off])
            tranprobseq = get(tranprobdict2,bvkeyls[off],missing)
            fwseg=forward(startprob,tranprobseq,dataprobseq)
            fwprob[off]=reduce(hcat,fwseg.fwprob)'
            fwlogl[off] = accumulate(+,log.(fwseg.fwscale)) .+ leftlogl[off]
        end
    end
    Dict("fwprob"=>fwprob,"fwlogl"=>fwlogl)
end

function calsegbwlogl(fwdict::AbstractDict,logbw::AbstractVector,t::Integer)
    fwprob=fwdict["fwprob"]
    fwlogl=fwdict["fwlogl"]
    loglls=[begin
        lp=logbw[i][t,:]
        lpmax=max(lp...)
        log(dot(exp.(lp .- lpmax), fwprob[i][end,:])) + lpmax + fwlogl[i][end]
    end for i=1:length(logbw)]
    sum(loglls)
end


function mergefwdict!(fwdict::AbstractDict,fwseg::AbstractDict,
    tseg::AbstractVector)
    fwprob = fwdict["fwprob"]
    fwlogl= fwdict["fwlogl"]
    noff = length(fwprob)
    for i=1:noff
        fwprob[i][tseg,:] = fwseg["fwprob"][i]
        fwlogl[i][tseg] = fwseg["fwlogl"][i]
    end
    fwdict
end

function reverseorder!(winsize::Integer, temperature::Real,
    snporder::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno)
    logbw = callogbackward(chrdose, deriveddose, bvpair,priorspace,priorprocess,
        doseerrls, polygeno; snporder=snporder)
    pri1 = first(values(priorprocess))
    tseq = findall(pri1.markerincl)
    snp = snporder[tseq[1]]
    dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
        chrdose,priorspace,polygeno)
    fwdict = calinitforward(tseq[1],dataprobls,bvpair,priorspace,priorprocess,
        polygeno,snporder=snporder)
    for kk=1:min(winsize-1,length(tseq))
        snp = snporder[tseq[kk+1]]
        dataprobls = calsitedataprobls(snp,doseerrls[snp],bvpair,deriveddose,
            chrdose,priorspace,polygeno)
        calnextforward!(fwdict,tseq[kk],tseq[kk+1],dataprobls,bvpair,priorspace,
            priorprocess,polygeno,snporder=snporder)
    end
    logl = sum(calindlogl(fwdict,logbw, tseq[1]))
    loglhis = Vector{Float64}()
    push!(loglhis,logl)
    accepthis = Vector{Bool}()
    for kk=1:length(tseq)-1
        kkmin=kk
        kkmax = min(kk+winsize-1,length(tseq))
        fwseg = getsegforward(false,tseq,kkmax,kkmax,fwdict,snporder,bvpair,
            priorspace,priorprocess,chrdose, deriveddose,doseerrls,polygeno)
        fwsegprop = getsegforward(true,tseq,kkmin,kkmax,fwdict,snporder,bvpair,
            priorspace,priorprocess,chrdose, deriveddose,doseerrls,polygeno)
        # logl = calsegbwlogl(fwseg,logbw, tseq[kkmax])
        proplogl = calsegbwlogl(fwsegprop,logbw, tseq[kkmax])
        # println("kk=",kk,",logl=",logl, "; proplogl=",proplogl)
        tseg = tseq[kkmin:kkmax]
        if temperature <=0
            isaccept = proplogl >= logl
        else
            isaccept = rand() < min(1,exp((proplogl - logl)/temperature))
        end
        if isaccept
            mergefwdict!(fwdict,fwsegprop,tseg)
            snporder[tseg] = snporder[reverse(tseg)]
            setsegrev!(priorprocess,tseg)
            logl = proplogl
        else
            mergefwdict!(fwdict,fwseg,[tseg[end]])
        end
        push!(accepthis,isaccept)
        push!(loglhis,logl)
    end
    snporder, priorprocess, loglhis, accepthis
end

function updateorder!(snporder::AbstractVector,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,
    doseerrls::AbstractVector, polygeno::PolyGeno;
    reversechr::Bool=false,
    maxwinsize::Integer=50,
    temperature::Real=0)
    if reversechr
        reverse!(snporder)
        reverse!(doseerrls)
        # chrdose=reverse(chrdose, dims=1)
        # deriveddose= Dict([key => reverse(val, dims=1) for (key, val) in deriveddose])
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    winsize=2
    nsnp = length(snporder)
    logbook=Vector()
    while true
        # println("winsize=",winsize)
        snporder,priorprocess,loglhis,accepthis = reverseorder!(winsize,temperature,
            snporder,bvpair,priorspace,priorprocess,
            chrdose, deriveddose, doseerrls, polygeno)
        accept = mean(accepthis[1:end-winsize+1])
        deltlogl = loglhis[end]-loglhis[1]
        temp = vcat(round.([loglhis[end], deltlogl],digits=2),[round(accept,digits=4)])
        push!(logbook,vcat([string("#reverse_w",winsize)],temp))
        if (abs(deltlogl)>1 && winsize<min(maxwinsize,nsnp/2))
            winsize += 1
        else
            break
        end
    end
    if reversechr
        reverse!(snporder)
        reverse!(doseerrls)
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    snporder,priorprocess,logbook
end

function checksnporder(snporder::AbstractVector,priorprocess::AbstractDict,
    markerid::AbstractVector)
    pri1=first(values(priorprocess))
    markerid[snporder] == pri1.markerid
end

function checktranprobseq(priorprocess::AbstractDict)
    pri1=first(values(priorprocess))
    incl = pri1.markerincl
    deltd = pri1.markerdeltd[incl][1:end-1]
    transeq=[getgametetran(i,pri1.nvalent) for i=deltd]
    transeq ≈ pri1.tranprobseq[incl][1:end-1]
end

# function offcorrect!(chr::Integer,snporder::AbstractVector,dataprobset::AbstractVector,
#     bvpair::AbstractVector,priorspace::AbstractDict,
#     priorprocess::AbstractDict,polygeno::PolyGeno;
#     callthreshold::Real=0.95)
#     # genoprob[off] is in the originas marker ordering
#     logl, genoprob =chrposteriordecode(snporder,dataprobset,bvpair,priorspace,
#         priorprocess,polygeno)
#     doseprob = calDoseprob([genoprob],polygeno.parentgeno[[chr]],
#         polygeno.offspringinfo,polygeno.designinfo,priorspace)
#     offkind = kindofgeno(polygeno.offspringgeno)
#     if offkind == "dosage"
#         calldose = [begin
#                 val,index=findmax(i)
#                 val> callthreshold ? index-1 : missing
#             end for i=doseprob[1]]
#         offerr = polygeno.offspringgeno[chr] .== calldose
#         offerr = offerr .=== false
#         @info string("#correction=",sum(offerr))
#         polygeno.offspringgeno[chr][offerr] = calldose[offerr]
#     else
#         # offkind == "probability"
#         @error "todo"
#     end
# end

function updatedscale(dataprobset::AbstractVector,bvpair::AbstractVector,
    snporder::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    # update dscale
    # priorprocess,loglscale,dscale = updatedscale(dataprobset,bvpair,snporder,
    #     priorspace,priorprocess,polygeno)
    # push!(logbook,["#update_dscale",loglscale,round(loglscale-logleps,digits=2),1.0])
    accuracygoal, precisiongoal, itmax = 2, 2, 20
    lowbound,upbound = log(1/10), log(10)
    # pri1=first(values(priorprocess))
    # chrlen = sum(pri1.markerdeltd)
     # - exp(x)*chrlen/5
    f(x)=logldscale(exp(x),dataprobset,bvpair,snporder,priorspace,
        priorprocess,polygeno)
    newx, loglscale, his=brentMax(f,lowbound,upbound,0,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)
    dscale=exp(newx)
    priorprocess=changedscale(priorprocess,dscale)
    priorprocess, loglscale, dscale
end

function logldscale(dscale::Real,dataprobset::AbstractVector,
    bvpair::AbstractVector,snporder::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    bvpairprop = [[i] for i=bvpair]
    priorprocess2= dscale==0 ? priorprocess : changedscale(priorprocess,dscale)
    logllist = calmarglogl(dataprobset, priorspace,priorprocess2,polygeno,
        bvpairprop,snporder=snporder)
    sum(vcat(logllist...))
end

function changedscale(inputpriorprocess::AbstractDict,dscale::Real)
    priorprocess = deepcopy(inputpriorprocess)
    for (strkey, pri) in priorprocess
        pri.markerdeltd .*= dscale
        tls = findall(pri.markerincl)
        pri.markerincl[end]==1 && pop!(tls)
        pri.tranprobseq[tls] = [getgametetran(pri.markerdeltd[t],
            pri.nvalent) for t=tls]
    end
    priorprocess
end

function stripchrend!(priorprocess;stripdis::Real=20)
    pri1=first(values(priorprocess))
    deltd = pri1.markerdeltd[1:end-1]
    nend = max(1,round(Int,0.05length(deltd)))
    ii=findfirst(x->x>stripdis/100.0,deltd[1:nend])
    if !isnothing(ii)
        for (strkey,pri) in priorprocess
            pri.markerincl[1:ii] .= 0
            pri.markerdeltd[1:ii] .= 0
        end
    end
    ii=findfirst(x->x>stripdis/100.0,reverse(deltd[end-nend:end]))
    if !isnothing(ii)
        for (strkey,pri) in priorprocess
            pri.markerincl[end-ii+1:end] .= 0
            pri.markerdeltd[end-ii:end] .= 0
        end
    end
end

function delsnpeps!(priorprocess::AbstractDict,doseerrls::AbstractVector,
    snporder::AbstractVector;maxdoseerr::Real=0.2)
    markerincl=first(values(priorprocess)).markerincl
    ls = findall((doseerrls .> maxdoseerr) .=== true)
    # doseerrls gives the eps for each original snps
    if !isempty(ls)
        doseerrls[ls] .= missing
        ii=[findfirst(x->x==i,snporder) for i=ls]
        markerincl[ii] .= false
        for (key, val) in priorprocess
             setmarkerincl!(val,markerincl)
        end
    end
end

function getskeleton(doseerrls::AbstractVector,snporder::AbstractVector,
    priorprocess::AbstractDict,skeletonsize::Integer)
    ttdoseerrls=doseerrls[snporder]
    pri1=first(values(priorprocess))
    excl = findall(.!pri1.markerincl)
    deltd=pri1.markerdeltd[1:end-1]
    loc =accumulate(+,deltd)
    pushfirst!(loc,0)
    bins = splitindex(loc)
    binls = [setdiff(collect(i),excl) for i=bins]
    binls = binls[length.(binls) .> 0]
    skeleton=[length(bin)>1 ? bin[argmin(ttdoseerrls[bin])] : bin[1] for bin=binls]
    mid=skeleton[2:end-1]
    s=sortperm(ttdoseerrls[mid])
    length(s) > skeletonsize-2 && (s=s[1:skeletonsize-2])
    skeleton2=vcat(skeleton[1],sort(mid[s]),skeleton[end])
    incl=findall(pri1.markerincl)
    incl[1] in skeleton2 || pushfirst!(skeleton2,incl[1])
    incl[end] in skeleton2 || push!(skeleton2,incl[end])
    skeleton2
end


function getskeleton_seg(doseerrls::AbstractVector,snporder::AbstractVector,
    priorprocess::AbstractDict,skeletonsize::Integer)
    ttdoseerrls=doseerrls[snporder]
    pri1=first(values(priorprocess))
    excl = findall(.!pri1.markerincl)
    deltd=pri1.markerdeltd
    loc =accumulate(+,deltd[1:end-1])
    pushfirst!(loc,0)
    bins = splitindex(loc)
    binls = [setdiff(collect(i),excl) for i=bins]
    binls = binls[length.(binls) .> 0]
    skeleton=[bin[argmin(ttdoseerrls[bin])] for bin=binls]
    if length(skeleton)<=skeletonsize
        incl=findall(pri1.markerincl)
        incl[1] in skeleton || pushfirst!(skeleton,incl[1])
        incl[end] in skeleton || push!(skeleton,incl[end])
        return skeleton
    end
    skeldeltd=deltd[[bin[end] for bin=binls]]
    skeldeltd[end]==0 || @error "unexpected last skeleton deltd"
    sum(skeldeltd) ≈ sum(deltd) || @error "unexpected sum skeleton deltd"
    skelloc=accumulate(+,circshift(skeldeltd,1))
    gridsize=skelloc[end]/(skeletonsize-2)
    grid=gridpartition(skelloc,gridsize)
    grid=grid[length.(grid) .> 0]
    grid=[skeleton[i] for i=grid]
    skeleton2=[bin[argmin(ttdoseerrls[bin])] for bin=grid]
    incl=findall(pri1.markerincl)
    incl[1] in skeleton2 || pushfirst!(skeleton2,incl[1])
    incl[end] in skeleton2 || push!(skeleton2,incl[end])
    ndiff = skeletonsize-length(skeleton2)
    if ndiff>0
        ii=setdiff(skeleton,skeleton2)
        s=sortperm(ttdoseerrls[ii])
        length(s) > ndiff && (s=s[1:ndiff])
        skeleton2=sort(vcat(skeleton2,ii[s]))
    end
    skeleton2
end


function getskeletonmap!(doseerrls::AbstractVector,snporder::AbstractVector,
    bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,
    chr::Integer,chrid::AbstractString,
    skeletonsize::Integer,verbose::Bool,io::IO)
    skeleton=getskeleton_seg(doseerrls,snporder,priorprocess,skeletonsize)
    skeletonprior = deepcopy(priorprocess)
    pri1=first(values(skeletonprior))
    if 0 in pri1.markerincl[skeleton]
        @error string("deleted markered cannot in skeleton")
    end
    incl = falses(length(pri1.markerincl))
    incl[skeleton] .= true
    for (strkey,pri) in skeletonprior
        setmarkerincl!(pri,incl)
    end
    s1=sum(first(values(priorprocess)).markerdeltd)
    s2=sum(first(values(skeletonprior)).markerdeltd)
    s1 ≈ s2 || @error string("inconsistent markerdeltd. s1=",s1,",s2=",s2)
    markerexcl=.!(first(values(skeletonprior)).markerincl)
    doseerrls[snporder[markerexcl]] .= missing
    fhaplo= getfhaplo(polygeno.parentgeno[chr])
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    chrdose = polygeno.offspringgeno[chr]
    for it=1:5
        startt = time()
        doseerrls,logleps = updatedoseerrls!(doseerrls, snporder,bvpair,
            chrdose,deriveddose,priorspace, skeletonprior,polygeno)
        logldis = updatedistance!(snporder,bvpair,priorspace,skeletonprior,
            chrdose, deriveddose,doseerrls, polygeno)
        pri1=first(values(skeletonprior))
        chrlen = round(100*sum(pri1.markerdeltd),digits=1)
        logl =round(logldis,digits=1)
        msg =string("chr=",chrid,
            # nprocs()>1 ? string(", p=",myid()) : "",
            ", it=",it,", logl=",logl,", len=",chrlen,
            ", skeleton=", length(skeleton),
            ", <err>=",round(mean(skipmissing(doseerrls)),digits=4),
            ", tused=",round(time()-startt,digits=1),"s",
            )
        printconsole(io,verbose,msg)
    end
    skeletonprior
end

function rescalemap!(priorprocess::AbstractDict,skeletonprior::AbstractDict)
    pri0=first(values(priorprocess))
    incl0 = findall(pri0.markerincl)
    deltd0=deepcopy(pri0.markerdeltd)
    pri1=first(values(skeletonprior))
    incl1 = findall(pri1.markerincl)
    deltd1=pri1.markerdeltd
    issubset(incl1,incl0) || @error string("unexpected skeleton incl")
    incl0[1] in incl1 || @error string("first marker must be in skeleton")
    incl0[end] in incl1 || @error string("last marker must be in skeleton")
    for i=1:length(incl1)-1
        s=incl1[i]:incl1[i+1]-1
        temp=sum(deltd0[s])
        scale = temp ≈ 0 ? 0 : sum(deltd1[s])/temp
        deltd0[s] .*= scale
    end
    # sum(deltd0) ≈ sum(deltd1) || @error string("wrong re-scale by skeleton map")
    for (strkey, pri) in priorprocess
        pri.markerdeltd = deltd0
        pri.tranprobseq = vcat([getgametetran(i,pri.nvalent) for i=deltd0[1:end-1]],missing)
    end
    priorprocess
end

function delOutlier!(polygeno::PolyGeno)
    outlier = polygeno.offspringinfo[!,:isoutlier]
    any(ismissing.(outlier)) && return polygeno
    if sum(outlier)>0
        incl = .! outlier
        polygeno.offspringgeno=[i[:,incl] for i=polygeno.offspringgeno]
        polygeno.offspringinfo=polygeno.offspringinfo[incl,:]
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
