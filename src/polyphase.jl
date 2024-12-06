function connected_polygenos(polygeno::PolyGeno)
    des = polygeno.designinfo
    ccparent = connected_parents(des)
    length(ccparent) ==1 && return [polygeno]
    ccpop = [findall(sum(Matrix(des[:,i .+ 1]),dims=2)[:,1] .> 0) for i= ccparent]
    ccpopid=[des[!,:population][i] for i=ccpop]
    offpop=polygeno.offspringinfo[!,:population]
    ccoff=[findall([i in j for i=offpop]) for j=ccpopid]
    [begin
        ppgeno = [i[:,ccparent[c]] for i=polygeno.parentgeno]
        offgeno = [i[:,ccoff[c]] for i=polygeno.offspringgeno]
        ppinfo = polygeno.parentinfo[ccparent[c],:]
        offinfo =  polygeno.offspringinfo[ccoff[c],:]
        cols = ccparent[c] .+ 1
        pushfirst!(cols,1)
        desinfo = polygeno.designinfo[ccpop[c],cols]
        PolyGeno(polygeno.markermap,ppgeno,offgeno,ppinfo,offinfo,desinfo,
            polygeno.delmarker,polygeno.correction)
    end for c=1:length(ccparent)]
end

function connected_parents(designinfo::DataFrame)
    np = size(designinfo,2)-1
    adj = zeros(Int, np,np)
    for i=1:size(designinfo,1)
        pp = findall(Vector(designinfo[i,2:end]) .> 0)
        if length(pp)==2
            adj[pp[1],pp[2]] =1
        end
    end
    adj=sign.(adj + adj')
    g = SimpleGraph(adj)
    cc = connected_components(g)
end


function printconsole(io::Union{Nothing,IO},verbose::Bool,msg::AbstractString)
    verbose && @info(msg)
    if !isnothing(io)
        write(io,string(msg,"\n"))
        flush(io)
    end
end


"""
    polyPhase(genofile,pedfile, delimchar=',', missingstring="NA",
        commentstring="#", keyargs...)

performs parental phasing and return phasedgeno::PolyGeno
with phasedgeno.parentgeno being phased.

# Positional arguments

`genofile::AbstractString`: filename for genotypic data file.

`pedfile::AbstractString`:  filename for population pedigree.

see [`readPolyGeno`](@ref) for the requirements of `genofile` and `pedfile`.

# Keyword arguments

`delimchar::AbstractChar=','`:  text delimiter.

`missingstring::AbstractString="NA"`: string code for missing value.

`commentstring::AbstractString="#"`: rows that begins with commentstring will be ignored.

see keyargs in polyPhase!(polygeno::PolyGeno, keyargs...)

"""
function polyPhase(genofile::AbstractString,pedfile::AbstractString;
    delimchar::AbstractChar=',', missingstring::AbstractString="NA",
    commentstring::AbstractString="#",
    doseerr::Real=0.01,seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,maxiter::Integer=30,
    minrun::Integer=3,maxrun::Integer=10,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    polygeno = readPolyGeno(genofile,pedfile,
        delimchar=delimchar,
        missingstring=missingstring,
        commentstring=commentstring,
        workdir=workdir,verbose=false)
    polyPhase!(polygeno, doseerr=doseerr,seqerr=seqerr,
        chrpairing_phase=chrpairing_phase,
        chrsubset=chrsubset,
        snpsubset=snpsubset,
        isparallel=isparallel,
        delmarker=delmarker,
        delsiglevel=delsiglevel,
        maxstuck=maxstuck,maxiter=maxiter,minrun=minrun,
        maxrun=maxrun,byparent=byparent,byneighbor=byneighbor,
        refhapfile=refhapfile,
        missingstring=missingstring,
        outstem=outstem,
        logfile=logfile,workdir=workdir,verbose=verbose)
end

"""
    polyPhase!(polygeno::PolyGeno, keyargs...)

performs parental phasing and return phasedgeno::PolyGeno
with  phasedgeno.parentgeno and polygeno.parentgeno being phased .

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`doseerr::Real=0.01`: genotypic error probability.

`seqerror::Real=0.001`: base sequencing error probability for GBS data.

`chrpairing_phase::Integer=22`: chromosome pairing in parental phasing, with 22 being only
bivalent formations and 44 being bi- and quadri-valent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosomes
to be considered, with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`delmarker::Bool=true`: if true, delete markers during parental phasing.

"delsiglevel::Real=0.05": significance level for deleting markers.

`maxstuck::Integer=5`: the max number of consecutive iterations that are rejected
in a phasing run.

`maxiter::Integer=30`: the max number of iterations in a phasing run."

`minrun::Integer=3`: if the min number of phasing runs that are at the same local maximimum or
have the same parental phases reaches minrun, phasing algorithm will stop before reaching the maxrun.

`maxrun::Integer=10`: the max number of phasing runs.

`byparent::Union{Nothing,Bool}=nothing`: if true, update parental phases
 parent by parent; if false, update parental phases one subpopulation by subpopulation.
 The nothing denotes that it is true if a connected component is a simple F1 cross,
and false otherwise.

`byneighbor::Union{Nothing,Bool}=nothing`: if ture, udpate the combination of bivalent
or multivalents in parents by their neighbors; if false, consider all the possible combinations.
The nothing denotes thtat it is true if max ploidy>=6, and false otherwise.

`refhapfile::Union{Nothing,AbstractString} = nothing`: reference haplotype file
for setting absolute parental phases. It has the same format as the input genofile,
except that parental genotypes are phased and offspring genotypes are ignored if they exist.

`missingstring::AbstractString="NA"`: string code for missing value.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyPhase!(polygeno::PolyGeno;
    doseerr::Real=0.01,seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,maxiter::Integer=30,
    minrun::Integer=3,maxrun::Integer=10,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    starttime = time()
    if logfile === nothing
        io=nothing
    else
        if typeof(logfile) <: AbstractString
            io=open(getabsfile(workdir,logfile), "w+")
            printpkgst(io,verbose,"PolyOrigin")
        else
            io=logfile
        end
    end
    printconsole(io,verbose,string("polyPhase!, logfile=", logfile, ", ", Dates.now()))
    minrun > maxrun && (minrun = maxrun)
    msg = string("list of option values: \n",
        "doseerr = ", doseerr, "\n",
        "seqerr = ",seqerr,"\n",
        "chrpairing_phase = ",chrpairing_phase, "\n",
        "byparent = ", byparent,"\n",
        "byneighbor = ", isnothing(byneighbor) ? "true for max(ploidy) >= 6" : byneighbor,"\n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",
        "isparallel = ", isparallel, "\n",
        "delmarker = ", delmarker, "\n",
        "delsiglevel = ", delsiglevel, "\n",
        "maxstuck = ",maxstuck, "\n",
        "maxiter = ",maxiter,"\n",
        "minrun = ",minrun, "\n",
        "maxrun = ",maxrun,"\n",
        "refhapfile = ",refhapfile,"\n",
        "missingstring = ",missingstring,"\n",
        "outstem = ",isnothing(outstem) ? "no output files" : outstem,"\n",
        "logfile = ",io,"\n",
        "workdir = ",workdir,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
    msg = string("data: #pop=", size(polygeno.designinfo,1),
        ", #parent=",size(polygeno.parentinfo,1),
        ", #offspring=",size(polygeno.offspringinfo,1),
        ", #chr=",length(polygeno.markermap),
        ", #marker=",sum(size.(polygeno.markermap,1)))
    printconsole(io,verbose,msg)
    if !isnothing(chrsubset)
        verbose && @info string("chrsubset=",chrsubset)
    end
    if !isnothing(snpsubset)
        verbose && @info string("snpsubset=",snpsubset)
    end
    getsubPolyGeno!(polygeno,chrsubset=chrsubset,snpsubset=snpsubset)
    rawDoseCall!(polygeno,seqerr=seqerr)
    # if kindofgeno(polygeno.parentgeno) == "probability" || kindofgeno(polygeno.offspringgeno) == "probability"
    #     outfile = string(outstem,"_rawdoseprob.csv")
    #     outputfile2 = getabsfile(workdir,outfile)
    #     savegenodata(outputfile2,polygeno)
    #     msg = string("genofile with raw dose probability: ", outfile)
    #     printconsole(io,verbose,msg)
    # end
    res=[]
    polygenols=connected_polygenos(polygeno)
    for subpolygeno = polygenols
        pp = subpolygeno.parentinfo[!,:individual]
        pp = [i<length(pp) ? string(pp[i],", ") : pp[i] for i=1:length(pp)]
        msg = string("phasing connected parents: ", string(pp...))
        printconsole(io,verbose,msg)
        if isnothing(byparent)
            # byparent2 = size(subpolygeno.designinfo,1) > size(subpolygeno.parentinfo,1)
            b1 = max(subpolygeno.parentinfo[!,:ploidy]...) <= 4
            b2 = size(subpolygeno.parentinfo,2) ==2
            b3 = size(subpolygeno.designinfo,1) ==1
            byparent2 = !( b1 && b2 && b3)
            printconsole(io,verbose,string("set byparent=",byparent2))
        else
            byparent2 = byparent
        end
        if isnothing(byneighbor)
            # byparent2 = size(subpolygeno.designinfo,1) > size(subpolygeno.parentinfo,1)
            byneighbor2 = max(subpolygeno.parentinfo[!,:ploidy]...) >=6
            # printconsole(io,verbose,string("set byneighbor=",byneighbor2))
        else
            byneighbor2 = byneighbor
        end
        phases =parentalphasing_allchr(subpolygeno,doseerr,chrpairing_phase,
            byparent2,byneighbor2,delmarker,delsiglevel,
            maxstuck,maxiter,minrun,maxrun,isparallel,io,verbose)
        push!(res,phases)
    end
    phasedict = merge(res...)
    hh = [phasedict[i] for i=polygeno.parentinfo[!,:individual]]
    nsnp = size.(first(hh),1)
    nchr = length(first(hh))
    hh2=[[h[chr][i,:]' for h=hh, i=1:nsnp[chr]]' for chr=1:nchr]
    polygeno.parentgeno = hh2
    isdel=[[all(ismissing.(unique(vcat(i...)))) for i=eachrow(j)] for j=hh2]
    delmarker = setparentphase!(polygeno,isdel)
    if size(delmarker,1)>=1
        msg = string("delete ", size(delmarker,1)," markers: ",
            join(sort(string.(delmarker[!,:marker])),", "))
        printconsole(io,verbose,msg)
    end
    if !isnothing(refhapfile)
        try
            setAbsPhase!(refhapfile,polygeno,
                workdir=workdir,io=io,verbose=verbose)
        catch err
            msg = "could not set absolute parental phases,"
            msg = string("warning: ", msg, " ", err)
            verbose && @warn msg
            printconsole(io,false,msg)
        end
    end
    if outstem !== nothing
        outfile = string(outstem,"_parentphased.csv")
        savegenodata(outfile,polygeno,
            missingstring=missingstring,workdir=workdir)
        msg = string("parentphased file: ", outfile)
        printconsole(io,verbose,msg)
    end
    printconsole(io,verbose,string("End, ", Dates.now(),", time elapsed = ",
        round(time()-starttime), " seconds by polyPhase!"))
    if typeof(logfile) <: AbstractString
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polygeno
end

function setparentphase!(polygeno::PolyGeno,isdel::AbstractVector)
    markermap = polygeno.markermap
    nchr = length(markermap)
    delmarker=vcat([markermap[ch][isdel[ch],[:marker,:chromosome,:position]] for ch=1:nchr]...)
    isnondel=[.!i for i=isdel]
    polygeno.offspringgeno = [polygeno.offspringgeno[i][isnondel[i],:]  for i=1:length(isnondel)]
    polygeno.parentgeno = [polygeno.parentgeno[i][isnondel[i],:]  for i=1:length(isnondel)]
    polygeno.markermap= [markermap[i][isnondel[i],:]  for i=1:length(isnondel)]
    polygeno.delmarker = vcat(polygeno.delmarker,delmarker)
    delmarker
end

function parentalphasing_allchr(polygeno::PolyGeno,
    doseerr::Real,chrpairing_phase::Integer,
    byparent::Bool,byneighbor::Bool,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,
    minrun::Integer,maxrun::Integer,
    isparallel::Bool,io::Union{Nothing,IO},verbose::Bool)
    nchr=length(polygeno.markermap)
    parenthaplo= Vector{Vector}(undef,nchr)
    if isparallel && nprocs()>1
        polygenols=[getsubPolyGeno(polygeno,chrsubset=[chr]) for chr=1:nchr]
        res = pmap(x->parentalphasing_chr(x,1,doseerr,chrpairing_phase,
            delmarker,delsiglevel,maxstuck,maxiter,
            minrun,maxrun,byparent,byneighbor,nothing,verbose),polygenols)
        for chr=1:nchr
            parenthaplo[chr], iobuffer = res[chr]
            write(io,String(take!(iobuffer)))
            flush(io)
            close(iobuffer)
        end
        @everywhere GC.gc()
    else
        for chr=1:nchr
            parenthaplo[chr]=first(parentalphasing_chr(polygeno,chr,doseerr,chrpairing_phase,
                delmarker,delsiglevel,maxstuck,maxiter,
                minrun,maxrun,byparent,byneighbor,io,verbose))
        end
    end
    parents = polygeno.parentinfo[!,:individual]
    Dict([parents[i] => [j[i] for j=parenthaplo] for i=1:length(parents)])
end

function parentalphasing_chr(polygeno::PolyGeno, chr::Integer,
    doseerr::Real,chrpairing_phase::Integer,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,
    minrun::Integer,maxrun::Integer,
    byparent::Bool, byneighbor::Bool,
    io::Union{Nothing,IO},verbose::Bool)
    starttime = time()
    isnothing(io) && (io=IOBuffer(append=true))
    priorspace = getpriorstatespace(polygeno,chrpairing_phase)
    chrdose=polygeno.offspringgeno[chr]
    nmarker = size(polygeno.markermap[chr],1)
    chrid = polygeno.markermap[chr][1,:chromosome]
    maxvalent = max(digits(chrpairing_phase)...)
    priorprocess = getpriorprocess(polygeno,chr,maxvalent)
    fhaploset, fhaploweight=calparenthaploset(polygeno,chr)
    phaselogl = Vector{Float64}(undef,0)
    phaseres = Vector{Vector}(undef,0)
    batchsize=minrun
    run=0
    while true
        for i = 1:batchsize
            run+=1
            # phaseres is a list of (fhaploindex, bvpair,logl, loglhis)
            res = parentalphasing_local(fhaploset,fhaploweight,doseerr,chrdose,
                priorspace,priorprocess,polygeno,delmarker,delsiglevel,
                maxstuck,maxiter,byparent,byneighbor,chrid,run,io,verbose)
            push!(phaseres,res)
            push!(phaselogl,round(res[3],digits=2))
        end
        # cunt by logl # ncount = count(x->x==max(phaselogl...),phaselogl)
        # count by phase
        best = argmax(phaselogl)
        haplist = [index2haplo(fhaploset,j[1]) for j=phaseres]
        dis = [sum(last(calpermhomolog(haplist[best],est))) for est=haplist]
        ncount = sum(dis .== 0)
        if (ncount >= minrun) || run>=maxrun
            msg = string("chr=",chrid,
                # nprocs()>1 ? string(", p=",myid()) : "",
                ", #marker=",nmarker,
                ", done, elapsed=", round(time()-starttime), " seconds",
                ", #mismatch(run",best," ~ runs) = ", dis)
            printconsole(io,verbose,msg)
            return haplist[best],io
        end
        batchsize= min(minrun - ncount,maxrun-run)
    end
end

function index2haplo(fhaploset::AbstractVector,haploindex::AbstractVector)
    [reduce(hcat,map((x,y)->ismissing(y) ? missings(length(x[1])) : x[y],
        fhaploset[i],haploindex[i]))' for i=1:length(fhaploset)]
end

function parentalphasing_local(fhaploset::AbstractVector,fhaploweight::AbstractVector,
    doseerr::Real,chrdose::AbstractMatrix,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,byparent::Bool,byneighbor::Bool,
    chrid::AbstractString,run::Integer,io::IO,verbose::Bool)
    # priorprocess = deepcopy(inputpriorprocess)
    # fhaploindex=[[rand(Categorical(i)) for i=j] for j=fhaploweight]
    fhaploindex=[Vector{Union{Missing,Int}}([rand(Categorical(i)) for i=j]) for j=fhaploweight]
    noff = size(polygeno.offspringinfo,1)
    bvpair = randinitbvpair(priorspace,polygeno)
    siblogl = repeat([-Inf],noff)
    popidls = keys(priorspace)
    doseerrls = [doseerr for i=1:size(chrdose,1)]
    if byneighbor
        updatebvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,priorspace,
            priorprocess,polygeno,popidls)
    else
        randbvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,priorspace,
            priorprocess,polygeno,popidls)
    end
    logl = sum(siblogl)
    loglhis =[logl]
    nstuck = 0
    for it=1:maxiter        
        startt = time()
        newfhaploindex,newbvpair,newsiblogl,newlogl=updatefhaplobvpair(fhaploindex,
            bvpair,siblogl,fhaploset,fhaploweight,doseerr,chrdose,
            priorspace,priorprocess,polygeno,byparent,byneighbor)
        logl = loglhis[end]
        accept = newlogl >= logl
        if accept
            newlogl == logl ? nstuck += (maxstuck-1) : nstuck=0
            fhaploindex,bvpair,siblogl,logl = newfhaploindex,newbvpair,newsiblogl,newlogl
        else
            nstuck +=1
            cond=all(map((x,y)->all(skipmissing(x .== y)),newfhaploindex,fhaploindex))
            cond && (nstuck += 1)
        end                        
        msg = string("chr=",chrid,", run=", run,", it=", it)                 
        if nstuck>0 && delmarker
            doseerr=first(updatedoseerr(doseerr,bvpair,fhaploindex,fhaploset,
                priorspace,priorprocess,chrdose,polygeno))
            doseerrls = [doseerr for i=1:size(chrdose,1)]
            # priorprocess is modified
            delsnps = polymarkerdel!(fhaploindex,fhaploset,bvpair,priorspace,
                priorprocess, chrdose, doseerrls, polygeno; delsiglevel)  
            ndel = length(delsnps)              
            if byneighbor
                updatebvpair!(bvpair,siblogl,fhaploindex,fhaploset,
                    doseerrls,chrdose,priorspace,priorprocess,polygeno,popidls)
            else
                randbvpair!(bvpair,siblogl,fhaploindex,fhaploset,
                    doseerrls,chrdose,priorspace,priorprocess,polygeno,popidls)
            end
            logl = round(sum(siblogl),digits=2)
            msg *= string(", logl=", round(logl,digits=1),                            
                ", stuck=", nstuck,                                 
                ", err=", round(doseerr,digits=4),                            
                ", tused=",round(time()-startt,digits=1),"s",
                ", #del=", ndel)
            ndel == 0 && (delmarker = false)               
        else
            msg *= string(", logl=",round(logl,digits=1),
                ", stuck=", nstuck,            
                ", err=", round(doseerr,digits=4),
                ", tused=",round(time()-startt,digits=1),"s") 
        end
        push!(loglhis,logl)                                                                                        
        printconsole(io,verbose,msg)            
        if nstuck >= maxstuck || it == maxiter
            break
        end
    end
    # keeping fhaploindex vlaue during local phasing
    # set fhaploindex to missing at deleted markers
    pri1=first(values(priorprocess))
    excl=.!(pri1.markerincl)
    for i=fhaploindex
        i[excl] .= missing
    end
    [fhaploindex, bvpair,loglhis[end], loglhis]
end

function logldoseerr(doseerr::Real,bvpair::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, deriveddose::AbstractDict,polygeno::PolyGeno)
    doseerrls = [doseerr for i=1:size(chrdose,1)]
    bvpairprop = [[i] for i=bvpair]
    logllist = calmarglogl(doseerrls, deriveddose, chrdose,
        priorspace,priorprocess,polygeno,bvpairprop)
    logl=sum(vcat(logllist...))
    logl
end

function updatedoseerr(doseerr::Real,bvpair::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    chrdose::AbstractMatrix, polygeno::PolyGeno)
    accuracygoal, precisiongoal, itmax = 3, 3, 20
    lowbound,upbound = log(10^(-6.0)), log(min(0.95,doseerr+0.25))
    # here x is Jaccobi factor for transformation eps->log(eps)
    fhaplo= getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    f(x)=logldoseerr(exp(x),bvpair,priorspace,priorprocess,
        chrdose,deriveddose,polygeno)
    x0=log(doseerr)
    newx, logleps, his=brentMax(f,lowbound,upbound,x0,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)
    # newx,logleps,accept =metroplos_norm(f,x0,maxiter=10,
    #     temperature=temperature,propstd=log(3))
    # @info "infer doseerror his: " his
    exp(newx),logleps
end

function updatefhaplobvpair(fhaploindex::AbstractVector,
    bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    doseerr::Real,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,byparent::Bool,byneighbor::Bool)
    if byparent
        nf= size(polygeno.designinfo,2)-1
        findexlist = 1:nf
    else
        npop = size(polygeno.designinfo,1)
        findexlist = [findall(Vector(polygeno.designinfo[popindex,2:end]) .> 0) for popindex=1:npop]
    end
    newfhaploindex = deepcopy(fhaploindex)
    newbvpair,newsiblogl = deepcopy(bvpair),deepcopy(siblogl)
    doseerrls = [doseerr for i=1:size(chrdose,1)]
    # tuse1=tuse2=0
    for findex=findexlist
        # println("findex=",findex)
        newfphase = randfphase(findex,newfhaploindex,
            fhaploset,fhaploweight,doseerr,
            newbvpair,chrdose,priorspace,priorprocess,polygeno)
        if byparent
            b = .!(ismissing.(newfphase))
            newfhaploindex[findex][b] .= newfphase[b]
        else
            for i=1:length(findex)
                b = .!(ismissing.(newfphase[i]))
                newfhaploindex[findex[i]][b] .= newfphase[i][b]
            end
        end
        # popidls = keys(priorspace)
        # println("findex=",findex,", torandbvpair")
        popidls = popfromfindex(findex,polygeno)
        if byneighbor
            updatebvpair!(newbvpair,newsiblogl,newfhaploindex,fhaploset,
                doseerrls,chrdose,priorspace,priorprocess,polygeno,popidls)
        else
            randbvpair!(newbvpair,newsiblogl,newfhaploindex,fhaploset,
                doseerrls,chrdose,priorspace,priorprocess,polygeno,popidls)
        end
        # println("\t[tfphase,tbvpair]=",[tuse1,tuse2])
    end
    newfhaploindex,newbvpair,newsiblogl,round(sum(newsiblogl),digits=2)
end

function randfphase(findex,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    doseerr::Real,bvpair::AbstractVector,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    isreverse = rand([true,false])
    if isreverse
        reverse!.(fhaploindex)
        reverse!.(fhaploset)
        reverse!.(fhaploweight)
        chrdose=reverse(chrdose, dims=1)
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    fwphaseset,fwlogpgeno,fwlogpost,bvdict=polyforward(findex,fhaploindex,
        fhaploset,fhaploweight,doseerr,bvpair,chrdose,priorspace,priorprocess,polygeno)
    fphases = first(polybackward(fwlogpgeno,fwlogpost,bvdict,priorprocess))
    fphases = map((x,y)->ismissing(y) ? missing : x[y],fwphaseset,fphases)
    if ndims(findex) == 0
        res = fphases
    elseif ndims(findex) == 1
        nfindex = length(findex)
        if nfindex==1
            res = [fphases]
        elseif nfindex==2
            hh0 = map((i,j)->[div(i-1,j)+1,rem(i-1,j)+1],fphases,length.(fhaploset[last(findex)]))
            hh = reduce(hcat,hh0)
            res=[hh[i,:] for i=1:size(hh,1)]
        else
            @error(string("too many #parents = ", nfindex, " in the ", popindex, "-th subpopulation"))
        end
    else
        @error(string("wrong dimension of findex = ", nfindex))
    end
    if isreverse
        ndims(findex) == 0 ? reverse!(res) : reverse!.(res)
        reverse!.(fhaploindex)
        reverse!.(fhaploset)
        reverse!.(fhaploweight)
        chrdose=reverse(chrdose, dims=1)
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    res
end

function randinitbvpair(priorspace::AbstractDict,polygeno::PolyGeno)
    [begin
        popid = polygeno.offspringinfo[off,:population]
        valents = priorspace[popid]["valent"]
        isnonmiss = .!ismissing.(valents)
        rand(LinearIndices(valents)[isnonmiss])
    end for off =1:size(polygeno.offspringinfo,1)]
end

function getbvpairprop(priorspace,polygeno)
    [begin
        popid = polygeno.offspringinfo[off,:population]
        valents = priorspace[popid]["valent"]
        isnonmiss = .!ismissing.(valents)
        LinearIndices(valents)[isnonmiss]
    end for off=1:size(polygeno.offspringinfo,1)]
end

function getbvpairprop(priorspace::AbstractDict,polygeno::PolyGeno,
    popid::AbstractString,dims::Integer, bvnow::Integer)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)
    cart = CartesianIndices(valents)
    linear = LinearIndices(valents) .* isnonmiss
    pos = cart[bvnow]
    if dims==1
        res = linear[pos[1],:]
        res = res[res .> 0]
    elseif dims==2
        res = linear[:,pos[2]]
        res = res[res .> 0]
    else
        @error string("wrong keyarg  dims =",dims)
    end
    res
end

function getbvpairprop(priorspace::AbstractDict,polygeno::PolyGeno,
    popid::AbstractString,loglgrid::AbstractMatrix)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)
    cart = CartesianIndices(valents)
    linear = LinearIndices(valents) .* isnonmiss
    # ntop = round(Int,length(loglgrid)^0.25)
    nb=priorspace[popid]["valentneighbor"]
    ntop=round(Int,mean(length.(nb[1])))
    isloglnonmiss = .!ismissing.(loglgrid)
    p = sortperm(loglgrid[isloglnonmiss],rev=true)
    (length(p) >= ntop) && (p = p[1:ntop])
    top = cart[isloglnonmiss][p]
    xtop = unique([i[1] for i=top])
    ytop = unique([i[2] for i=top])
    res=reshape(linear[xtop,ytop],:)
    res[res .> 0]
end

function uploglgrid!(loglgrid::AbstractMatrix,bvnow::Integer,
    dataprob::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popid::AbstractString,
    snporder::AbstractVector)
    # itmax=ceil(Int, sqrt(size(loglgrid,1)))
    itmax = size(loglgrid,1)
    nstuckmax=ceil(Int, itmax/5)
    nstuck=0
    loglnow = -Inf
    for it=1:itmax
        dims = rem(it,3)
        if dims == 0
            bvprop = getbvpairprop(priorspace,polygeno,popid,loglgrid)
        else
            bvprop = priorspace[popid]["valentneighbor"][bvnow][dims]
        end
        loglls = upoffmarglogl!(loglgrid,dataprob, popid, priorspace,priorprocess,bvprop,snporder)
        logl = max(loglls...)
        bv = rand(bvprop[loglls .== logl])
        # println("   it=",it,",(bv,logl)=",(bv,logl))
        (bv == bvnow && logl ≈ loglnow) ? nstuck +=1 : nstuck = 0
        bvnow = bv
        loglnow = logl
        nstuck>=nstuckmax && break
    end
    bvnow
end

function randbvprop(priorspace::AbstractDict,polygeno::PolyGeno,
    popid::AbstractString,loglgrid::AbstractMatrix)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)
    # cart = CartesianIndices(valents)
    linear = LinearIndices(valents)
    bool = ismissing.(loglgrid) .* isnonmiss
    res = linear[bool]
    length(res) == 0 ? missing : rand(res)
end

function randbvprop(bvinput::Integer,loglinput::Real,dataprob::AbstractVector,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popid::AbstractString,
    snporder::Union{Nothing,AbstractVector}=nothing)
    loglgrid = Matrix{Union{AbstractFloat,Missing}}(missing,size(priorspace[popid]["valent"]))
    minrun = 3
    maxrun = 10
    loglhis = [-Inf]
    bvhis = [0]
    isnothing(snporder) && (snporder=1:length(dataprob))
    for run=1:maxrun
        bvnow = randbvprop(priorspace,polygeno,popid,loglgrid)
        ismissing(bvnow) && break
        bv = uploglgrid!(loglgrid,bvnow,dataprob,priorspace,priorprocess,
            polygeno,popid,snporder)
        logl = loglgrid[bv]
        push!(bvhis,bv)
        push!(loglhis,logl)
        nrepeat =sum(loglhis .== max(loglhis...))        
        nrepeat == minrun && break
    end
    bv,logl=last(bvhis), last(loglhis)
    if logl < loglinput
        bv = uploglgrid!(loglgrid,bvinput,dataprob,priorspace,priorprocess,
            polygeno,popid,snporder)
        logl = loglgrid[bv]
    end
    bv, logl
end

function updatebvpair!(bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    doseerrls::AbstractVector,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popidls;
    snporder::Union{Nothing,AbstractVector}=nothing,
    show_progress::Bool=false)
    # startt=time()
    fhaplo=getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno,popidls)
    noff = sum([in(i, popidls) for i in polygeno.offspringinfo[!,:population]])
    progress = Progress(noff;
        enabled = show_progress,
        desc=string("Update chrpairing..."))
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        ploidy = polygeno.offspringinfo[first(offls),:ploidy]
        for off=offls
            # println("off=", off,",t=",round(time()-startt,digits=2))
            offdose = chrdose[:,off]
            popid = polygeno.offspringinfo[off,:population]
            dataprob = caldataprob(offdose,ploidy,deriveddose[popid],doseerrls)
            bvpair[off],siblogl[off] = randbvprop(bvpair[off],siblogl[off],
                dataprob,priorspace,priorprocess,polygeno,popid,snporder)
            next!(progress)
        end
    end
    markerincl = first(values(priorprocess)).markerincl
    tseq = findall(.!markerincl)
    if length(tseq)>0
        isologl = calisomarkerlogl(tseq,fhaploindex,fhaploset,chrdose,doseerrls,
            bvpair,priorspace,priorprocess,polygeno)
        offls = findall([i in popidls for i=polygeno.offspringinfo[!,:population]])
        siblogl[offls] .+= isologl[offls]
    end
    bvpair,siblogl
end

function randbvpair!(bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    doseerrls::AbstractVector,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popidls;
    isrand::Bool=true,
    snporder::Union{Nothing,AbstractVector}=nothing)
    fhaplo=getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno,popidls)
    bvpairprop = getbvpairprop(priorspace,polygeno)
    logllist = calmarglogl(doseerrls,deriveddose,chrdose,priorspace,
        priorprocess,polygeno,bvpairprop;snporder)
    offls = findall(.!(ismissing.(logllist)))
    # problist, siblogl[offls]=normalizelogl(logllist[offls])
    problist = first(normalizelogl(logllist[offls]))
    if isrand
        bvpair0=[rand(Categorical(i)) for i = problist]
    else
        bvpair0=[rand(findall(i .== max(i...))) for i = problist]
    end
    bvpair[offls]=map((x,y)->x[y],bvpairprop[offls],bvpair0)
    siblogl[offls]=[max(i...) for i = logllist[offls]]
    # bvpair0=[rand(findall(i .== max(i...))) for i = logllist[offls]]
    # siblogl[offls]=map((x,y)->x[y], logllist[offls],bvpair0)
    # bvpair[offls]=map((x,y)->x[y],bvpairprop[offls],bvpair0)
    markerincl = first(values(priorprocess)).markerincl
    tseq = findall(.!markerincl)
    if length(tseq)>0
        isologl = calisomarkerlogl(tseq,fhaploindex,fhaploset,chrdose,doseerrls,
            bvpair,priorspace,priorprocess,polygeno)
        offls = findall([i in popidls for i=polygeno.offspringinfo[!,:population]])
        siblogl[offls] .+= isologl[offls]
    end
    bvpair,siblogl
end

function calisomarkerlogl(tseq::AbstractVector,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,chrdose::AbstractMatrix,
    doseerrls::AbstractVector,bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    sum(begin
        # println(["t=",t, "; ",fhaploindex[i][t] for i=1:length(fhaploindex)])
        # markerfphase= [i[t] for i=fhaploset]
        # markerfphase=[begin
        #     j=fhaploindex[i][t]
        #     j2= ismissing(j) ? (1:length(fhaploset[i][t])) : [j]
        #     fhaploset[i][t][j2]
        # end for i=1:length(fhaploindex)]
        markerfphase=[fhaploset[i][t][[fhaploindex[i][t]]] for i=1:length(fhaploindex)]
        markeroffdose = chrdose[t,:]
        calisomarkerlogl(markerfphase,markeroffdose,doseerrls,bvpair,priorspace,priorprocess,polygeno)
    end for t=tseq)
end

function calisomarkerlogl(markerfphase::AbstractVector,markeroffdose::AbstractVector,
    doseerrls::AbstractVector,bvpair::AbstractVector,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    if length(markerfphase)==1
        fhaplolist = markerfphase[1]
    else
        fhaplolist = getouterlist(markerfphase...)
    end
    res= zeros(length(markeroffdose),length(fhaplolist))
    popidls = polygeno.designinfo[!,:population]
    for kk=1:length(fhaplolist)
        deriveddose = getderiveddose(fhaplolist[kk]',priorspace,polygeno,popidls)
        for popid=popidls
            offls = findall(polygeno.offspringinfo[!,:population] .== popid)
            ploidy = polygeno.offspringinfo[first(offls),:ploidy]
            condstates = priorspace[popid]["condstate"]
            keyls = priorspace[popid]["valentkey"]
            res[offls,kk]=[begin
                bv = bvpair[off]
                ddose = deriveddose[popid][:, condstates[bv]]
                dataprob = only(caldataprob([markeroffdose[off]],ploidy,ddose,doseerrls))
                startprob = only(getstartprobls(priorprocess,[keyls[bv]]))
                log(dot(startprob,dataprob))
            end for off=offls]
        end
    end
    maximum(res,dims=2)[:,1]
end

function polyforward(findex,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    doseerr::Real,bvpair::AbstractVector,chrdose::AbstractMatrix,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno)
    # startt=time()
    deriveddose = getderiveddose(findex,fhaploindex,fhaploset,priorspace,polygeno,bvpair)
    bvdict = getbvdict(findex,bvpair,priorspace,polygeno)
    phaseweight=getphaseweight(findex,fhaploweight)
    nseq=size(chrdose,1)
    fwlogpost = Vector{Union{Missing,Vector}}(missing,nseq)
    fwlogpgeno = Vector{Union{Missing,Vector}}(missing,nseq)
    fwphaseset = Vector{Union{Missing,Vector}}(missing,nseq)
    # markerincl is the same for all keys of priorprocess
    tseq=findall(first(values(priorprocess)).markerincl)
    t=tseq[1]
    dataprob = calsitedataprob(t,findex,doseerr,bvpair,deriveddose,chrdose,priorspace,polygeno)
    # orde4ring of keys(bvdict) = popfromfindex, the same as that for dataprob
    startprob=[getstartprobls(priorprocess, bvdict[id][2]) for id=keys(bvdict)]
    # size(dataprobls[k][ind],1)is the same for any ind=1,...
    fwlogpost[t]=[[Vector{Float64}(log.(startprob[k][ind] .* dataprob[k][ind][g,:])) for
        g=1:size(dataprob[k][1],1),ind=1:length(dataprob[k])] for k=1:length(startprob)]
    fwlogpost[t], fwlogpgeno[t], fwphaseset[t]  = normfwpost(fwlogpost[t],phaseweight[t])
    for marker=2:length(tseq)
        t=tseq[marker]
        tbef=tseq[marker-1]
        dataprob = calsitedataprob(t,findex,doseerr,bvpair,deriveddose,chrdose,
                    priorspace,polygeno)
        tranprob=[gettranprobls(priorprocess, bvdict[id][2],tbef) for id=keys(bvdict)]
        fwlogpost[t] =[begin
            ls4=[begin
                ls = log.(exp.(reduce(hcat,fwlogpost[tbef][k][:,ind])') * tranprob[k][ind])
                # logsumexp([-Inf]) = -Inf
                ls2 = [logsumexp(ls[:,i] + fwlogpgeno[tbef])  for i=1:size(ls,2)]
                ls3 = log.(dataprob[k][ind])
                ls4 = [ls3[i,:] + ls2 for i=1:size(ls3,1)]
                ls4
            end for ind=1:length(dataprob[k])]
            reduce(hcat, ls4)
        end for k=1:length(tranprob)]
        fwlogpost[t], fwlogpgeno[t], fwphaseset[t]  = normfwpost(fwlogpost[t],phaseweight[t])
    end
    fwphaseset,fwlogpgeno,fwlogpost,bvdict
end

function getphaseweight(findex, fhaploweight::AbstractVector)
    if typeof(findex) <:AbstractVector
        if length(findex)>1
            phaseweight0 = reduce(hcat,fhaploweight[findex])
            phaseweight = [kron(i...) for i=eachrow(phaseweight0)]
        elseif length(findex)==1
            phaseweight = fhaploweight[findex[1]]
        else
            @error string("wrong findex = ",findex)
        end
    elseif typeof(findex) <:Integer
        phaseweight = fhaploweight[findex]
    else
        @error string("wrong findex = ",findex)
    end
    phaseweight
end

function polybackward(fwlogpgeno::AbstractVector,fwlogpost::AbstractVector,
    bvdict::AbstractDict,priorprocess::AbstractDict)
    nseq = length(fwlogpgeno)
    phase = Vector{Union{Missing,Int}}(missing,nseq)
    orig = Vector{Union{Missing,Vector}}(missing,nseq)
    tseq = findall(.!ismissing.(fwlogpgeno))
    tend = tseq[end]
    pphase = exp.(fwlogpgeno[tend] .- logsumexp(fwlogpgeno[tend]))
    phase[tend]=rand(Categorical(pphase))
    orig[tend]=[begin
        logpost=(fwlogpost[tend][k])[phase[tend],:]
        [rand(Categorical(exp.(i .- logsumexp(i)))) for i=logpost]
    end for k=1:length(fwlogpost[tend])]
    for marker=length(tseq)-1:-1:1
        t=tseq[marker]
        tafter=tseq[marker+1]
        tranprob=[gettranprobls(priorprocess, bvdict[id][2],t) for id=keys(bvdict)]
        # #phases (=ng) depends on marker t, but not the subpopulation indnex k,
        ng = size(fwlogpost[t][1],1)
        weight = [ones(0) for i=1:ng]
        post=[[[begin
            ls = log.(tranprob[k][ind][:,orig[tafter][k][ind]])
            ls += fwlogpost[t][k][g,ind]
            w=logsumexp(ls)
            push!(weight[g],w)
            exp.(ls .- w)
        end for ind=1:length(tranprob[k])] for k=1:length(tranprob)] for g=1:ng]
        pphase = sum.(weight) + fwlogpgeno[t]
        pphase = exp.(pphase .- logsumexp(pphase))
        # phase[t]=rand(findall(pphase .== max(pphase...)))
        phase[t]=rand(Categorical(pphase))
        orig[t]=[[rand(Categorical(i)) for i=j] for j=post[phase[t]]]
    end
    phase,orig
end

function normalizelogl(logllist)
    # logllist is a list of list of logl for all possible pairing
    # (bivalent and multi-vavelent)for all offspring
    # equally probable for each of 9 or 16 the chromosome pairing
    # prior(V_0|H)=1/9 for bvModel, or 1/16 for fullModel for tetraploid
    logllist = [i .- log(length(i)) for i=logllist]
    siblogl=[logsumexp(i) for i=logllist]
    problist = [exp.(logllist[i] .- siblogl[i]) for i=1:length(logllist)]
    problist, siblogl
end

function normfwpost(sitelogpost::AbstractVector,weight::AbstractVector)
    logpost=reduce(hcat,sitelogpost)
    sibscale=logsumexp.(logpost)
    # cal logpgeno
    logpgeno=sum(sibscale,dims=2)[:,1]
    logpgeno .+= log.(weight)
    logpgeno .-= logsumexp(logpgeno)
    phaseset = findall(logpgeno .- max(logpgeno...) .> log(10^(-10.)))
    logpgeno = logpgeno[phaseset]
    # up sitelogpost
    map!((x,y)->x .- y ,logpost,logpost,sibscale)
    ls=size.(sitelogpost,2)
    pushfirst!(ls,0)
    cumsum!(ls,ls)
    sitelogpost = [logpost[phaseset,ls[i]+1:ls[i+1]] for i=1:length(ls)-1]
    sitelogpost, logpgeno, phaseset
end



function getbvdict(findex,bvpair::AbstractVector,
    priorspace::AbstractDict,polygeno::PolyGeno)
    popidls=popfromfindex(findex,polygeno)
    # using ordereddict so that keys(bvdict) = popidls
    # the ordering will be the same as startprob and transtionprob
    OrderedDict([begin
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        bvkeys = priorspace[popid]["valentkey"][bvpair[offls]]
        popid => (offls,bvkeys)
    end for popid = popidls])
end
