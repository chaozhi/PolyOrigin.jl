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
    if io!=nothing
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
    epsilon::Real=0.01,seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=false,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,maxiter::Integer=30,
    minrun::Integer=3,maxrun::Integer=10,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    polygeno = readPolyGeno(genofile,pedfile,
        delimchar=delimchar,
        missingstring=missingstring,
        commentstring=commentstring,
        workdir=workdir,verbose=false)
    polyPhase!(polygeno, epsilon=epsilon,seqerr=seqerr,
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

`epsilon::Real=0.01`: genotypic error probability.

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

`isparallel::Bool=false`: if true, multicore computing over chromosomes.

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

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyPhase!(polygeno::PolyGeno;
    epsilon::Real=0.01,seqerr::Real=0.001,
    chrpairing_phase::Integer=22,
    byparent::Union{Nothing,Bool}=nothing,
    byneighbor::Union{Nothing,Bool}=nothing,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=false,
    delmarker::Bool=true,
    delsiglevel::Real=0.05,
    maxstuck::Integer=5,maxiter::Integer=30,
    minrun::Integer=3,maxrun::Integer=10,
    refhapfile::Union{Nothing,AbstractString} = nothing,
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
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
    printconsole(io,verbose,string("PolyOrigin, polyPhase, logfile=", logfile, ", ", Dates.now()))
    minrun > maxrun && (minrun = maxrun)
    msg = string("list of option values: \n",
        "epsilon = ", epsilon, "\n",
        "seqerr = ",seqerr,"\n",
        "chrpairing_phase = ",chrpairing_phase, "\n",
        "byparent = ", byparent == nothing ? "true for simple biparental cross" : byparent,"\n",
        "byneighbor = ", byneighbor == nothing ? "true for max(ploidy) >= 6" : byneighbor,"\n",
        "chrsubset = ", chrsubset == nothing ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", snpsubset == nothing ? "all markers" : snpsubset,"\n",
        "isparallel = ", isparallel, "\n",
        "delmarker = ", delmarker, "\n",
        "delsiglevel = ", delsiglevel, "\n",
        "maxstuck = ",maxstuck, "\n",
        "maxiter = ",maxiter,"\n",
        "minrun = ",minrun, "\n",
        "maxrun = ",maxrun,"\n",
        "refhapfile = ",refhapfile,"\n",
        "missingstring = ",missingstring,"\n",
        "outstem = ",outstem == nothing ? "no output files" : outstem,"\n",
        "logfile = ",io,"\n",
        "workdir = ",workdir,"\n",
        "verbose = ",verbose)
    printconsole(io,false,msg)
    msg = string("data: #pop=", size(polygeno.designinfo,1),
        ", #parent=",size(polygeno.parentinfo,1),
        ", #offspring=",size(polygeno.offspringinfo,1),
        ", #chr=",length(polygeno.markermap),
        ", #marker=",sum(size.(polygeno.markermap,1)))
    printconsole(io,verbose,msg)
    if chrsubset!=nothing
        verbose && @info string("chrsubset=",chrsubset)
    end
    if snpsubset!=nothing
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
        if byparent == nothing
            # byparent2 = size(subpolygeno.designinfo,1) > size(subpolygeno.parentinfo,1)
            byparent2 = !(size(polygeno.designinfo,1) ==1 && size(polygeno.parentinfo,1) ==2)
            printconsole(io,verbose,string("set byparent=",byparent2))
        else
            byparent2 = byparent
        end
        if byneighbor == nothing
            # byparent2 = size(subpolygeno.designinfo,1) > size(subpolygeno.parentinfo,1)
            byneighbor2 = max(subpolygeno.parentinfo[!,:ploidy]...) >=6
            # printconsole(io,verbose,string("set byneighbor=",byneighbor2))
        else
            byneighbor2 = byneighbor
        end
        phases =parentalphasing_allchr(subpolygeno,epsilon,chrpairing_phase,
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
    if refhapfile != nothing
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
        round(time()-starttime), " seconds by polyPhase"))
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
    epsilon::Real,chrpairing_phase::Integer,
    byparent::Bool,byneighbor::Bool,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,
    minrun::Integer,maxrun::Integer,
    isparallel::Bool,io::Union{Nothing,IO},verbose::Bool)
    nchr=length(polygeno.markermap)
    parenthaplo= Vector{Vector}(undef,nchr)
    if isparallel && nprocs()>1
        polygenols=[getsubPolyGeno(polygeno,chrsubset=[chr]) for chr=1:nchr]
        res = pmap(x->parentalphasing_chr(x,1,epsilon,chrpairing_phase,
            delmarker,delsiglevel,maxstuck,maxiter,
            minrun,maxrun,byparent,byneighbor,nothing,verbose),polygenols)
        for chr=1:nchr
            parenthaplo[chr], iobuffer = res[chr]
            write(io,String(take!(iobuffer)))
            flush(io)
            close(iobuffer)
        end
    else
        for chr=1:nchr
            parenthaplo[chr]=first(parentalphasing_chr(polygeno,chr,epsilon,chrpairing_phase,
                delmarker,delsiglevel,maxstuck,maxiter,
                minrun,maxrun,byparent,byneighbor,io,verbose))
        end
    end
    parents = polygeno.parentinfo[!,:individual]
    Dict([parents[i] => [j[i] for j=parenthaplo] for i=1:length(parents)])
end

function parentalphasing_chr(polygeno::PolyGeno, chr::Integer,
    epsilon::Real,chrpairing_phase::Integer,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,
    minrun::Integer,maxrun::Integer,
    byparent::Bool, byneighbor::Bool,
    io::Union{Nothing,IO},verbose::Bool)
    starttime = time()
    io===nothing && (io=IOBuffer(append=true))
    priorspace = getpriorstatespace(polygeno,chrpairing_phase)
    chrdose=polygeno.offspringgeno[chr]
    nmarker = size(polygeno.markermap[chr],1)
    chrid = polygeno.markermap[chr][1,:chromosome]
    priorprocess = getpriorprocess(polygeno,chr,chrpairing_phase)
    fhaploset, fhaploweight=calparenthaploset(polygeno,chr)
    phaselogl = Vector{Float64}(undef,0)
    phaseres = Vector{Vector}(undef,0)
    batchsize=minrun
    run=0
    while true
        for i = 1:batchsize
            run+=1
            # phaseres is a list of (fhaploindex, bvpair,logl, loglhis)
            res = parentalphasing_local(fhaploset,fhaploweight,epsilon,chrdose,
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
            msg = string("#chr=",chrid,
                nprocs()>1 ? string(", procs=",myid()) : "",
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
    [hcat(map((x,y)->ismissing(y) ? missings(length(x[1])) : x[y],
        fhaploset[i],haploindex[i])...)' for i=1:length(fhaploset)]
end

function parentalphasing_local(fhaploset::AbstractVector,fhaploweight::AbstractVector,
    epsilon::Real,chrdose::AbstractMatrix,
    priorspace::AbstractDict,inputpriorprocess::AbstractDict,polygeno::PolyGeno,
    delmarker::Bool,delsiglevel::Real,
    maxstuck::Integer,maxiter::Integer,byparent::Bool,byneighbor::Bool,
    chrid::AbstractString,run::Integer,io::IO,verbose::Bool)
    priorprocess = deepcopy(inputpriorprocess)
    # fhaploindex=[[rand(Categorical(i)) for i=j] for j=fhaploweight]
    fhaploindex=[Vector{Union{Missing,Int}}([rand(Categorical(i)) for i=j]) for j=fhaploweight]
    noff = size(polygeno.offspringinfo,1)
    bvpair = randinitbvpair(priorspace,polygeno)
    siblogl = zeros(noff)
    popidls = keys(priorspace)
    if byneighbor
        updatebvpair!(bvpair, siblogl,fhaploindex,fhaploset,epsilon,chrdose,priorspace,
            priorprocess,polygeno,popidls)
    else
        randbvpair!(bvpair, siblogl,fhaploindex,fhaploset,epsilon,chrdose,priorspace,
            priorprocess,polygeno,popidls)
    end
    logl = sum(siblogl)
    # msg = string("#chr=",chrid,", run=",run, ", it=0, logl=",round(logl,digits=2))
    # printconsole(io,verbose,msg)
    # https://stackoverflow.com/questions/24492429/julia-controlling-where-cursor-is-in-printed-output
    # verbose && print(string(msg, "\r")) # "\r" go the first column
    # http://julia-programming-language.2336112.n4.nabble.com/How-to-flush-output-of-print-td19401.html
    # print(io, "\u1b[1G")   # go to first column
    # print(io, "\u1b[K")    # clear the rest of the line
    # verbose && print("\u1b[K",msg,"\u1b[1G")
    loglhis =[logl]
    nstuck =0
    for it=1:maxiter
        newfhaploindex,newbvpair,newsiblogl,newlogl=updatefhaplobvpair(fhaploindex,
            bvpair,siblogl,fhaploset,fhaploweight,epsilon,chrdose,
            priorspace,priorprocess,polygeno,byparent,byneighbor)
        logl = loglhis[end]
        accept = newlogl >= logl
        if accept
            newlogl == logl ? nstuck += maxstuck : nstuck=0
            fhaploindex,bvpair,siblogl,logl = newfhaploindex,newbvpair,newsiblogl,newlogl
        else
            nstuck +=1
            cond=all(map((x,y)->all(skipmissing(x .== y)),newfhaploindex,fhaploindex))
            cond && (nstuck += 1)
        end
        push!(loglhis,logl)
        incl = first(values(priorprocess)).markerincl
        ndel= length(incl)-sum(incl)
        msg = string("#chr=",chrid,
            nprocs()>1 ? string(", procs=",myid()) : "",
            ", run=", run,", it=", it,
            ", logl=",round(logl,digits=2),
            ", epsilon=", round(epsilon,digits=4),
            ", nstuck=", nstuck,
            ", ndel=", ndel
            )
        if nstuck >= maxstuck || it == maxiter
            # verbose && print("\u1b[1K")
            printconsole(io,verbose,msg)
            break
        else
            if nstuck==1 && delmarker
                epsilon=first(updateepsilon(epsilon,bvpair,fhaploindex,fhaploset,chrdose,
                    priorspace,priorprocess,polygeno))
                dataprobset = caldataprobset(fhaploindex,fhaploset,epsilon,
                    chrdose,priorspace,polygeno)
                # priorprocess is modified
                polymarkerdel!(fhaploindex,fhaploset,dataprobset,bvpair,
                    priorspace, priorprocess,polygeno,delsiglevel=delsiglevel)
                delmarker = false
            end
            printconsole(io,verbose,msg)
            # verbose && print("\u1b[K",msg,"\u1b[1G")
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

function inferepsilonls!(epsilonls::AbstractVector,chrdose::AbstractMatrix,
    bvpair::AbstractVector,fhaploindex::AbstractVector,fhaploset::AbstractVector,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno)
    snporder = 1:size(chrdose,1)
    dataprobset = caldataprobset(fhaploindex,fhaploset,epsilonls,
        chrdose,priorspace,polygeno)
    fhaplo= getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno)
    updateepsilonls!(epsilonls, snporder,dataprobset,bvpair,
        chrdose,deriveddose,priorspace, priorprocess,polygeno)
end

function loglepsilon(epsilon::Real,bvpair::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    chrdose::AbstractMatrix, priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    dataprobset = caldataprobset(fhaploindex,fhaploset,epsilon,
        chrdose,priorspace,polygeno)
    bvpairprop = [[i] for i=bvpair]
    logllist = calmarglogl(dataprobset, priorspace,priorprocess,polygeno,bvpairprop)
    logl=sum(vcat(logllist...))
    logl, dataprobset
end

function updateepsilon(epsilon::Real,bvpair::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    accuracygoal, precisiongoal, itmax = 2, 2, 20
    lowbound,upbound = log(10^(-6.0)), log(1.0-10^(-6.0))
    # here x is Jaccobi factor for transformation eps->log(eps)
    f(x)=first(loglepsilon(exp(x),bvpair,fhaploindex,fhaploset,chrdose,
        priorspace,priorprocess,polygeno))
    x0=log(epsilon)
    newx, logleps, his=brentMax(f,lowbound,upbound,x0,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)
    # newx,logleps,accept =metroplos_norm(f,x0,maxiter=10,
    #     temperature=temperature,propstd=log(3))
    exp(newx),logleps
end

function updatefhaplobvpair(fhaploindex::AbstractVector,
    bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    epsilon::Real,chrdose::AbstractMatrix,priorspace::AbstractDict,
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
    for findex=findexlist
        # println("findex=",findex)
        newfphase = randfphase(findex,newfhaploindex,
            fhaploset,fhaploweight,epsilon,
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
                epsilon,chrdose,priorspace,priorprocess,polygeno,popidls)
        else
            randbvpair!(newbvpair,newsiblogl,newfhaploindex,fhaploset,
                epsilon,chrdose,priorspace,priorprocess,polygeno,popidls)
        end
    end
    newfhaploindex,newbvpair,newsiblogl,round(sum(newsiblogl),digits=2)
end

function randfphase(findex,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    epsilon::Real,bvpair::AbstractVector,chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno)
    isreverse = rand([true,false])
    if isreverse
        reverse!.(fhaploindex)
        reverse!.(fhaploset)
        reverse!.(fhaploweight)
        chrdose=reverse(chrdose, dims=1)
        # priorprocess = Dict([key => reverseprior(val) for (key, val) in priorprocess])
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    fwphaseset,fwlogpgeno,fwlogpost,bvdict=polyforward(findex,fhaploindex,
        fhaploset,fhaploweight,epsilon,bvpair,chrdose,priorspace,priorprocess,polygeno)
    fphases = first(polybackward(fwlogpgeno,fwlogpost,bvdict,priorprocess))
    fphases = map((x,y)->ismissing(y) ? missing : x[y],fwphaseset,fphases)
    if ndims(findex) == 0
        res = fphases
    elseif ndims(findex) == 1
        nfindex = length(findex)
        if nfindex==1
            res = [fphases]
        elseif nfindex==2
            hh = map((i,j)->[div(i-1,j)+1,rem(i-1,j)+1],fphases,length.(fhaploset[last(findex)]))
            hh = hcat(hh...)
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
        # priorprocess = Dict([key => reverseprior(val) for (key, val) in priorprocess])
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
    dataprob::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popid::AbstractString)
    itmax=10
    nstuck=0
    loglnow = -Inf
    for it=1:itmax
        dims = rem(it,3)
        if dims == 0
            bvprop = getbvpairprop(priorspace,polygeno,popid,loglgrid)
        else
            # bvprop = getbvpairprop(priorspace,polygeno,popid,dims,bvnow)
            bvprop = priorspace[popid]["valentneighbor"][bvnow][dims]
        end
        # loglls=loglgrid[bvprop]
        loglls = upoffmarglogl!(loglgrid,dataprob, popid, priorspace,priorprocess,bvprop)
        logl = max(loglls...)
        bv = rand(bvprop[loglls .== logl])
        # println("   it=",it,",(bv,logl)=",(bv,logl))
        (bv == bvnow && logl == loglnow) ? nstuck +=1 : nstuck = 0
        bvnow = bv
        loglnow = logl
        nstuck>=2 && break
    end
    bvnow
end

function randbvprop(priorspace::AbstractDict,polygeno::PolyGeno,
    popid::AbstractString,loglgrid::AbstractMatrix)
    valents = priorspace[popid]["valent"]
    isnonmiss = .!ismissing.(valents)
    cart = CartesianIndices(valents)
    linear = LinearIndices(valents)
    bool = ismissing.(loglgrid) .* isnonmiss
    res = linear[bool]
    length(res) == 0 ? missing : rand(res)
end

function randbvprop(bvnow::Integer,loglnow::Real,dataprob::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popid::AbstractString)
    loglgrid = Matrix{Union{AbstractFloat,Missing}}(missing,size(priorspace[popid]["valent"]))
    minrun = 3
    maxrun = 10
    loglhis = [-Inf]
    bvhis = [0]
    for run=1:maxrun
        bvnow = randbvprop(priorspace,polygeno,popid,loglgrid)
        ismissing(bvnow) && break
        bv = uploglgrid!(loglgrid,bvnow,dataprob,priorspace,priorprocess,
            polygeno,popid)
        logl = loglgrid[bv]
        push!(bvhis,bv)
        push!(loglhis,logl)
        nrepeat =sum(loglhis .== max(loglhis...))
        nrepeat == minrun && break
        # println("run=",run, ",bv=",bv, ",logl=", logl,",nrepeat=",nrepeat)
    end
    bv,logl=last(bvhis), last(loglhis)
    if logl < loglnow
        bv = uploglgrid!(loglgrid,bv,dataprob,priorspace,priorprocess,
            polygeno,popid)
        logl = loglgrid[bv]
    end
    bv, logl
end

function updatebvpair!(bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    epsilon::Union{Real,AbstractVector},chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popidls)
    # startt=time()
    fhaplo=getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno,popidls)
    for popid = popidls
        offls = findall(polygeno.offspringinfo[!,:population] .== popid)
        for off=offls
            # println("off=", off,",t=",round(time()-startt,digits=2))
            offdose = chrdose[:,off]
            popid = polygeno.offspringinfo[off,:population]
            ploidy = polygeno.offspringinfo[off,:ploidy]
            dataprob = caldataprob(offdose,popid,ploidy,deriveddose,epsilon)
            bvpair[off],siblogl[off] = randbvprop(bvpair[off],siblogl[off],dataprob,priorspace,priorprocess,polygeno,popid)
        end
    end
    markerincl = first(values(priorprocess)).markerincl
    tseq = findall(.!markerincl)
    if length(tseq)>0
        isologl = calisomarkerlogl(tseq,fhaploindex,fhaploset,chrdose,epsilon,
            bvpair,priorspace,priorprocess,polygeno)
        offls = findall([i in popidls for i=polygeno.offspringinfo[!,:population]])
        siblogl[offls] .+= isologl[offls]
    end
    bvpair,siblogl
end

function randbvpair!(bvpair::AbstractVector,siblogl::AbstractVector,
    fhaploindex::AbstractVector,fhaploset::AbstractVector,
    epsilon::Union{Real,AbstractVector},chrdose::AbstractMatrix,priorspace::AbstractDict,
    priorprocess::AbstractDict,polygeno::PolyGeno,popidls;isrand::Bool=true)
    fhaplo=getfhaplo(fhaploindex,fhaploset)
    deriveddose = getderiveddose(fhaplo,priorspace,polygeno,popidls)
    bvpairprop = getbvpairprop(priorspace,polygeno)
    logllist = calmarglogl(epsilon,deriveddose,chrdose,priorspace,
        priorprocess,polygeno,bvpairprop)
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
        isologl = calisomarkerlogl(tseq,fhaploindex,fhaploset,chrdose,epsilon,
            bvpair,priorspace,priorprocess,polygeno)
        offls = findall([i in popidls for i=polygeno.offspringinfo[!,:population]])
        siblogl[offls] .+= isologl[offls]
    end
    bvpair,siblogl
end

function calisomarkerlogl(tseq::AbstractVector,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,chrdose::AbstractMatrix,
    epsilon::Union{Real,AbstractVector},bvpair::AbstractVector,priorspace::AbstractDict,
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
        calisomarkerlogl(markerfphase,markeroffdose,epsilon,bvpair,priorspace,priorprocess,polygeno)
    end for t=tseq)
end

function calisomarkerlogl(markerfphase::AbstractVector,markeroffdose::AbstractVector,
    epsilon::Union{Real,AbstractVector},bvpair::AbstractVector,priorspace::AbstractDict,
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
                pri = priorprocess[keyls[bv]]
                dataprob0 = caldataprob([markeroffdose[off]],popid,ploidy,deriveddose,epsilon)
                dataprob = dataprob0[1,condstates[bv]]
                startprob = pri.startprob
                log(dot(startprob,dataprob))
            end for off=offls]
        end
    end
    maximum(res,dims=2)[:,1]
end

function polyforward(findex,fhaploindex::AbstractVector,
    fhaploset::AbstractVector,fhaploweight::AbstractVector,
    epsilon::Real,bvpair::AbstractVector,chrdose::AbstractMatrix,
    priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno)
    #
    # startt=time()
    # println("startt time = ", startt)
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
    dataprob = calsitedataprob(t,findex,epsilon,bvpair,deriveddose,chrdose,priorspace,polygeno)
    # orde4ring of keys(bvdict) = popfromfindex, the same as that for dataprob
    startprob=[[priorprocess[i].startprob for i=bvdict[id][2]] for id=keys(bvdict)]
    # size(dataprobls[k][ind],1)is the same for any ind=1,...
    fwlogpost[t]=[[log.(startprob[k][ind] .* dataprob[k][ind][g,:]) for
        g=1:size(dataprob[k][1],1),ind=1:length(dataprob[k])] for k=1:length(startprob)]
    fwlogpost[t], fwlogpgeno[t], fwphaseset[t]  = normfwpost(fwlogpost[t],phaseweight[t])
    for marker=2:length(tseq)
        t=tseq[marker]
        tbef=tseq[marker-1]
        dataprob = calsitedataprob(t,findex,epsilon,bvpair,deriveddose,chrdose,
                    priorspace,polygeno)
        tranprob=[[priorprocess[i].tranprobseq[tbef] for i=bvdict[id][2]] for id=keys(bvdict)]
        fwlogpost[t] =[begin
            ls4=[begin
                ls = log.(exp.(hcat(fwlogpost[tbef][k][:,ind]...)') * tranprob[k][ind])
                ls2 = [logsumexp(ls[:,i] + fwlogpgeno[tbef])  for i=1:size(ls,2)]
                ls3 = log.(dataprob[k][ind])
                [ls3[i,:] + ls2 for i=1:size(ls3,1)]
            end for ind=1:length(dataprob[k])]
            # ls4=[begin
            #     ls = [exp(fwlogpgeno[t-1][g]) .* exp.(fwlogpost[t-1][k][g,ind])
            #             for g=1:length(fwlogpgeno[t-1])]
            #     ls2 = sum(hcat(ls...)' * tranprob[k][ind],dims=1)
            #     ls3 = log.(dataprob[k][ind])
            #     ls3 .+= repeat(log.(ls2),size(ls3,1))
            #     [ls3[i,:] for i=1:size(ls3,1)]
            # end for ind=1:length(dataprob[k])]
            hcat(ls4...)
        end for k=1:length(tranprob)]
        fwlogpost[t], fwlogpgeno[t], fwphaseset[t]  = normfwpost(fwlogpost[t],phaseweight[t])
    end
    fwphaseset,fwlogpgeno,fwlogpost,bvdict
end

function getphaseweight(findex, fhaploweight::AbstractVector)
    if typeof(findex) <:AbstractVector
        if length(findex)>1
            phaseweight0 = hcat(fhaploweight[findex]...)
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
    phase[tend]=rand(Categorical(exp.(fwlogpgeno[tend])))
    orig[tend]=[begin
        logpost=(fwlogpost[tend][k])[phase[tend],:]
        [rand(Categorical(exp.(i))) for i=logpost]
    end for k=1:length(fwlogpost[tend])]
    for marker=length(tseq)-1:-1:1
        t=tseq[marker]
        tafter=tseq[marker+1]
        tranprob=[[priorprocess[i].tranprobseq[t] for i=bvdict[id][2]] for id=keys(bvdict)]
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
    logpost=hcat(sitelogpost...)
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
