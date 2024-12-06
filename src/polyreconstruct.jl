"""
    polyReconstruct!(polygeno::PolyGeno, keyargs...)

performs ancestral inference from polygeno and return polyancestry::PolyAncestry.

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`epsilon::Real=0.01`: genotypic error probability for offspring phasing.

`seqerror::Real=0.001`: base sequencing error probability for GBS data.

`chrpairing::Integer=44`: chromosome pairing in offspring decoding, with 22 being only
bivalent formations and 44 being bivalent and quadrivalent formations.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosomes
to be considered, with nothing denoting all chromosomes.
Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of markers
to be considered, with nothing denoting all markers. within a chromosome, marker
index starts from 1, and marker indices that are larger than the number of markers
within the chromosome are deleted.

`isparallel::Bool=false`: if true, multicore computing over chromosomes.

`correctthreshold::AbstractFloat=0.15`: a candidate marker is selected for
parental error correction if the fraction of offspring genotypic error >= correctthreshold.

`isplot::Bool=false`: if true, plot haploprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyReconstruct!(phasedgeno::PolyGeno;
    epsilon::Real=0.01,seqerr::Real=0.001,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=false,
    correctthreshold::AbstractFloat=0.15,
    missingstring::AbstractString="NA",
    isplot::Bool=false,
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    verbose::Bool=true)
    starttime = time()
    if kindofgeno(phasedgeno.parentgeno) != "phasedgeno"
        @error "parental genotypes are not not phased "
    end
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
    getsubPolyGeno!(phasedgeno,chrsubset=chrsubset,snpsubset=snpsubset)
    rawDoseCall!(phasedgeno,seqerr=seqerr)
    printconsole(io,verbose,string("PolyOrigin, polyReconstruct!, ", Dates.now()))
    if chrsubset!=nothing
        verbose && @info string("chrsubset = ",chrsubset)
    end
    if snpsubset!=nothing
        verbose && @info string("snpsubset = ",snpsubset)
    end
    msg = string("list of option values: \n",
        "epsilon = ", epsilon, "\n",
        "seqerr = ",seqerr,"\n",
        "chrpairing = ",chrpairing, "\n",
        "chrsubset = ", chrsubset == nothing ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", snpsubset == nothing ? "all markers" : snpsubset,"\n",
        "isparallel = ", isparallel, "\n",
        "correctthreshold = ", correctthreshold, "\n",
        "missingstring = ", missingstring, "\n",
        "outstem = ", outstem == nothing ? "not save to files" : outstem,"\n",
        "logfile = ",logfile,"\n",
        "workdir = ", workdir,"\n",
        "verbose = ",verbose)
    printconsole(io,false,msg)
    priorspace = getpriorstatespace(phasedgeno,chrpairing)
    statespace = getstatespace(priorspace)
    nchr = length(phasedgeno.markermap)
    genoprob = Vector{Vector}(undef,nchr)
    valentprob = Vector{Vector}(undef,nchr)
    correction = Vector{DataFrame}(undef,nchr)
    if isparallel && nprocs()>1
        phasedgenols=[getsubPolyGeno(phasedgeno,chrsubset=[chr]) for chr=1:nchr]
        resmap = pmap(x->chrreconstruct(x,1,
            epsilon,chrpairing,priorspace,correctthreshold,nothing,verbose),phasedgenols)
        for chr=1:nchr
            valentprob[chr],genoprob[chr],correction[chr] = resmap[chr][1]
            iobuffer = resmap[chr][2]
            if io === nothing
                close(iobuffer)
            else
                write(io,String(take!(iobuffer)))
                flush(io)
            end
        end
    else
        for chr=1:nchr
            valentprob[chr],genoprob[chr],correction[chr] = first(chrreconstruct(phasedgeno,chr,
                epsilon,chrpairing,priorspace,correctthreshold,io,verbose))
        end
    end
    changedf = vcat(correction...)
    if size(changedf,1) > 0
        i = size(phasedgeno.correction,1) == 0 ? 1 : phasedgeno.correction[end,:round]+1
        insertcols!(changedf,1,:round =>i*ones(Int, size(changedf,1)))
        phasedgeno.correction = vcat(phasedgeno.correction,changedf)
    end
    polyancestry = PolyAncestry(phasedgeno.markermap,phasedgeno.parentgeno,
        phasedgeno.parentinfo,phasedgeno.offspringinfo,
        phasedgeno.designinfo,phasedgeno.delmarker,phasedgeno.correction,
        statespace,valentprob,genoprob,nothing)
    sethaploprob!(polyancestry)
    # set outliers
    detectOutlier!(polyancestry,tukeyfence=3,minprob=0.6)
    phasedgeno.offspringinfo[!,:isoutlier] = deepcopy(polyancestry.offspringinfo[!,:isoutlier])
    if outstem != nothing
        outlier0=phasedgeno.offspringinfo[!,:isoutlier]
        outlier = skipmissing(outlier0)
        isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
        noutlier == 0 && printconsole(io,verbose,string("no outlier offspring"))
        if noutlier > 0
            offoutid = phasedgeno.offspringinfo[outlier0,:individual]
            msg = string(noutlier, " outlier offspring: ",join(offoutid,", "))
            printconsole(io,verbose,msg)
        end
        if size(phasedgeno.correction,1)>0
            outfile = string(outstem,"_parentphased_corrected.csv")
            savegenodata(outfile,phasedgeno,
                missingstring=missingstring, workdir=workdir)
            msg = string("parentcorrected file: ", outfile)
            printconsole(io,verbose,msg)
        end
        outfile =string(outstem,"_polyancestry.csv")
        savePolyAncestry(outfile, polyancestry,
            missingstring=missingstring, workdir=workdir)
        msg = string("polyancestry file: ", outfile)
        printconsole(io,verbose,msg)
        # genoprob file is a concise version
        outfile =string(outstem,"_genoprob.csv")
        savegenoprob(outfile, polyancestry,
            missingstring=missingstring, workdir=workdir)
        msg = string("genoprob file: ", outfile, ", a concise polyancestry file")
        printconsole(io,verbose,msg)
        outfile = string(outstem,"_postdoseprob.csv")
        doseprob = [calchrdoseprob(genoprob[chr],
            phasedgeno.parentgeno[chr],phasedgeno.offspringinfo,
            phasedgeno.designinfo,priorspace) for chr=1:nchr]
        doseprob = PolyGeno(phasedgeno.markermap,phasedgeno.parentgeno,
            doseprob, phasedgeno.parentinfo,
            phasedgeno.offspringinfo,phasedgeno.designinfo,
            phasedgeno.delmarker,phasedgeno.correction)
        savegenodata(outfile,doseprob,
            missingstring=missingstring, workdir=workdir)
        msg = string("postdoseprob file: ", outfile)
        printconsole(io,verbose,msg)
        if isplot
            probprobdir=saveProbPlot(polyancestry,outstem=outstem,workdir=workdir)
            msg = string("probplot in folder: ", probprobdir)
            printconsole(io,verbose,msg)
        end
    end
    printconsole(io,verbose,string("End, ", Dates.now(), ", time elapsed = ",
        round(time()-starttime), " seconds by polyReconstruct!"))
    if typeof(logfile) <: AbstractString
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    polyancestry
end

function chrreconstruct(phasedgeno::PolyGeno,chr::Integer,epsilon::Real,
    chrpairing::Integer,priorspace::AbstractDict,correctthreshold::Real,
    io::Union{IO,Nothing},verbose::Bool)
    io===nothing && (io=IOBuffer(append=true))
    dict2group = getdict2group(priorspace,phasedgeno)
    priorprocess = getpriorprocess(phasedgeno,chr,chrpairing)
    epsilonls = [epsilon for i=1:size(phasedgeno.markermap[chr],1)]
    bvpair = randinitbvpair(priorspace,phasedgeno) # independent of chr
    idlen=max([length(i[1,:chromosome]) for i=phasedgeno.markermap]...)
    chrid=lpad(phasedgeno.markermap[chr][1,:chromosome],idlen)
    itmax=5
    res=Vector(undef,2)
    correctiondf = DataFrame()
    for it=1:itmax
        # modify epsilonls and bvpair
        inferepsilonls!(chr,epsilonls,bvpair,priorspace,priorprocess,phasedgeno)
        # println(describe(epsilonls))
        decode=chrposteriordecode(chr, epsilonls,dict2group,
            priorspace,priorprocess,phasedgeno,verbose)
        res[1] = [i[1] for i=decode]
        res[2] = [i[2] for i=decode]
        # ["marker", "chromosome","parent", "old_genotype", "new_genotype","old_nerr", "new_nerr"]
        changedf=parentCorrect!(phasedgeno,chr,res[2],priorspace,
            correctthreshold=correctthreshold)
        correctiondf = vcat(correctiondf, changedf)
        # println("chr=", chr, ", correctiondf=",correctiondf)
        nchange= size(changedf,1)
        # println(changedf)
        avgepsilon = round(mean(epsilonls),digits=4)
        msg = string("#chr=", chrid,
            nprocs()>1 ? string(", procs=",myid()) : "",
            ", it=",it, ", <eps>=", avgepsilon)
        if nchange==0
            msg = string(msg, ", no_error_correction")
            printconsole(io,verbose,msg)
            break
        else
            # (:ndosediff, :nphasediff, :nalleleflip, :nerrdiff)
            changesum=getchangesum(changedf)
            msg = string(msg, ", #dose_diff=", changesum.ndosediff,
                ", #phase_diff=",changesum.nphasediff,
                ", #error_diff=",round(changesum.nerrdiff,digits=1))
            printconsole(io,verbose,msg)
        end
    end
    push!(res, correctiondf)
    res, io
end

# posterio decoding

function getfhaplo(parentgeno::AbstractMatrix)
    hcat([vcat(parentgeno[i,:]...) for i=1:size(parentgeno,1)]...)'
end

# organizing outputs as PolyProb

function getstatespace(priorspace::AbstractDict)
    statespace = deepcopy(priorspace)
    for key=keys(statespace)
        for i=["valentkey", "condstate", "state"]
            delete!(statespace[key],i)
        end
    end
    statespace
end

function getdict2group(priorspace::AbstractDict,polygeno::PolyGeno)
    Dict([begin
        condstates = priorspace[popid]["condstate"]
        groupstates = priorspace[popid]["groupstate"]
        states = priorspace[popid]["state"]
        nstate = length(states)
        sortedstates = [sort(i) for i=states]
        tran2group=sparse([i==j for i=sortedstates, j=groupstates])
        popid => tran2group
    end for popid = polygeno.designinfo[!,:population]])
end

function offposteriordecode(dataprob::AbstractMatrix,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvlist::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    loglmax = [-Inf]
    res=[]
    bvkeys = keyls[bvlist]
    keyset = unique(bvkeys)
    for strkey =keyset
        subbvlist = bvlist[bvkeys .== strkey]
        # println("strkey=",strkey, ",#bv=",length(subbvlist))
        prior=priorprocess[strkey]
        if ismissing(prior)
            # startprob,tranprobseq = valentkey2priorproces(strkey,rfseq)
            @error string("todo")
        else
            startprob = prior.startprob
            tranprobseq = Vector{Matrix{Float64}}(prior.tranprobseq[1:end-1])
        end
        for bv=subbvlist
            dataprobseq = [dataprob[i,condstates[bv]] for i=1:size(dataprob,1)]
            fwprob,fwscale = forward(startprob,tranprobseq,dataprobseq)
            logl = calloglike(fwscale)
            if logl- loglmax[1] > -20.0
                bwprob=backward(tranprobseq,dataprobseq,fwscale)
                posteriorprob=[fwprob[i] .* bwprob[i] for i=1:size(dataprobseq,1)]
                push!(res,[bv,logl, posteriorprob])
            end
            loglmax[1] = max(logl, loglmax[1])
            if bv % 100 == 0
                bool = findall([i[2]-loglmax[1] < -20.0 for i=res])
                deleteat!(res,bool)
            end
        end
    end
    bool = findall([i[2]-loglmax[1] < -20.0 for i=res])
    deleteat!(res,bool)
    permutedims(hcat(res...))
end

function getgenoprob(decode::AbstractMatrix,popid::AbstractString,
    priorspace::AbstractDict,tran2group::AbstractMatrix)
    logl = decode[:,2]
    logmax = logsumexp(logl)
    weight = round.(exp.(logl .- logmax),digits=4)
    bool = findall(weight .> 0)
    p = bool[sortperm(weight[bool])]
    logl = round.(logl[p],digits=3)
    weight = weight[p]
    vvpos = decode[p, 1]
    prob = decode[p,3]
    prob2=[weight[i] .* sparse(round.(hcat(prob[i]...),digits=5)') for i=1:length(weight)]
    offres = spzeros(size(prob2[1],1),size(tran2group,1))
    condstates = priorspace[popid]["condstate"]
    for i=1:length(vvpos)
        offres[:,condstates[vvpos[i]]] .+= prob2[i]
    end
    genoprob = round.(offres * tran2group,digits=3)
    vvpos, logl, weight = reverse.([vvpos, logl, weight])
    valentprob = [vvpos logl weight]
    valentprob, genoprob
end

function chrposteriordecode(chr::Integer,epsilon::Union{Real,AbstractVector},
    dict2group::AbstractDict,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    phasedgeno::PolyGeno,verbose::Bool)
    # idlen=max([length(i[1,:chromosome]) for i=phasedgeno.markermap]...)
    # chrid=lpad(phasedgeno.markermap[chr][1,:chromosome],idlen)
    chrdose=phasedgeno.offspringgeno[chr]
    fhaplo= getfhaplo(phasedgeno.parentgeno[chr])
    deriveddose = getderiveddose(fhaplo,priorspace,phasedgeno)
    noff = size(chrdose,2)
    bvprop = getbvpairprop(priorspace,phasedgeno)
    # cd("C:\\Chaozhi\\Workspace\\JuliaWorkspace\\Workspace_Polyploid\\PolyOrigin\\examples\\dihaploid")
    # serialize("temp.test",[bvprop,priorspace,phasedgeno])
    postdecode = [begin
        # msg = string("#chr=", chrid,", haplotype reocnstruction offspring = ",off," out of ",noff)
        # verbose && print("\u1b[K",msg,"\u1b[1G")
        offdose = chrdose[:,off]
        popid = phasedgeno.offspringinfo[off,:population]
        ploidy = phasedgeno.offspringinfo[off,:ploidy]
        tran2group = dict2group[popid]
        dataprob = caldataprob(offdose,popid,ploidy,deriveddose,epsilon)
        decode = offposteriordecode(dataprob, popid, priorspace,priorprocess,bvprop[off])
        valentprob, genoprob = getgenoprob(decode,popid, priorspace,tran2group)
        # valentprob2 = hcat(off*ones(Int,size(valentprob,1),1),valentprob)
        [valentprob, genoprob]
    end for off=1:noff]
    # verbose && print("\u1b[1G\u1b[K")
    postdecode
end

function inferepsilonls!(chr::Integer,epsilonls::AbstractVector,bvpair::AbstractVector,
   priorspace::AbstractDict,priorprocess::AbstractDict,polygeno::PolyGeno;
   byneighbor::Bool=false)
    chrdose=polygeno.offspringgeno[chr]
    fphase = polygeno.parentgeno[chr]
    fhaploset=[[[i] for i=j] for j=eachcol(fphase)]
    fhaploindex=[ones(Int, length(j)) for j=fhaploset]
    siblogl = zeros(length(bvpair))
    popidls = keys(priorspace)
    if byneighbor
        updatebvpair!(bvpair, siblogl,fhaploindex,fhaploset,epsilonls,chrdose,
            priorspace,priorprocess,polygeno,popidls)
    else
        randbvpair!(bvpair, siblogl,fhaploindex,fhaploset,epsilonls,chrdose,
            priorspace,priorprocess,polygeno,popidls,isrand=false)
    end
    snporder = 1:size(chrdose,1)
    dataprobset = caldataprobset(fhaploindex,fhaploset,epsilonls,chrdose,priorspace,polygeno)
    updateepsilonls!(epsilonls, snporder,dataprobset,bvpair,
        priorspace, priorprocess,polygeno,chrindex=chr)
end


function ancestrycall(condprob::AbstractVector;minprob::Real= 0.0)
    [begin
        chrbest = [[begin
            val,index0 = findmax(i)
            val > minprob ? index0 : missing
        end for i=eachrow(j)] for j=condprob[ch]]
        hcat(chrbest...)
    end for ch=1:length(condprob)]
end

function calnumrecom(offancestry::AbstractVector)
    nrecom=[[length(PolyOrigin.splitindex(collect(skipmissing(i))))-1 for i=eachcol(a)]
            for a=offancestry]
    nrecom2 = hcat(nrecom...)
    sum(nrecom2,dims=2)[:,1]
end

function detectOutlier!(polyancestry::PolyAncestry;
    minprob::Real= 0.6,tukeyfence::Real=3)
    bestancestry=ancestrycall(polyancestry.genoprob,minprob=minprob)
    nrecom=calnumrecom(bestancestry)
    anscombe = 2 * (sqrt.(nrecom .+ 3.0/8))
    q1,q2,q3=quantile(anscombe,[0.25,0.5,0.75])
    upbound=q3+tukeyfence*(q3-q1)
    outlier=anscombe .> upbound
    polyancestry.offspringinfo[!,:isoutlier] = outlier
end
