"""
    polyReconstruct!(polygeno::PolyGeno, keyargs...)

performs ancestral inference from polygeno and return polyancestry::PolyAncestry.

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`doseerr::Real=0.01`: genotypic error probability for offspring phasing.

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

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`correctthreshold::AbstractFloat=0.15`: a candidate marker is selected for
parental error correction if the fraction of offspring genotypic error >= correctthreshold.

`isinfererror::Bool=false`: if true, infer dosage error rate per marker. 

`isplot::Bool=false`: if true, plot haploprob for all offspring and save in
the folder "outstem_plots".

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyReconstruct!(phasedgeno::PolyGeno;
    doseerr::Real=0.01,seqerr::Real=0.001,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    isparallel::Bool=true,
    isparalleloffspring::Bool=true, 
    correctthreshold::AbstractFloat=0.15,
    byneighbor::Union{Nothing,Bool}=nothing,
    isinfererror::Bool=false, 
    missingstring::AbstractString="NA",
    isplot::Bool=false,
    nplot_subpop::Integer=10, 
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
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
    printconsole(io,verbose,string("polyReconstruct!, ", Dates.now()))
    if !isnothing(chrsubset)
        verbose && @info string("chrsubset = ",chrsubset)
    end
    if !isnothing(snpsubset)
        verbose && @info string("snpsubset = ",snpsubset)
    end
    isparallel,isparalleloffspring= reset_parallel(isparallel,isparalleloffspring, 
        isinfererror,length(phasedgeno.markermap); io, verbose)
    msg = string("list of option values: \n",
        "doseerr = ", doseerr, "\n",
        "seqerr = ",seqerr,"\n",
        "chrpairing = ",chrpairing, "\n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",
        "isparallel = ", isparallel, "\n",
        "isparalleloffspring = ", isparalleloffspring, "\n",
        "correctthreshold = ", correctthreshold, "\n",
        "byneighbor = ", byneighbor, "\n",
        "isinfererror = ", isinfererror, "\n",
        "missingstring = ", missingstring, "\n",
        "isplot = ", isplot, "\n",
        "nplot_subpop = ", nplot_subpop, "\n",
        "outstem = ", isnothing(outstem) ? "not save to files" : outstem,"\n",
        "logfile = ",logfile,"\n",
        "workdir = ", workdir,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
    if isnothing(byneighbor)
        byneighbor = max(phasedgeno.parentinfo[!,:ploidy]...) >=6
        verbose && @info string("set byneighbor = ",byneighbor)
    end
    msg = string("data: #pop=", size(phasedgeno.designinfo,1),
        ", #parent=",size(phasedgeno.parentinfo,1),
        ", #offspring=",size(phasedgeno.offspringinfo,1),
        ", #chr=",length(phasedgeno.markermap),
        ", #marker=",sum(size.(phasedgeno.markermap,1)))
    printconsole(io,verbose,msg)
    priorspace = getpriorstatespace(phasedgeno,chrpairing)
    statespace = getstatespace(priorspace)
    nchr = length(phasedgeno.markermap)
    genoprob = Vector{Vector}(undef,nchr)
    valentprob = Vector{Vector}(undef,nchr)
    correction = Vector{DataFrame}(undef,nchr)
    tused = @elapsed if isparallel && nprocs()>1
        phasedgenols=[getsubPolyGeno(phasedgeno,chrsubset=[chr]) for chr=1:nchr]
        resmap = pmap(x->chrreconstruct(x,1,doseerr,chrpairing,
            priorspace, correctthreshold,byneighbor,isinfererror, 
            false, nothing,verbose),phasedgenols)
        @everywhere GC.gc()
        for chr=1:nchr
            valentprob[chr],genoprob[chr],correction[chr] = resmap[chr][1]
            iobuffer = resmap[chr][2]
            if isnothing(io)
                close(iobuffer)
            else
                write(io,String(take!(iobuffer)))
                flush(io)
            end
        end
    else
        for chr=1:nchr
            valentprob[chr],genoprob[chr],correction[chr] = first(chrreconstruct(phasedgeno,chr,
                doseerr,chrpairing,
                priorspace,correctthreshold,byneighbor, isinfererror, 
                isparalleloffspring, io,verbose))
        end
    end
    msg = string("tused=",round(tused,digits=1), "s for all ", nchr, " chromosomes")
    printconsole(io, verbose, msg)
    # printconsole(io,verbose,"here0_changedf")
    changedf = vcat(correction...)
    if size(changedf,1) > 0
        i = size(phasedgeno.correction,1) == 0 ? 1 : phasedgeno.correction[end,:round]+1
        insertcols!(changedf,1,:round =>i*ones(Int, size(changedf,1)))
        phasedgeno.correction = vcat(phasedgeno.correction,changedf)
    end
    # printconsole(io,verbose,"here1_polyancestry")
    polyancestry = PolyAncestry(phasedgeno.markermap,phasedgeno.parentgeno,
        phasedgeno.parentinfo,phasedgeno.offspringinfo,
        phasedgeno.designinfo,phasedgeno.delmarker,phasedgeno.correction,
        statespace,valentprob,genoprob,nothing)
    sethaploprob!(polyancestry)
    # set outliers
    # printconsole(io,verbose,"here2_outlier")
    detectOutlier!(polyancestry,tukeyfence=3,minprob=0.5)
    phasedgeno.offspringinfo[!,:isoutlier] = deepcopy(polyancestry.offspringinfo[!,:isoutlier])    
    if !isnothing(outstem)
        # outliers
        outlier0=phasedgeno.offspringinfo[!,:isoutlier]
        outlier = skipmissing(outlier0)
        isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
        noutlier == 0 && printconsole(io,verbose,string("no outlier offspring"))
        if noutlier > 0
            offoutid = phasedgeno.offspringinfo[outlier0,:individual]
            msg = string(noutlier, " outlier offspring: ",join(offoutid,", "))
            printconsole(io,verbose,msg)
        end                
        # visualize results
        plot_valent_DR_recom!(polyancestry; 
            isplot, minprob=0.5, fontsize=14, chrpairing, 
            workdir, outstem, io,verbose
        )        
        # save correction 
        if size(phasedgeno.correction,1)>0
            outfile = string(outstem,"_parentphased_corrected.csv")
            savegenodata(outfile,phasedgeno,
                missingstring=missingstring, workdir=workdir)
            msg = string("parentcorrected file: ", outfile)
            printconsole(io,verbose,msg)
        end

        # printconsole(io,verbose,"here3_save")
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
            probprobdir=saveProbPlot(polyancestry; nplot_subpop, outstem,workdir)
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

function chrreconstruct(phasedgeno::PolyGeno,chr::Integer,doseerr::Real,
    chrpairing::Integer,priorspace::AbstractDict, correctthreshold::Real,
    byneighbor::Bool,isinfererror::Bool,
    isparalleloffspring::Bool, io::Union{IO,Nothing},verbose::Bool)
    isnothing(io) && (io=IOBuffer(append=true))        
    maxvalent = max(digits(chrpairing)...)
    priorprocess = getpriorprocess(phasedgeno,chr,maxvalent)
    dict2group = getdict2group(priorspace,phasedgeno)    
    nmarker = size(phasedgeno.markermap[chr],1)
    doseerrls = [doseerr for i=1:nmarker]
    if isinfererror
        bvpair = randinitbvpair(priorspace,phasedgeno) # independent of chr
    end
    idlen=max([length(i[1,:chromosome]) for i=phasedgeno.markermap]...)
    chrid=lpad(phasedgeno.markermap[chr][1,:chromosome],idlen)
    itmax=20
    res=Vector(undef,2)
    correctiondf = DataFrame()    
    for it=1:itmax
        startt = time()        
        if isinfererror
            # modify doseerrls and bvpair
            inferdoseerrls!(doseerrls,bvpair,priorspace,priorprocess,byneighbor, 
                chr, phasedgeno)
        end
        # println(describe(doseerrls))
        decode=chrposteriordecode(chr, doseerrls,dict2group,
            priorspace,priorprocess,phasedgeno,isparalleloffspring,verbose)
        res[1] = [i[1] for i=decode]
        res[2] = [i[2] for i=decode]
        # printconsole(io,verbose,string("it=",it," after decode"))
        # ["marker", "chromosome","parent", "old_genotype", "new_genotype","old_nerr", "new_nerr"]
        changedf=parentCorrect!(phasedgeno,chr,res[2],priorspace,
            correctthreshold=correctthreshold)
        # printconsole(io,verbose,string("it=",it," after correction"))
        correctiondf = vcat(correctiondf, changedf)
        # println("chr=", chr, ", correctiondf=",correctiondf)
        nchange= size(changedf,1)
        # println(changedf)
        avgdoseerr = round(mean(doseerrls),digits=4)
        msg = string("chr=", chrid, ", #markers=",nmarker, 
            ", it=",it, ", <err>=", avgdoseerr)
        msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
        if nchange==0
            msg = string(msg, msgtused, ", no_error_correction")
            printconsole(io,verbose,msg)
            break
        else
            # (:ndosediff, :nphasediff, :nalleleflip, :nerrdiff)
            changesum=getchangesum(changedf)
            msg = string(msg, ", #dose_diff=", changesum.ndosediff,
                ", #phase_diff=",changesum.nphasediff,
                ", #error_diff_progeny=",round(changesum.nerrdiff,digits=1),
                msgtused
                )
            printconsole(io,verbose,msg)
        end
    end
    push!(res, correctiondf)
    res, io
end

# posterio decoding
function getfhaplo(parentgeno::AbstractMatrix)
    # parentgeno is phased genotypes without missing
    reduce(hcat,[Vector{Int8}(vcat(parentgeno[i,:]...)) for i=1:size(parentgeno,1)])'
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

function offposteriordecode(dataprob::AbstractVector,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvlist::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    loglmax = -Inf
    res=[]
    bvkeys = keyls[bvlist]
    keyset = unique(bvkeys)
    mindiff = log(1e-4)
    for strkey =keyset
        subbvlist = bvlist[bvkeys .== strkey]
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, strkey)
        for bv=subbvlist
            dataprobseq = [i[condstates[bv]] for i = dataprob[markerincl]]
            fwprob,fwscale = forward(startprob,tranprobseq,dataprobseq)
            logl = calloglike(fwscale)
            if logl- loglmax >= mindiff
                bwprob=backward(tranprobseq,dataprobseq,fwscale)
                posteriorprob=[fwprob[i] .* bwprob[i] for i=1:size(dataprobseq,1)]
                push!(res,[bv,logl, posteriorprob])
            end
            loglmax = max(logl, loglmax)
            if bv % 100 == 0
                bool = findall([i[2]-loglmax< mindiff for i=res])
                deleteat!(res,bool)
            end
        end
    end
    bool = findall([i[2]-loglmax < mindiff for i=res])
    deleteat!(res,bool)
    permutedims(reduce(hcat,res))
end

function getgenoprob(decode::AbstractMatrix,popid::AbstractString,
    priorspace::AbstractDict,tran2group::AbstractMatrix)
    logl = decode[:,2]
    logmax = logsumexp(logl)
    weight = exp.(logl .- logmax)
    maxw = maximum(weight)
    weight[weight .< 0.01*maxw] .= 0.0    
    weight = round.(weight ./ sum(weight), digits=3)
    bool = findall(weight .> 0)
    p = bool[sortperm(weight[bool])]
    logl = round.(logl[p],digits=3)
    weight = weight[p]
    vvpos = decode[p, 1]
    prob = decode[p,3]
    # prob2=[weight[i] .* sparse(round.(hcat(prob[i]...),digits=5)') for i=1:length(weight)]
    prob2=[weight[i] .* reduce(hcat,prob[i])' for i=1:length(weight)]
    offres = zeros(size(prob2[1],1),size(tran2group,1))
    condstates = priorspace[popid]["condstate"]
    for i=1:length(vvpos)
        offres[:,condstates[vvpos[i]]] .+= prob2[i]
    end
    genoprob = offres * tran2group
    vvpos, logl, weight = reverse.([vvpos, logl, weight])
    valentprob = [vvpos logl weight]
    valentprob, sparse(round.(genoprob,digits=3))
end

function chrposteriordecode(chr::Integer,doseerrls::AbstractVector,
    dict2group::AbstractDict,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    phasedgeno::PolyGeno,
    isparalleloffspring::Bool, verbose::Bool)
    # idlen=max([length(i[1,:chromosome]) for i=phasedgeno.markermap]...)
    # chrid=lpad(phasedgeno.markermap[chr][1,:chromosome],idlen)
    chrdose=phasedgeno.offspringgeno[chr]
    fhaplo= getfhaplo(phasedgeno.parentgeno[chr])    
    deriveddose = getderiveddose(fhaplo,priorspace,phasedgeno)
    bvprop = getbvpairprop(priorspace,phasedgeno)    
    noff = size(chrdose,2)    
    progress = Progress(noff;
        enabled = verbose,
        desc=string("Posterior decoding..."))
    mymap = isparalleloffspring ? progress_pmap : progress_map
    offinfo = phasedgeno.offspringinfo
    postdecode = mymap((x,y,z)->decode_per_offspring(x, y,z, offinfo, dict2group, 
        doseerrls, deriveddose, priorspace, priorprocess), 
        1:noff, bvprop, eachcol(chrdose);
        progress=progress
    ) 
    postdecode
end

function decode_per_offspring(off::Integer, offbvprop::AbstractVector, offdose::AbstractVector,     
    offinfo::AbstractDataFrame, dict2group::AbstractDict,
    doseerrls::AbstractVector, deriveddose::AbstractDict, 
    priorspace::AbstractDict,priorprocess::AbstractDict)    
    popid = offinfo[off,:population]
    ploidy = offinfo[off,:ploidy]
    tran2group = dict2group[popid]
    dataprob = caldataprob(offdose,ploidy,deriveddose[popid],doseerrls)
    decode = offposteriordecode(dataprob, popid, priorspace,priorprocess,offbvprop)
    valentprob, genoprob = getgenoprob(decode,popid, priorspace,tran2group)    
    [valentprob, genoprob]
end

function inferdoseerrls!(doseerrls::AbstractVector,bvpair::AbstractVector,
   priorspace::AbstractDict,priorprocess::AbstractDict, byneighbor::Bool,
   chr::Integer, polygeno::PolyGeno)
    chrdose=polygeno.offspringgeno[chr]
    fphase = polygeno.parentgeno[chr]
    fhaploset=[[[i] for i=j] for j=eachcol(fphase)]
    fhaploindex=[ones(Int, length(j)) for j=fhaploset]
    siblogl = repeat([-Inf],length(bvpair))
    popidls = keys(priorspace)
    if byneighbor
        updatebvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,
            priorspace,priorprocess,polygeno,popidls,show_progress=true)
    else
        randbvpair!(bvpair, siblogl,fhaploindex,fhaploset,doseerrls,chrdose,
            priorspace,priorprocess,polygeno,popidls,isrand=false)
    end    
    snporder = 1:size(chrdose,1)
    updatedoseerrls!(doseerrls, snporder,bvpair,priorspace, priorprocess,
        chr, polygeno)   
    nothing 
end


function ancestrycall(condprob::AbstractVector; minprob::Real= 0.5)
    [begin
        chrbest = [[begin
            val,index0 = findmax(i)
            val > minprob ? index0 : missing
        end for i=eachrow(j)] for j=condprob[ch]]
        reduce(hcat,chrbest)
    end for ch=1:length(condprob)]
end

function calnumrecom(chrancestry::AbstractMatrix)
    [begin 
        ls = collect(skipmissing(i))
        isempty(ls) ? 0 : length(splitindex(ls))-1             
    end for i=eachcol(chrancestry)]
end

function detectOutlier!(polyancestry::PolyAncestry;
    minprob::Real= 0.5,tukeyfence::Real=3)
    bestancestry=ancestrycall(polyancestry.genoprob,minprob=minprob)    
    nrecom = sum(calnumrecom.(bestancestry))
    anscombe = 2 * (sqrt.(nrecom .+ 3.0/8))
    q1,q2,q3=quantile(anscombe,[0.25,0.5,0.75])
    upbound=q3+tukeyfence*(q3-q1)
    outlier=anscombe .> upbound
    polyancestry.offspringinfo[!,:isoutlier] = outlier
    nothing
end
