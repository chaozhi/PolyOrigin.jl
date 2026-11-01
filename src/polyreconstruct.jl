"""
    polyreconstruct!(polygeno::PolyGeno, keyargs...)

performs ancestral inference from polygeno and return polyancestry::PolyAncestry.

# Positional arguments

`polygeno::PolyGeno`: a struct that stores genotypic data and pedigree info.

# Keyword arguments

`doseerror::Real=0.01`: genotypic error probability for offspring phasing.

`baseerror::Real=0.001`: base sequencing error probability for GBS data.

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

`isplot::Bool=true`: if true, plot haploprob for all offspring and save in
the folder "outstem_plots".

'trueorigfile::Union{Nothing,AbstractString} = nothing': true ancestral origin file produced by PolyoriginSim. 

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.
If nothing, no output files.

`logfile::Union{Nothing,AbstractString,IO}= (outstem===nothing ? nothing : string(outstem,".log"))`:
log file or IO for writing log. If nothing, no log file.

`workdir::AbstractString = pwd()`: directory for reading and writing files.

`verbose::Bool=true`: if true, print messages on console.

"""
function polyreconstruct!(phasedgeno::PolyGeno;
    doseerror::Real=0.01,
    baseerror::Real=0.001,
    chrpairing::Integer=44,
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,    
    bindis::Real=1e-3,     
    iscorrectparent::Union{Nothing,Bool}=nothing,     
    isimputefirst::Bool = true, 
    isparallel::Union{Nothing,Bool}=nothing, 
    isparallelind::Union{Nothing,Bool}=nothing,     
    missingstring::AbstractString="NA",    
    isplot::Bool=true,    
    nplot_subpop::Integer=5, 
    tukeyfence::Real=3, 
    samplesize::Integer=1000,
    trueorigfile::Union{Nothing,AbstractString}=nothing, 
    outstem::Union{Nothing,AbstractString}="outstem",    
    outext::AbstractString = ".vcf.gz", 
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    workdir::AbstractString = pwd(),
    tempdirectory::AbstractString = tempdir(),
    verbose::Bool=true)    
    starttime = time()    
    io = set_logfile_begin(logfile, workdir, "polyreconstruct!"; verbose)     
    if isimputefirst 
        msg_imputefirst = " (reconstruct using imputed offspring)"
    else 
        msg_imputefirst = " (reconstruct before imputing offspring)"
    end
    msg = string("list of option values: \n",
        "doseerror = ", doseerror, "\n",
        "baseerror = ",baseerror,"\n",
        "chrpairing = ",chrpairing, "\n",
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",             
        "bindis = ", bindis, "cM\n",           
        "iscorrectparent = ", iscorrectparent, "\n",
        "isimputefirst = ", isimputefirst, msg_imputefirst, "\n",
        "isparallel = ", isparallel, "\n",
        "isparallelind = ", isparallelind, "\n",        
        "missingstring = ", missingstring, "\n",
        "isplot = ", isplot, "\n",
        "nplot_subpop = ", nplot_subpop, "\n",
        "tukeyfence = ", tukeyfence, "\n",
        "samplesize = ", samplesize, "\n",
        "trueorigfile = ", trueorigfile, "\n",
        "outstem = ", isnothing(outstem) ? "not save to files" : outstem,"\n",
        "outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "workdir = ", workdir,"\n",
        "tempdirectory = ", tempdirectory,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
    info_polygeno(phasedgeno; io, verbose)    
    isnothing(chrsubset) || printconsole(io, verbose, string("chrsubset=",chrsubset))    
    isnothing(snpsubset) || printconsole(io, verbose, string("snpsubset=",snpsubset))    
    getsubpolygeno!(phasedgeno,chrsubset=chrsubset,snpsubset=snpsubset)    
    if isnothing(isparallel) || isnothing(isparallelind) 
        if nprocs()>1
            nmarker = sum([size(i,1) for i in phasedgeno.markermap])            
            nchr = length(phasedgeno.markermap)
            if nchr == 1
                isparallel = false                
            else
                isparallel = isnothing(isparallel) ? (nmarker <= 10^4) : isparallel                 
            end
            isparallelind = !isparallel
        else
            isparallel = false
            isparallelind = false
        end
        printconsole(io,verbose,string("reset isparallel=",isparallel))
        printconsole(io,verbose,string("reset isparallelind=",isparallelind))
    end
    if isnothing(iscorrectparent) 
        iscorrectparent = true
        printconsole(io, verbose, string("reset iscorrectparent = ", iscorrectparent))
    end    

    # rawDoseCall!(phasedgeno,baseerror=baseerror)    
    priorspace = get_priorstatespace(phasedgeno.polyped,chrpairing)
    statespace = getstatespace(priorspace)
    nchr = length(phasedgeno.markermap)    
    outfilels = ["" for _ in 1:nchr]
    jltempdir = mktempdir(tempdirectory; prefix="jl_polyreconstruct_", cleanup=true)	
    tused = @elapsed if isparallel && nprocs()>1
        phasedgenols=[getsubpolygeno(phasedgeno,chrsubset=[chr]) for chr=1:nchr]
        resconstructls = pmap(x->chrreconstruct!(x,1; bindis, 
            chrpairing,doseerror, baseerror, iscorrectparent, isimputefirst, 
            isparallelind=false,jltempdir, io=nothing,verbose),phasedgenols)   # phasedgenols is not modified when using pmap, 
        @everywhere GC.gc()                
        for chr=1:nchr
            newgeno = resconstructls[chr][3]
            phasedgeno.parentgeno[chr] .= newgeno.parentgeno[1]
            phasedgeno.offspringgeno[chr] .= newgeno.offspringgeno[1]
            phasedgeno.markermap[chr] = newgeno.markermap[1]
            outfilels[chr] = resconstructls[chr][1]
            iobuffer = resconstructls[chr][2]            
            if isnothing(io)
                close(iobuffer)
            else
                write(io,String(take!(iobuffer)))
                flush(io)
            end
        end        
    else
        for chr=1:nchr
            outfilels[chr] = first(chrreconstruct!(phasedgeno,chr; bindis, 
                chrpairing,doseerror, baseerror, iscorrectparent, isimputefirst,  
                isparallelind, jltempdir, io,verbose))
            @everywhere GC.gc()
        end
    end    
    msg = string("tused=",round(tused,digits=1), "s for all ", nchr, " chromosomes")
    printconsole(io, verbose, msg)

    
    valentprob = Vector{Vector}(undef,nchr)
    diploprob = Vector{Vector}(undef,nchr)        
    viterbiloglike = Vector{Vector}(undef,nchr)
    viterbipath = Vector{Matrix{Int64}}(undef,nchr)
    parentvalentprob = Vector{Any}(undef,nchr)
    parentdiploprob = Vector{Any}(undef,nchr)
    parentviterbiloglike = Vector{Any}(undef,nchr)    
    parentviterbipath = Vector{Any}(undef,nchr)
    correction = Vector{DataFrame}(undef,nchr)
    chridls = [i[1,:linkagegroup] for i in phasedgeno.markermap]
    for chr in eachindex(chridls)
        jldopen(outfilels[chr]) do file
            # res: offspring_valentprob, offspring_diploprob, offspring_viterbi_loglike, ofspring_viterbi_path, 
            #      correctiondf, parenyt_valentprob, parent_diploprob, parent_viterbi_loglike, parent_viterbi_path, 
            res = file[string(chridls[chr])]
            valentprob[chr],diploprob[chr], viterbiloglike[chr], vpath, correction[chr] = res[1:5]            
            viterbipath[chr] = reduce(hcat,vpath)
            chrnmarker = size(viterbipath[chr],1)
            parentvalentprob[chr],parentdiploprob[chr], parentviterbiloglike[chr],  parentvpath = res[6:9]                   
            parentviterbipath[chr] = isnothing(parentvpath) ? nothing : reduce(hcat, [isnothing(i) ? repeat([missing],chrnmarker) : i for i in parentvpath])            
        end
    end        
    if all(isnothing.(parentvalentprob)) 
        parentvalentprob = nothing
    end
    if all(isnothing.(parentdiploprob)) 
        parentdiploprob = nothing
    end    
    changedf = reduce(vcat, correction)
    if size(changedf,1) > 0
        if haskey(phasedgeno.misc,"correction")
            i = phasedgeno.misc["correction"][end,:run]+1
        else
            i = 1
        end
        insertcols!(changedf,1,:run =>i*ones(Int, size(changedf,1)))
        pushmisc!(phasedgeno,"correction"=>changedf)
    end    
    if !isnothing(outstem)
        # save  offspring_imputed genofile    
        outfile = string(outstem,"_imputed"*outext)         
        savegenodata(outfile,phasedgeno,
            missingstring=missingstring, workdir=workdir)
        msg = string("imputed genofile: ", outfile)
        printconsole(io,verbose,msg)   
    end    
    set_phased_offspring!(phasedgeno, viterbipath, statespace) # offspring geno from GP to GT_phased
    polyancestry = PolyAncestry(phasedgeno.polyped, phasedgeno.markermap,
        phasedgeno.parentgeno,phasedgeno.offspringgeno, 
        statespace, valentprob, parentvalentprob,
        diploprob,nothing, nothing, parentdiploprob, nothing, nothing, 
        viterbipath, parentviterbipath, phasedgeno.misc)
    setgenoprob!(polyancestry)
    sethaploprob!(polyancestry)
    # set outliers    
    detectoutlier!(polyancestry; isviterbi=true, tukeyfence)
    phasedgeno.polyped.offspringinfo[!,:isoutlier] = deepcopy(polyancestry.polyped.offspringinfo[!,:isoutlier])


    if !isnothing(outstem)        
        outlier0=phasedgeno.polyped.offspringinfo[!,:isoutlier]
        outlier = skipmissing(outlier0)
        isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
        noutlier == 0 && printconsole(io,verbose,string("no outlier offspring"))
        if noutlier > 0
            offoutid = phasedgeno.polyped.offspringinfo[outlier0,:individual]
            msg = string(noutlier, " outlier offspring: ",join(offoutid,", "))
            @warn msg
            printconsole(io,false,"Warn: "*msg)
        end        

        # save  offspring_phased genofile    
        outfile = string(outstem,"_phased"*outext)           
        savegenodata(outfile,phasedgeno,
            missingstring=missingstring, workdir=workdir)
        msg = string("phased genofile: ", outfile)
        printconsole(io,verbose,msg)   
        
        # visualize valent results        
        plot_valent_DR_recom!(polyancestry; 
            isplot, isviterbi=true, minprob=0.5, fontsize=14,  
            workdir, outstem, io,verbose
        )        
        # save ancestry
        outfile =string(outstem,"_polyancestry.csv.gz")
        savepolyancestry(outfile, polyancestry,
            missingstring=missingstring, workdir=workdir)
        msg = string("polyancestry file: ", outfile)
        printconsole(io,verbose,msg)
    end        
    msg = "it might take a while to produce joinedancestry!"
    printconsole(io,verbose,msg)
    tused = @elapsed joinedancestry = JoinedPolyAncestry(polyancestry; samplesize = samplesize, isparallel);         
    printconsole(io, verbose, string("tused = ", round(tused,digits=1), "s"))
    if !isnothing(outstem)
        outfile =string(outstem,"_joinedancestry.csv.gz")    
        try 
            savejoinedancestry(outfile, joinedancestry,
                missingstring=missingstring, workdir=workdir)
            msg = string("joinedancestry file: ", outfile)
            printconsole(io,verbose,msg)
        catch err
            msg = string("saving joinedancestry failed! ", err)
            @error msg
            printconsole(io,false,"Error: "*msg)                                    
        end
        # genoprob file is a concise version
        # outfile =string(outstem,"_genoprob.csv.gz")
        # savegenoprob(outfile, polyancestry,
        #     missingstring=missingstring, workdir=workdir)
        # msg = string("genoprob file: ", outfile, ", a concise polyancestry file")
        # printconsole(io,verbose,msg)       
        if isplot        
            if isnothing(trueorigfile)
                trueorig = nothing
            else                
                trueorig = readtruegeno!(trueorigfile,joinedancestry; workdir,verbose=false)    
            end
            for isviterbi  in [false]
                saveprobplot(joinedancestry; trueorig, isviterbi, nplot_subpop, outstem,workdir)                            
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    set_logfile_end(logfile, io, tused,"polyreconstruct!"; verbose)
    polyancestry, joinedancestry
end

function reset_errorls(errls::AbstractVector, errrate::Real)
    [ismissing(i) ? errrate : i for i in errls]    
end

function set_phased_offspring!(phasedgeno::PolyGeno, viterbipath::AbstractVector, statespace::AbstractDict)
    nchr = length(phasedgeno.markermap) 
    for chr = 1:nchr        
        offgeno = phasedgeno.offspringgeno[chr]
        parentgeno = phasedgeno.parentgeno[chr]
        chrvpath = viterbipath[chr]
        for (fam,famspace) in statespace
            fam_pgeno = parentgeno[:,famspace["parentindex"]]
            fam_phaplo = eachcol(reduce(hcat,[reduce(vcat,i) for i in eachrow(fam_pgeno)]))            
            statels = famspace["state"]
            offls = famspace["offspringindex"]
            for off =  offls                
                pathls = statels[chrvpath[:,off]]
                offgeno[:,off] .= map((x,y)->x[y], fam_phaplo, pathls)
            end
        end
        df = phasedgeno.markermap[chr]
        df[!,"offspringformat"] = ["GT_phased" for _ in 1:size(df,1)]
    end
    phasedgeno
end

# offfspring are imputed into GP
function chrreconstruct!(phasedgeno::PolyGeno,chr::Integer;    
    chrpairing::Integer, 
    bindis::Real,   
    doseerror::Real,
    baseerror::Real,        
    iscorrectparent::Bool, 
    isimputefirst::Bool, 
    isparallelind::Bool, 
    jltempdir::AbstractString, 
    io::Union{IO,Nothing},
    verbose::Bool)
    io===nothing && (io=IOBuffer(append=true))
    nmarker = size(phasedgeno.markermap[chr],1)
    doseerrorls = reset_errorls(phasedgeno.markermap[chr][!,:doseerror],doseerror)
    baseerrorls = reset_errorls(phasedgeno.markermap[chr][!,:baseerror],baseerror)    
    allelicbiasls = reset_errorls(phasedgeno.markermap[chr][!,:allelicbias],1.0)    
    overdispersionls = reset_errorls(phasedgeno.markermap[chr][!,:overdispersion],0.0)
    doseerrorls = allequal(doseerrorls) ? doseerrorls[1] : doseerrorls
    baseerrorls = allequal(baseerrorls) ? baseerrorls[1] : baseerrorls
    allelicbiasls  = allequal(allelicbiasls) ? allelicbiasls[1] : allelicbiasls 
    overdispersionls = allequal(overdispersionls) ? overdispersionls[1] : overdispersionls  

    priorspace = get_priorstatespace(phasedgeno.polyped,chrpairing)
    polyped = phasedgeno.polyped
    # dict2group = getdict2group(priorspace,polyped)    
    dict2group = nothing # NOT transforming diploprob into genoprob
    maxvalent = max(digits(chrpairing)...)
    priorprocess = get_priorprocess(phasedgeno,chr,maxvalent; bindis)

    # idlen=max([length(i[1,:linkagegroup]) for i=phasedgeno.markermap]...)
    chrid = phasedgeno.markermap[chr][1,:linkagegroup]        
    chroffgeno = get_chroffgeno(phasedgeno,chr)      
    chroffformat = phasedgeno.markermap[chr][!,:offspringformat]
    chrmultiallele = phasedgeno.markermap[chr][!,:multiallele]    
    chrparentgeno= getchrparentgeno(phasedgeno.parentgeno[chr])
    derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)
    targetls = ["doseerror"]
    if !all(occursin.("GT", unique(chroffformat))) 
        append!(targetls, ["baseerror", "allelicbias","overdispersion"])
    end
    
    # decoding for leaf-offspring    
    itmax=20
    res=Vector(undef,4)
    correctiondf = DataFrame()    
    for it=1:itmax
        startt = time()        
        msg = string("chr=", chrid, ", #markers=",nmarker)
        msg *= get_msg_err(targetls, priorprocess; 
            doseerror = doseerrorls, baseerror = baseerrorls, 
            allelicbias = allelicbiasls, 
            overdispersion = overdispersionls)        
        chrdecode = chrhmmdecode(chroffgeno, derivedgeno,priorspace, priorprocess, 
            polyped, chroffformat, chrmultiallele; dict2group, 
            doseerror = doseerrorls, baseerror = baseerrorls, 
            allelicbias = allelicbiasls, overdispersion = overdispersionls,
            ismidparent =false, isinclviterbi = true, isparallelind, jltempdir, verbose)
        for k in 1:4
            res[k] = [i[k] for i=chrdecode]
        end                
        if iscorrectparent               
            correct_minerrdiff = min(3,max(1, round(Int,size(chroffgeno,2)/5)))      
            changelogl, changedf = parentcorrect(chrparentgeno, chroffgeno, chroffformat, chrmultiallele, 
                priorspace, priorprocess, polyped; 
                doseerror=doseerrorls, baseerror=baseerrorls, 
                allelicbias=allelicbiasls, overdispersion = overdispersionls, 
                phaseerror_midparent = nothing, bvpair = nothing, callthreshold = 0.7, 
                minerrdiff = correct_minerrdiff, 
            )               
            isaccept = changelogl > 0 
            if isaccept 
                correct_chrphalo!(chrparentgeno, changedf)    
                correctiondf = vcat(correctiondf, changedf)                
            end
            ncorrect = size(changedf,1)                                   
            msg_correct = string(", #correct=", ncorrect)
            ncorrect > 0 && (msg_correct *= isaccept ? "(accepted)" : "(rejected)")
            if ncorrect > 0 && isaccept
                derivedgeno = get_derivedgeno(chrparentgeno,chrmultiallele, priorspace)
            end                                
        else
            isaccept = false
            ncorrect = 0
            msg_correct = ", iscorrectparent=false"
        end         
        if ncorrect == 0 || it == itmax || !isaccept
            phasedgeno.parentgeno[chr] = Matrix{Any}([ismissing(i) ? i : string.(i) for i in chrparentgeno])
            phasedgeno.offspringgeno[chr] .= cal_condprob2doseprob(chrdecode, chrparentgeno, priorspace,priorprocess)                        
            chroffformat .= "GP"               
            if isimputefirst
                # baseerror/allelicbias/od are not requred for GP format 
                chrdecode = chrhmmdecode(phasedgeno.offspringgeno[chr], derivedgeno,priorspace, priorprocess, 
                    polyped, chroffformat, chrmultiallele; dict2group, 
                    doseerror=1e-3, ismidparent =false, isinclviterbi = true, isparallelind, jltempdir, verbose)                            
                for k in 1:4
                    res[k] = [i[k] for i=chrdecode]
                end
            end            
            msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
            msg *= string(msgtused, msg_correct)            
            printconsole(io,verbose,msg)            
            break
        else
            # (:ndosediff, :nphasediff, :nalleleflip, :nerrdiff)
            changesum=getchangesum(changedf)
            msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
            msg *= string(msgtused, msg_correct)            
            msg *= string(", #dose_diff=", changesum.ndosediff,
                ", #phase_diff=",changesum.nphasediff,
                ", #error_diff_progeny=",round(changesum.nerrdiff,digits=1),
                )
            printconsole(io,verbose,msg)            
        end        
    end
    push!(res, correctiondf)

    # decoding for intermidiate parent
    nmidparent = sum(.!polyped.parentinfo[!,:isfounder])
    if nmidparent == 0
        append!(res, [nothing for _ in 1:4])
    else
        startt = time()
        chrformat = repeat(["GT_phased"],size(chrparentgeno,1))        
        chrdecode = chrhmmdecode(chrparentgeno, derivedgeno,priorspace, priorprocess, 
                polyped, chrformat, chrmultiallele; 
                # baseerror/allelicbias/od are not requred for GT_phased
                # doseerror is used to for phasing error rate of GT_phased; too small phased error rates might result in NaN or Inf in posterior decoding                 
                doseerror = 0.05,                                 
                isparallelind = isparallelind && (nmidparent > nprocs() > 1), 
                ismidparent =true, isinclviterbi = true, jltempdir, verbose)    
        for k in 1:4
            push!(res, [isnothing(i) ? i : i[k] for i=chrdecode]) # parent  valentprob
        end        
        msgtused = string(", tused=",round(time() - startt, digits=1),"s" )
        msg = string("chr=", chrid, ", #markers=",nmarker, ", midparent decoding", msgtused)
        printconsole(io, verbose, msg)
    end
    # save results    
    tempid = tempname(jltempdir,cleanup=true)		
    outfile = string(tempid, "_chrreconstruct_",chrid, ".jld2")
    jldopen(outfile,"w") do file
        # res: offspring_valentprob, offspring_genoprob, offspring_viterbi_loglike, ofspring_viterbi_path, 
        #      correctiondf, parenyt_valentprob, parent_genoprob, parent_viterbi_loglike, parent_viterbi_path, 
        file[string(chrid)] = res
    end
    outfile, io,phasedgeno
end

function cal_condprob2doseprob(resdecode::AbstractVector, chrphaplo::AbstractMatrix,
    priorspace::AbstractDict,priorprocess::AbstractDict; 
    snporder::Union{Nothing,AbstractVector}=nothing)
    dict_d2g = get_state2groupstate(priorspace)
    nmarker = size(chrphaplo,1)
    noff = length(resdecode)    
    isnothing(snporder)  && (snporder = 1:nmarker)
    markerincl = first(values(priorprocess)).markerincl # markerincl is the same for all popid
    markerindices = snporder[markerincl]
    doseprob= Matrix{Union{Missing, Vector{Float64}}}(missing,noff,nmarker)
    for popid in keys(priorspace)
        poppriorspace = priorspace[popid]
        groupstates = poppriorspace["groupstate"]
        ploidy = poppriorspace["offspringploidy"]
        offls = poppriorspace["offspringindex"]
        pls = poppriorspace["parentindex"]
        phaplo = reduce(hcat,reduce.(vcat, eachrow(chrphaplo[markerindices,pls])))
        d2g = dict_d2g[popid]
        for off in offls
            diploprob = resdecode[off][2]
            size(diploprob, 1) == length(markerindices) || @error "inconsistent markerindices"
            if size(diploprob, 2) == size(d2g,1)
                genoprob = (diploprob * d2g)'            
            elseif size(diploprob, 2) == size(d2g,2)
                genoprob = diploprob'
            else
                error(string("unexpected size of condprob: ",size(diploprob)))
            end
            nallele = 2 # TODO for multiallelic markers
            nunphased = get_nunphased(ploidy,nallele)
            for t in eachindex(markerindices)                    
                I, V = findnz(genoprob[:,t])                
                prob = zeros(nunphased) 
                for i in eachindex(I, V)
                    d = unphased_genocode_2_index(phaplo[groupstates[I[i]], t])
                    prob[d] += V[i]
                end
                doseprob[off, markerindices[t]] = prob
            end
        end
    end
    permutedims(doseprob)    
end

function get_state2groupstate(priorspace::AbstractDict)    
    Dict([begin        
        groupstates = priorspace[popid]["groupstate"]
        states = priorspace[popid]["state"]        
        sortedstates = [sort(i) for i=states]
        tran2group=sparse([i==j for i=sortedstates, j=groupstates])
        popid => tran2group
    end for popid = keys(priorspace)])
end

function getdict2group(priorspace::AbstractDict,polyped::PolyPed)
    Dict([begin
        condstates = priorspace[popid]["condstate"]
        groupstates = priorspace[popid]["groupstate"]
        states = priorspace[popid]["state"]
        nstate = length(states)
        sortedstates = [sort(i) for i=states]
        tran2group=sparse([i==j for i=sortedstates, j=groupstates])
        popid => tran2group
    end for popid = polyped.designinfo[!,:family]])
end

function chrhmmdecode(chrgeno::AbstractMatrix, 
    derivedgeno::AbstractDict,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed,    
    chrformat::AbstractVector, 
    chrmultiallele::AbstractVector; 
    dict2group::Union{Nothing,AbstractDict} = nothing,
    ismidparent::Bool,     
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}=0.001, 
    allelicbias::Union{Real,AbstractVector}=1.0,        
    overdispersion::Union{Real,AbstractVector}=0.0,
    bvpairprop::Union{Nothing, NamedTuple}=nothing,
    probdigits::Integer=4,
    isinclviterbi::Bool, 
    isparallelind::Bool, 
    jltempdir::AbstractString, 
    verbose::Bool)    
    if isnothing(bvpairprop)
        bvpairprop = get_prop_bvpair(priorspace,polyped)
    end    
    indinfo = ismidparent ? polyped.parentinfo : polyped.offspringinfo    
    if ismidparent                
        indls = findall(.!indinfo[!,:isfounder])
        pgenols = [phased_genocode_2_index.(chrgeno[:,i],chrmultiallele) for i in indls]
        bvpropls =  bvpairprop.parent[indls]
    else
        indls = collect(1:size(indinfo,1))
        bvpropls =  bvpairprop.offspring
        pgenols = eachcol(chrgeno)
    end    
    nmarker = size(chrgeno,1)    
    doseerrorls = isa(doseerror, Real) ? doseerror * ones(nmarker) : doseerror
    baseerrorls = isa(baseerror, Real) ? baseerror * ones(nmarker) : baseerror
    allelicbiasls = isa(allelicbias, Real) ? allelicbias * ones(nmarker) : allelicbias
    overdispersionls = isa(overdispersion, Real) ? overdispersion * ones(nmarker) : overdispersion
    if isparallelind                
        progress = Progress(length(indls);
            enabled = verbose,    
            desc=string(ismidparent ? "midparent" : "offspring", " decoding...")
        )   
        tempid = tempname(jltempdir,cleanup=true)		        
        outfilels = [string(tempid, "_hmmdecode_", ismidparent ? "parent" : "offspring", i, ".jld2") for i in indls]    
        resls = try 
            progress_pmap((x,y,z,w)->decode_per_ind(x, y,z, indinfo, 
                chrformat, chrmultiallele, doseerrorls, baseerrorls, allelicbiasls, overdispersionls,
                derivedgeno, priorspace, priorprocess; probdigits, dict2group, outfile = w,isinclviterbi), 
                indls, bvpropls, pgenols,outfilels;
                progress=progress
            ) 
            [jldopen(outfile, "r") do file
                [file["valentprob"],file["genoprob"], file["viterbi_loglike"], file["viterbi_path"]]
            end for outfile in outfilels]        
        finally
            rm.(outfilels; force=true)
        end        
    else        
        resls = map((x,y,z)->decode_per_ind(x, y,z, indinfo, 
            chrformat, chrmultiallele, doseerrorls, baseerrorls, allelicbiasls, overdispersionls,
            derivedgeno, priorspace, priorprocess; probdigits,dict2group, isinclviterbi), 
            indls, bvpropls, pgenols
        )     
        
    end
    if ismidparent
        resls2 = Vector{Any}(nothing, size(indinfo,1)) 
        resls2[indls] .= resls
        resls2
    else
        resls
    end
end

function decode_per_ind(ind::Integer, indbvprop::Union{Missing, AbstractVector}, 
    indgeno::AbstractVector,     
    indinfo::AbstractDataFrame,     
    chrformat::AbstractVector, 
    chrmultiallele::AbstractVector, 
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector,    
    allelicbiasls::AbstractVector, 
    overdispersionls::AbstractVector,    
    derivedgeno::AbstractDict, 
    priorspace::AbstractDict,
    priorprocess::AbstractDict;
    dict2group::Union{Nothing,AbstractDict} = nothing,
    probdigits::Integer, 
    isinclviterbi::Bool, 
    outfile::Union{Nothing,AbstractString}=nothing)        
    popid = indinfo[ind,:family]
    ploidy = indinfo[ind,:ploidy]    
    tran2group = isnothing(dict2group) ? nothing : dict2group[popid]    
    nstate = length(priorspace[popid]["state"])
    dataprobls = [ones(nstate) for _ in 1:length(indgeno)]        
    cal_ind_dataprobls!(dataprobls, indgeno,chrformat, chrmultiallele, derivedgeno[popid]; 
        ploidy,doseerror=doseerrorls, baseerror=baseerrorls, 
        allelicbias=allelicbiasls, overdispersion=overdispersionls)    
    decode = indposteriordecode(dataprobls, popid, priorspace,priorprocess,indbvprop)    
    valentprob, condprob = getcondprob(decode,popid, priorspace,tran2group; probdigits)  
    if isinclviterbi
        # valentprob columns: vvpos, logl, weight; row is ordered by decreasing logl
        pos = argmax(valentprob[:,2])
        pos == 1 || @warn "unexpected row ordering of valentprob" valentprob
        indbv = valentprob[pos,1]
        viterbi_loglike, viterbi_path = indviterbidecode(dataprobls, popid,priorspace, priorprocess, indbv)
        res = (valentprob, condprob, viterbi_loglike, viterbi_path)
    else
        res = (valentprob, condprob)
    end
    if isnothing(outfile)
        res
    else
        if isinclviterbi
            jldsave(outfile; valentprob = res[1], genoprob=res[2], viterbi_loglike = res[3], viterbi_path = res[4])
        else
            jldsave(outfile; valentprob = res[1], genoprob=res[2])
        end
        outfile
    end
end


function indviterbidecode(dataprobls::AbstractVector,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    indbv::Integer)
    spacespace = priorspace[popid]["condstate"][indbv]    
    strkey = priorspace[popid]["valentkey"][indbv]    
    markerincl, startprob, tranprobseq = strkey2prior(priorprocess, strkey)
    dataprobseq = [i[spacespace] for i = dataprobls[markerincl]]
    resviterbi = viterbidecode(startprob,tranprobseq,dataprobseq)
    (; logprob, path) = resviterbi    
    path .= spacespace[path]
    logprob,path
end

function indposteriordecode(dataprobls::AbstractVector,popid::AbstractString,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvlist::AbstractVector)
    condstates = priorspace[popid]["condstate"]
    keyls = priorspace[popid]["valentkey"]
    loglmax = -Inf
    res=[]
    bvkeys = keyls[bvlist]
    keyset = unique(bvkeys)
    mindiff = log(1e-6) # smaller than 1e-(probdigits)
    for strkey =keyset
        subbvlist = bvlist[bvkeys .== strkey]
        markerincl, startprob, tranprobseq = strkey2prior(priorprocess, strkey)
        for bv=subbvlist
            dataprobseq = [i[condstates[bv]] for i = dataprobls[markerincl]]
            fwprob,fwscale = forward(startprob,tranprobseq,dataprobseq)
            logl = calloglike(fwscale)
            # NaN >= mindiff: false
            if logl- loglmax >= mindiff || isnan(logl)
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
    loglls = [i[2] for i in res]
    bool = [i-loglmax < mindiff for i=loglls]
    if all(bool) || isempty(loglls) || any(isnan.(loglls))
        @warn string("unexpected list of logl = ",loglls, ", loglmax = ",loglmax)
    end
    deleteat!(res,bool)
    permutedims(reduce(hcat,res))
end

function getcondprob(decode::AbstractMatrix,popid::AbstractString,
    priorspace::AbstractDict,tran2group::Union{Nothing,AbstractMatrix}; probdigits::Integer=4)
    logl = decode[:,2]
    logmax = logsumexp(logl)    
    weight = round.(exp.(logl .- logmax), digits=probdigits+1)
    p = sortperm(logl, rev=true)
    logl = round.(logl[p],digits=probdigits+1)
    weight = weight[p]
    vvpos = decode[p, 1]
    prob = decode[p,3]
    # prob2=[weight[i] .* sparse(round.(hcat(prob[i]...),digits=5)') for i=1:length(weight)]
    prob2=[weight[i] .* reduce(hcat,prob[i])' for i= eachindex(weight)]    
    nmarker = size(prob2[1],1)
    nstate = length(priorspace[popid]["state"])
    offres = zeros(nmarker,nstate)
    condstates = priorspace[popid]["condstate"]
    for i=eachindex(vvpos)
        offres[:,condstates[vvpos[i]]] .+= prob2[i]
    end
    condprob = isnothing(tran2group) ? offres : (offres * tran2group)        
    valentprob = [vvpos logl weight]
    valentprob, sparse(round.(condprob,digits=probdigits))
end

# posterio decoding
function get_phaplo(parentgeno::AbstractMatrix)
    # parentgeno is phased genotypes without missing
    reduce(hcat,[Vector{Int8}(vcat(parentgeno[i,:]...)) for i=1:size(parentgeno,1)])'
end

# organizing outputs as PolyProb

function getstatespace(priorspace::AbstractDict)
    statespace = deepcopy(priorspace)    
    statespace
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

function detectoutlier!(polyancestry::PolyAncestry;
    isviterbi::Bool=true,
    minprob::Real= 0.5,
    tukeyfence::Real=3)
    noff = size(polyancestry.polyped.offspringinfo,1)
    noff < 3 && return nothing # too few offspring to detect outlier
    if isviterbi
        bestancestry = polyancestry.viterbipath
    else
        bestancestry=ancestrycall(polyancestry.genoprob,minprob=minprob)    
    end
    nrecom = sum(calnumrecom.(bestancestry))
    anscombe = 2 * (sqrt.(nrecom .+ 3/8))
    q1,q2,q3=quantile(anscombe,[0.25,0.5,0.75])
    # upbound=q3+tukeyfence*(q3-q1)
    # Carling median rule: 
    nind = length(anscombe)
    nind == noff || @error "unexpected length of anscombe: ", length(anscombe), ", expected: ", noff
    ck = (17.63*nind - 23.64)/(7.74*nind-3.71) # correspoinding to tukeyfence=1.5
    ck *= tukeyfence/1.5 
    upbound = q2+ck*(q3-q1)
    # Sigma rule
    # sigma = StatsBase.std(anscombe)
    # mu = StatsBase.mean(anscombe)
    # upbound3 = mu + tukeyfence*1.349*sigma
    outlier= anscombe .> upbound
    polyancestry.polyped.offspringinfo[!,:isoutlier] = outlier
    nothing
end