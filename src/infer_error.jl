function get_outerlier_indices(ls::AbstractVector; 
	tukeyfence::Real=3, 
    side::String="both")    
    sidels = ["lower", "upper","both"]
    in(side, sidels) || @warn string(side, " must be one of ", sidels)
    q1,q2,q3=quantile(ls,[0.25,0.5,0.75])
	if side in ["lower","both"]
        lowbound= q1-tukeyfence*(q3-q1)
		res = findall(ls .< lowbound) 
	else
		res = Int[]
	end
	if side in ["upper","both"]
        upbound=q3+tukeyfence*(q3-q1)
		append!(res,findall(ls .> upbound))
	end
	res
end

function cal_logprior_error(targetls::AbstractVector,
    doseerrorls::Union{Real,AbstractVector}, 
    baseerrorls::Union{Real,AbstractVector}, 
    allelicbiasls::Union{Real,AbstractVector},         
    overdispersionls::Union{Real,AbstractVector},   
    errorprior::NamedTuple; 
    markerls::Union{Nothing,AbstractVector}=nothing)
    res = 0.0
    errlsls = (doseerrorls,allelicbiasls,overdispersionls,baseerrorls)
    isvecls = [isa(i,AbstractVector) for i in errlsls]
    vartargetls = intersect(["doseerror","allelicbias", "overdispersion","baseerror"][isvecls],targetls)
    if isnothing(markerls) && any(isvecls)             
        markerls = eachindex(errlsls[findfirst(isvecls)])                
        nmarker = length(markerls)
    else
        nmarker = 1
    end    
    doseerrorls2,allelicbiasls2,overdispersionls2,baseerrorls2 = [isa(i,AbstractVector) ? i : i*ones(nmarker) for i in errlsls]
    for target in targetls
        if target == "doseerror"          
            prior_err = errorprior.doseerror
        elseif target == "baseerror"
            prior_err = errorprior.baseerror
        elseif target == "allelicbias"        
            prior_err = errorprior.allelicbias
        elseif target == "overdispersion"
            prior_err = errorprior.overdispersion
        else
            @error string("unknown error target: ",target)
        end              
        if in(target, vartargetls)
            for snp in markerls
                res += cal_logprior_errormarker(target, prior_err, doseerrorls2[snp], baseerrorls2[snp],
                    allelicbiasls2[snp],overdispersionls2[snp])
            end
        else
            res += cal_logprior_errormarker(target, prior_err, first(doseerrorls2), first(baseerrorls2),first(allelicbiasls2), first(overdispersionls2))        
        end
    end
    res
end

function cal_logprior_error(targetls::AbstractVector, priorprocess::AbstractDict; 
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector,         
    overdispersionls::AbstractVector, 
    errorprior::NamedTuple, 
    snporder::Union{Nothing,AbstractVector}=nothing)
    markerincl = first(values(priorprocess)).markerincl
    markerls = isnothing(snporder) ? findall(markerincl) : collect(eachindex(markerincl))
    cal_logprior_error(targetls, doseerrorls, baseerrorls, allelicbiasls, overdispersionls, errorprior; markerls)
end


function cal_logprior_errormarker(target::AbstractString, prior_err::Distribution, 
    doseerrormarker::Real, baseerrormarker::Real, allelicbiasmarker::Real,
    overdispersionmarker::Real)
    if target == "doseerror"
        doseerrormarker = max(1e-5,doseerrormarker)
        x = log(doseerrormarker/(1-doseerrormarker))
    elseif target == "baseerror"
        baseerrormarker = max(1e-5,baseerrormarker)
        x = log(baseerrormarker/(1-baseerrormarker))    
    elseif target == "allelicbias"
        x = max(-10.0,log(allelicbiasmarker))    
    elseif target == "overdispersion"
        x = max(-10.0,log(overdispersionmarker))    
    else
        @error string("unknown genoerror type: ",target)
    end       
    logpri, logjacobi = cal_logpri_logjacobi_err(target, prior_err, x)
    logpri + logjacobi
end


function infer_errorls!(chrparentgeno::AbstractMatrix, # input parent geno
    chrparentformat::AbstractVector, 
    phaploindex::AbstractVector,phaploset::AbstractVector,    
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed, 
    targetls::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    errorprior::NamedTuple, 
    errorsoftthresh::NamedTuple,
    tukeyfence::Real, 
    errorhardthresh::NamedTuple,
    isdellargeerr::Bool=true, 
    bindis::Real,
    errorbound::Real=5, 
    errormaxit::Integer=20,
    snporder::Union{Nothing,AbstractVector}=nothing,        
    temperature::Real=0.0)
    nmarker = size(chroffgeno,1)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nmarker) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nmarker) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nmarker) : allelicbias      
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nmarker) : overdispersion
    chrphaplo = get_phaplo(phaploindex,phaploset)
    derivedgeno = get_derivedgeno(chrphaplo,chrmultiallele, priorspace)   
    infer_errorls!(chrparentgeno, chrparentformat, 
        chrphaplo, derivedgeno, chroffgeno, chroffformat, chrmultiallele,     
        bvpair, priorspace, priorprocess, polyped, targetls;     
        phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, errorprior, 
        errorsoftthresh,tukeyfence, errorhardthresh, bindis,
        isdellargeerr, snporder, temperature, errorbound, errormaxit, 
    )    
end

function infer_avgerror(targetls::AbstractVector, 
    phaploindex::AbstractVector,phaploset::AbstractVector,    
    chrparentgeno::AbstractMatrix, # input parent geno
    chrparentformat::AbstractVector, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed;    
    phaseerror_midparent::Union{Nothing,Real},    
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    errorprior::NamedTuple, 
    errorbound::Real=5, 
    errormaxit::Integer=20,
    snporder::Union{Nothing,AbstractVector}=nothing)
    isempty(targetls) && return (doseerror, baseerror, allelicbias, overdispersion)
    phaplo = PolyOrigin.get_phaplo(phaploindex,phaploset)
    derivedgeno = PolyOrigin.get_derivedgeno(phaplo,chrmultiallele, priorspace)    
    for target in targetls        
        est = infer_avgerror(target, chrparentgeno, chrparentformat, phaplo, derivedgeno, 
            chroffgeno, chroffformat, chrmultiallele, 
            bvpair, priorspace, priorprocess,polyped;     
            phaseerror_midparent, snporder, errorprior, errorbound, errormaxit, 
            doseerror, baseerror, allelicbias, overdispersion, 
        )    
        if target == "doseerror"
            doseerror = est
        elseif target == "baseerror"
            baseerror = est
        elseif target == "allelicbias"
            allelicbias = est   
        elseif target == "overdispersion"
            overdispersion = est    
        end
    end
    doseerror, baseerror, allelicbias, overdispersion
end

function infer_avgerror(target::AbstractString, 
    chrparentgeno::AbstractMatrix, # input parent geno
    chrparentformat::AbstractVector, 
    chrphaplo::AbstractMatrix, # estimated phased gneotypes
    derivedgeno::AbstractDict, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed;    
    phaseerror_midparent::Union{Nothing,Real},    
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector},  
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    errorprior::NamedTuple, 
    errorbound::Real=5,
    errormaxit::Integer=20, 
    skeletonsize::Union{Nothing,Integer}=200, 
    snporder::Union{Nothing,AbstractVector}=nothing)
    accuracygoal, precisiongoal, itmax = 3, 3,errormaxit
    chrparent_hiddenunphased = unphased_genocode_2_index.(join.(sort.(chrphaplo),"/"))
    function f(x)
        doseerror2, baseerror2, allelicbias2, overdispersion2 = doseerror, baseerror, allelicbias, overdispersion
        if target == "doseerror"
            doseerror2 = 1/(1+exp(-x))  # inverse of logit transformation     
            prior_err = errorprior.doseerror
        elseif target == "baseerror"            
            baseerror2 = 1/(1+exp(-x))      
            prior_err = errorprior.baseerror                        
        elseif target == "allelicbias"
            allelicbias2 = exp(x)     
            prior_err = errorprior.allelicbias   
        elseif target == "overdispersion"
            overdispersion2 = exp(x)
            prior_err = errorprior.overdispersion
        else
            @error string("unknown genoerror type: ",target)
        end
        logpri, logjacobi = cal_logpri_logjacobi_err(target, prior_err, x)
        # calculate logl for offspring
        loglls = calmarglogl_bvpair(chrphaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
            priorspace, priorprocess, bvpair, polyped; phaseerror_midparent, snporder, skeletonsize, 
            exclparent = true,  
            doseerror=doseerror2, baseerror=baseerror2, allelicbias=allelicbias2, overdispersion=overdispersion2)       
        logl_off = sum_margloglls(loglls)  #   only account offpring data, excluding founders & midparents
        # calculate log for parents        
        logl_parent = cal_logl_parent(chrparentgeno,chrparent_hiddenunphased,
            chrparentformat,chrmultiallele,priorprocess,polyped;
            doseerror=doseerror2, baseerror=baseerror2, allelicbias=allelicbias2, overdispersion=overdispersion2)                   
        logl_off + logl_parent + logpri + logjacobi 
    end            
    if target == "doseerror"
        xstart = doseerror <= 1e-3 ? log(1e-3) : log(doseerror/(1-doseerror)) 
    elseif target == "baseerror"
        xstart = baseerror <= 1e-5 ? log(1e-5) : log(baseerror/(1-baseerror))    
    elseif target == "allelicbias"
        xstart = min(10.0,max(-10.0,log(allelicbias)))
    elseif target == "overdispersion"
        xstart = overdispersion <= 1e-3 ? log(1e-3) : max(-10,log(overdispersion))
    else
        @error string("unknown genoerror type: ",target)
    end 
    lowbound, upbound = max(-10.0,xstart-errorbound), min(10.0,xstart+errorbound)
    xstart = min(max(lowbound,xstart),upbound)
    x, fmax, his = PolyOrigin.brentMax(f,lowbound,upbound,xstart,
        precisiongoal=precisiongoal,accuracygoal=accuracygoal,maxiter=itmax)        
    est =  in(target, ["allelicbias","overdispersion"]) ? exp(x) : 1.0/(1.0+exp(-x)) 
    # @info "test" target his
    est    
end



function cal_logl_parent(chrparentgeno::AbstractMatrix, 
    chrparent_hiddenunphased::AbstractMatrix,
    chrparentformat::AbstractVector, 
    chrmultiallele::AbstractVector,            
    priorprocess::AbstractDict,
    polyped::PolyPed;
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector},  
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector},     
    snporder::Union{Nothing,AbstractVector}=nothing)
    logl_parent = 0.0
    pgenols = eachrow(chrparentgeno)
    hiddenls = eachrow(chrparent_hiddenunphased)
    ploidyls = polyped.parentinfo[!,:ploidy]     
    markerincl = first(values(priorprocess)).markerincl
    snps = isnothing(snporder) ? findall(markerincl) : collect(eachindex(markerincl))
    # eachindex(pgenols,hiddenls,chrparentformat,chrmultiallele)        
    for j in snps 
        geno_j = pgenols[j] 
        hidden_j = hiddenls[j]        
        logl_parent += sum(get_loglike_hidden.(geno_j, hidden_j, ploidyls; 
            format = chrparentformat[j],
            nallele = chrmultiallele[j],
            doseerror = isa(doseerror,AbstractVector) ? doseerror[j] : doseerror, 
            baseerror = isa(baseerror,AbstractVector) ? baseerror[j] : baseerror, 
            allelicbias = isa(allelicbias,AbstractVector) ? allelicbias[j] : allelicbias, 
            overdispersion = isa(overdispersion,AbstractVector) ? overdispersion[j] : overdispersion             
        ))

    end
    logl_parent
end



function infer_errorls!(chrparentgeno::Union{Nothing,AbstractMatrix}, # input parent geno
    chrparentformat::Union{Nothing,AbstractVector}, 
    chrphaplo::AbstractMatrix,
    derivedgeno::AbstractDict, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,    
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed, 
    targetls::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector,         
    overdispersionls::AbstractVector,
    errorprior::NamedTuple, 
    errorsoftthresh::NamedTuple,
    tukeyfence::Real, 
    errorhardthresh::NamedTuple, 
    isdellargeerr::Bool=true, 
    bindis::Real,
    errorbound::Real=5,
    errormaxit::Integer=20, 
    snporder::Union{Nothing,AbstractVector}=nothing,        
    temperature::Real=0.0)
    nsnp = size(chroffgeno,1)
    nsnp == length(doseerrorls) == length(baseerrorls) == length(allelicbiasls) || @error string("inconsistent #markers")
    isnothing(snporder) && (snporder = 1:nsnp)    
    loglls = []
    delmarkerls = []               
    parentploidyls = polyped.parentinfo[!,:ploidy]    
    chrparent_hiddenunphased = unphased_genocode_2_index.(join.(sort.(chrphaplo),"/"))
    for target in targetls        
        if target == "doseerror"          
            prior_err = errorprior.doseerror
        elseif target == "baseerror"
            prior_err = errorprior.baseerror
        elseif target == "allelicbias"        
            prior_err = errorprior.allelicbias
        elseif target == "overdispersion"
            prior_err = errorprior.overdispersion
        else
            @error string("unknown error target: ",target)
        end        
        # to exclude parents, revise:  cal_logbackward & get_bvdict & cal_sitedataprob_multiphase, 
        # length of startprobls or tranprobls is determined by bvdict        
        reslogbw = cal_logbackward(chrphaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
            bvpair, priorspace, priorprocess; 
            phaseerror_midparent,doseerror=doseerrorls, baseerror=baseerrorls, 
            allelicbias=allelicbiasls, overdispersion=overdispersionls, snporder,
            exclparent = true, 
        )
        # forward              
        markerincl=first(values(priorprocess)).markerincl
        tseq = findall(markerincl)    
        snp = snporder[tseq[1]]
        dataprobls = cal_sitedataprob_singlephase(snp,chrphaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],
            exclparent = true, 
        )
        bvdict = get_bvdict_exclparent(bvpair, priorspace)        
        startprobls=get_startprobls(priorprocess, bvdict)
        fwprob_logl = cal_initforward(startprobls, dataprobls)        
        for kk in eachindex(tseq)                        
            occursin("GT_", chroffformat[snporder[tseq[kk]]]) && in(target,["baseerror","allelicbias","overdispersion"]) && continue
            if kk >= 3
                snp = snporder[tseq[kk-1]]
                dataprobls = cal_sitedataprob_singlephase(snp,chrphaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
                    phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
                    allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
                    site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],
                    exclparent = true, 
                )
                tranprobls = get_tranprobls(priorprocess,bvdict,tseq[kk-2])
                # input fwprob_logl refers to tseq[kk-2], output fw_prob_logl refers to tseq[kk-1]
                cal_nextforward!(fwprob_logl, dataprobls, tranprobls)    
            end                        
            tpre = kk==1 ? 0 : tseq[kk-1]
            tnow = tseq[kk]
            logbw_tnow = reslogbw[tnow];
            snp = snporder[tnow]            
            errmarker = infer_errmarker(tpre, tnow, target, prior_err,          
                fwprob_logl, logbw_tnow, 
                chrparentgeno, chrparentformat, parentploidyls, 
                chrphaplo, chrparent_hiddenunphased, 
                derivedgeno,  chroffgeno,chroffformat, chrmultiallele, 
                bvdict, bvpair,priorspace,priorprocess; 
                phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
                allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
                snporder, temperature,errorbound, errormaxit
            )     
            if target == "doseerror"
                doseerrorls[snp] = errmarker                
            elseif target == "baseerror"
                baseerrorls[snp] = errmarker                            
            elseif target == "allelicbias"
                allelicbiasls[snp] = errmarker                
            elseif target == "overdispersion"
                overdispersionls[snp] = errmarker
            else
                @error string("unknown genoerror type: ",target)
            end 
        end
        # calculate logl
        kk = length(tseq) 
        snp = snporder[tseq[kk]]
        dataprobls = cal_sitedataprob_singlephase(snp,chrphaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion  = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],
            exclparent = true, 
        )
        tranprobls = get_tranprobls(priorprocess,bvdict,tseq[kk-1])        
        cal_nextforward!(fwprob_logl, dataprobls, tranprobls)    
        logl = sum(calindlogl(fwprob_logl,reslogbw, tseq[end]))
        push!(loglls,logl)
        # delete markers with too large error rates
        if !isinf(tukeyfence)
            snpincl = snporder[first(values(priorprocess)).markerincl]		        
            if target == "doseerror"
                errls = doseerrorls[snpincl] # either logit or log transformation does not work well for outlier detection of error rates.
                delsnps = get_outerlier_indices(errls; side="upper",tukeyfence)                                            
                delsnps = delsnps[errls[delsnps] .> errorsoftthresh.doseerror]            
                union!(delsnps,findall(errls .> errorhardthresh.doseerror))
                delsnps = snpincl[delsnps]
                # (isdellargeerr && !isempty(delsnps)) && println("large doseerror: ", round.(doseerrorls[delsnps],digits=3))                        
            elseif target == "baseerror"
                errls = baseerrorls[snpincl] # either logit or log transformation does not work well for outlier detection of error rates.
                delsnps = get_outerlier_indices(errls; side="upper",tukeyfence)
                delsnps = delsnps[errls[delsnps] .> errorsoftthresh.baseerror]
                union!(delsnps,findall(errls .> errorhardthresh.baseerror))
                delsnps = snpincl[delsnps]
                # (isdellargeerr && !isempty(delsnps)) && println("large baseerror: ", round.(baseerrorls[delsnps],digits=3))                        
            elseif target == "overdispersion"            
                errls = log.(overdispersionls[snpincl])
                delsnps = get_outerlier_indices(errls; side="upper",tukeyfence)
                delsnps = delsnps[errls[delsnps] .> log(errorsoftthresh.overdispersion)]
                union!(delsnps,findall(errls .> log(errorhardthresh.overdispersion)))
                delsnps = snpincl[delsnps]
                # (isdellargeerr && !isempty(delsnps)) && println("large overdispersion: ", round.(overdispersionls[delsnps],digits=2))                        
            else
                #  target == "allelicbias"            
                biasls = allelicbiasls[snpincl]
                side = "both"
                errls = log.(biasls)            
                outindices = get_outerlier_indices(errls; side,tukeyfence)
                absoft = abs(log(errorsoftthresh.allelicbias))
                delsnps_soft = findall(@. errls > absoft || errls < -1*absoft)
                delsnps = intersect(outindices, delsnps_soft) 
                abhard = abs(log(errorhardthresh.allelicbias))
                delsnps_hard = findall(@. errls > abhard || errls < -1*abhard)
                union!(delsnps, delsnps_hard) 
                delsnps = snpincl[delsnps]
                # (isdellargeerr && !isempty(delsnps)) && println("small/large bias: ", round.(allelicbiasls[delsnps],digits=3))                        
            end
            push!(delmarkerls, delsnps)
        end
    end  
    delmarkers = unique(reduce(vcat, delmarkerls;init=[]))
    ndel = length(delmarkers)
    # println("delsnps=",delsnps, ",err=", round.(doseerrorls[delsnps],digits=3))
    if !isempty(delmarkers) && isdellargeerr
        setpriorprocess!(priorprocess,snporder, delmarkers,bindis)
        checktranprobseq(priorprocess,bindis)
    end    
    loglls,ndel
end

function logit(x)
    log(x/(1-x))
end

function get_msg_err(targetls::AbstractVector, priorprocess::AbstractDict; 
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},    
    overdispersion::Union{Real,AbstractVector},
    snporder::Union{Nothing,AbstractVector}=nothing)
    msg = ""
    snpincl0 = first(values(priorprocess)).markerincl
    snpincl = isnothing(snporder) ? findall(snpincl0) : snporder[snpincl0]
    if in("doseerror",targetls)
        if isa(doseerror,AbstractVector)
            mdoseerr = round(mean(doseerror[snpincl]),digits=4)
            msg *= string(", <doseerror>=",mdoseerr)
        else
            msg *= string(", doseerror=",round(doseerror,digits=4))
        end
    end    
    if in("allelicbias",targetls)      
        if isa(allelicbias,AbstractVector)  
            mab = round(median(allelicbias[snpincl]),digits=2)
            msg *= string(", <allelicbias>=",mab)
        else
            msg *= string(", allelicbias=",round(allelicbias,digits=2))
        end
    end
    if in("overdispersion",targetls)      
        if isa(overdispersion,AbstractVector)  
            od = round(mean(overdispersion[snpincl]),digits=3)
            msg *= string(", <overdispersion>=",od)
        else
            msg *= string(", overdispersion=",round(overdispersion, digits=3))
        end
    end
    if in("baseerror",targetls)
        if isa(baseerror,AbstractVector)
            mseqerr = round(mean(baseerror[snpincl]),digits=5)
            msg *= string(", <baseerror>=",mseqerr)
        else
            msg *= string(", baseerror=",round(baseerror,digits=5))
        end
    end
    msg
end

function cal_logpri_logjacobi_err(target::AbstractString, prior_err::Distribution,x::Real)    
    if target in ["doseerror", "baseerror"]
        p =  1/(1+exp(-x)) 
        logpri = logpdf(prior_err,p)                        
        # estimation refers to transformed space
        logjacobi = x - 2*log(1+exp(x)) # dp/dx = exp(x)/(1+exp(x))^2                          
    else
        # target == "allelicbias", "overdispersion"    
        p = exp(x)
        logpri = logpdf(prior_err,p)                 
        logjacobi = x   # dp/dx = exp(x) # estimation refers to transformed space       
    end    
    # logjacobi = 0 # remove jacobian
    logpri, logjacobi
end


function infer_errmarker(tpre::Integer,tnow::Integer,
    target::AbstractString,  
    prior_err::Union{Nothing,Distribution},        
    fwprob_logl::AbstractVector,
    logbw_tnow::AbstractVector,  
    chrparentgeno::Union{Nothing,AbstractMatrix}, # input parent geno
    chrparentformat::Union{Nothing,AbstractVector}, 
    parentploidyls::AbstractVector,
    chrphaplo::AbstractMatrix, 
    chrparent_hiddenunphased::AbstractMatrix,         
    derivedgeno::AbstractDict,             
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    bvdict::AbstractDict, 
    bvpair::NamedTuple,        
    priorspace::AbstractDict,
    priorprocess::AbstractDict;    
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Real, 
    baseerror::Real, 
    allelicbias::Real,        
    overdispersion::Real, 
    errorbound::Real,
    errormaxit::Integer, 
    snporder::AbstractVector,
    temperature::Real=0.0)    
    bwprobls = [exp.(i .- maximum(i)) for i in logbw_tnow]
    fwprobls_tpre = fwprob_logl[1]
    fwprobls_tnow = get_fwprobls_tnow(tpre, fwprobls_tpre,bvdict, priorprocess)
    snp = snporder[tnow]    
    if !isnothing(chrparentgeno) && !isnothing(chrparentformat)
        pformat = chrparentformat[snp]
        pnallele = chrmultiallele[snp]    
        pgeno_snp = view(chrparentgeno, snp,:)
        phidden_snp = view(chrparent_hiddenunphased,snp,:)
    end
    function loglfun(x::Real)
        doseerror2, baseerror2, allelicbias2, overdispersion2 = doseerror, baseerror, allelicbias, overdispersion
        if target == "doseerror"
            doseerror2 = 1/(1+exp(-x))  # inverse of logit transformation     
        elseif target == "baseerror"
            baseerror2 = 1/(1+exp(-x))                              
        elseif target == "allelicbias"
            allelicbias2 = exp(x)        
        elseif target == "overdispersion"
            overdispersion2 = exp(x)
        else
            @error string("unknown genoerror type: ",target)
        end
        logpri, logjacobi = cal_logpri_logjacobi_err(target, prior_err, x)
        logl_off = calloglmarker(fwprobls_tnow, bwprobls, snp,chrphaplo, derivedgeno, bvpair,priorspace,
            chroffgeno,chroffformat, chrmultiallele, phaseerror_midparent,doseerror2, 
            baseerror2, allelicbias2, overdispersion2, 
        )         
        if isnothing(chrparentgeno) || isnothing(chrparentformat)
            logl_parent = 0
        else
            logl_parent = sum(get_loglike_hidden.(pgeno_snp, phidden_snp, parentploidyls; 
                format = pformat, nallele = pnallele, 
                doseerror=doseerror2, baseerror=baseerror2, allelicbias=allelicbias2, overdispersion=overdispersion2))  
        end
        logl_off + logl_parent + logpri + logjacobi        
    end    
    if target == "doseerror"
        xstart = doseerror <= 1e-3 ? log(1e-3) : log(doseerror/(1-doseerror)) 
    elseif target == "baseerror"
        xstart = baseerror <= 1e-5 ? log(1e-5) : log(baseerror/(1-baseerror))    
    elseif target == "allelicbias"
        xstart = min(10.0,max(-10.0,log(allelicbias)))
    elseif target == "overdispersion"
        xstart = overdispersion <= 1e-3 ? log(1e-3) : max(-10,log(overdispersion))
    else
        @error string("unknown genoerror type: ",target)
    end 
    lowbound, upbound = max(-10.0,xstart-errorbound), min(10.0,xstart+errorbound)
    xstart = min(max(lowbound,xstart),upbound)
    if temperature ≈ 0
        # always using brent
        accuracygoal, precisiongoal, maxiter = 3, 3, errormaxit
        res= brentMax(loglfun,lowbound,upbound, xstart; precisiongoal,accuracygoal,maxiter)            
        x = res[1]
    else
        constraint = x-> lowbound < x < upbound
        res2 = metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5)
        x = res2[end][1]        
    end
    est =  in(target, ["allelicbias", "overdispersion"]) ? exp(x) : 1.0/(1.0+exp(-x))     
    est
end



function calloglmarker(fwprobls_tnow::AbstractVector, 
    bwprobls::AbstractVector, 
    snp::Integer,            
    chrphaplo::AbstractMatrix, 
    derivedgeno::AbstractDict,     
    bvpair::NamedTuple,        
    priorspace::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    phaseerror_midparent::Union{Nothing,Real},
    doseerrormarker::Real,
    baseerrormarker::Real,
    allelicbiasmarker::Real, 
    overdispersionmarker::Real)
    dataprobls = cal_sitedataprob_singlephase(snp,chrphaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
        phaseerror_midparent,doseerror=doseerrormarker, baseerror = baseerrormarker, 
        allelicbias = allelicbiasmarker, overdispersion = overdispersionmarker,
        site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],
        exclparent = true, 
    )
    sum([log(sum(fwprobls_tnow[off] .* dataprobls[off] .* bwprobls[off])) for off in eachindex(fwprobls_tnow,dataprobls, bwprobls)])    
end

function get_fwprobls_tnow(tpre::Integer,fwprobls_tpre::AbstractVector,
    bvdict::AbstractDict, priorprocess::AbstractDict)
    if tpre == 0
        fwprobls_tnow=get_startprobls(priorprocess, bvdict)
    else
        fwprobls_tnow = Vector{Vector{Float64}}(undef, length(fwprobls_tpre))
        tranprobls = get_tranprobls(priorprocess,bvdict,tpre)
        for off in eachindex(fwprobls_tnow, tranprobls)
            tranprob = tranprobls[off]    
            if any(isnothing.(tranprob))
                fwprobls_tnow[off] = copy(fwprobls_tpre[off])
            else
                fwprobls_tnow[off] = kronvec(tranprob[1]',tranprob[2]',fwprobls_tpre[off])        
            end
        end
    end
    fwprobls_tnow
end

