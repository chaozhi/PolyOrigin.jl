
function cal_singlelogl(phaplo::AbstractMatrix, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    priorspace::AbstractDict,priorprocess::AbstractDict,
    bvpair::NamedTuple; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},     
    overdispersion::Union{Real,AbstractVector},
    markerls::Union{Nothing, AbstractVector} = nothing)
    nmarker = size(chroffgeno,1)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nmarker) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nmarker) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nmarker) : allelicbias      
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nmarker) : overdispersion
    bvdict = get_bvdict(bvpair, priorspace)
    startprobls = get_startprobls(priorprocess,bvdict)
    if isnothing(markerls)  
        markerls = 1:nmarker    
    end    
    singlelogl = zeros(length(startprobls),nmarker)
    for snp in markerls    
        dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent, doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
        )
        singlelogl[:, snp] .= map((x,y) -> log(dot(x,y)), startprobls, dataprobls)
    end 
    singlelogl
end

function cal_logbackward(phaplo::AbstractMatrix, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,    
    derivedgeno::AbstractDict,
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector}, 
    overdispersion::Union{Real,AbstractVector}, 
    exclparent::Bool=false,
    snporder::Union{Nothing,AbstractVector}=nothing)
    nseq = size(chroffgeno,1)
    isnothing(snporder) && (snporder = 1:nseq)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nseq) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nseq) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nseq) : allelicbias        
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nseq) : overdispersion
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)   
    reslogbw = Vector{Any}(undef, length(markerincl))
    bvdict = exclparent ? get_bvdict_exclparent(bvpair, priorspace) : get_bvdict(bvpair, priorspace)               
    reslogbw[tseq[end]] = [zeros(length(i)) for i in get_startprobls(priorprocess, bvdict)]            
    for kk in length(tseq)-1:-1:1                
        tnow = tseq[kk]            
        snp = snporder[tseq[kk+1]]
        dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],             
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],
            exclparent
        )
        tranprobls = get_tranprobls(priorprocess, bvdict,tnow);
        logbwprob = deepcopy(reslogbw[tseq[kk+1]])
        cal_nextlogbackward!(logbwprob, dataprobls, tranprobls)                    
        reslogbw[tnow] = logbwprob
    end    
    reslogbw
end

# input logbwprob refers to tback, and output logbwprob refers to tnow
# dataprobls refers to tback not tnow!, and tranprobls gives tranprob from tnow to tback for each individual
# tnow < tback backwardly
function cal_nextlogbackward!(logbwprob::AbstractVector, dataprobls::AbstractVector, tranprobls::AbstractVector)        
    for ind in eachindex(logbwprob, dataprobls, tranprobls)
        bw = logbwprob[ind]
        logpmax = maximum(bw)
        @. bw = exp(bw - logpmax) * dataprobls[ind] 
        tranprob = tranprobls[ind]
        if !any(isnothing.(tranprob))
            bw .= kronvec(tranprob[1],tranprob[2], bw) 
        end
        @. bw = log(bw) + logpmax    
    end
    logbwprob
end

function cal_initforward(startprobls::AbstractVector, dataprobls::AbstractVector)
    fwprob = [similar(i) for i in dataprobls]
    fwlogl = zeros(length(dataprobls))
    for ind in eachindex(startprobls, dataprobls)
        fwprob[ind] .= startprobls[ind] .* dataprobls[ind]
        scale = sum(fwprob[ind])
        fwprob[ind] ./= scale
        fwlogl[ind] = log(scale)
    end
    [fwprob, fwlogl]
end

function cal_nextforward!(fwprob_logl::AbstractVector, dataprobls::AbstractVector, tranprobls::AbstractVector)
    fwprob, fwlogl = fwprob_logl
    for ind in eachindex(fwprob, fwlogl, dataprobls,tranprobls)
        dataprob = dataprobls[ind]
        tranprob = tranprobls[ind]
        if any(isnothing.(tranprob))
            fwprob[ind] .= fwprob[ind] .* dataprob
        else
            # fwprob[ind] .= (tranprob' * fwprob[ind]) .* dataprob
            fwprob[ind] .= kronvec(tranprob[1]',tranprob[2]',fwprob[ind]) .* dataprob
        end
        scale = sum(fwprob[ind])
        fwprob[ind] ./= scale
        fwlogl[ind] += log(scale)
    end
    fwprob_logl
end

function calindlogl(fwprob_logl::AbstractVector,logbwprob::AbstractVector)
    fwprob, fwlogl = fwprob_logl
    [begin 
        lpmax = maximum(logbwprob[ind])            
        lp = exp.(logbwprob[ind] .- lpmax)             
        log(dot(lp, fwprob[ind])) + lpmax + fwlogl[ind]
    end for ind in eachindex(logbwprob, fwprob, fwlogl)]
end

function calindlogl(fwprob_logl::AbstractVector,reslogbw::AbstractVector, t::Integer)
    logbwprob = reslogbw[t] 
    calindlogl(fwprob_logl,logbwprob)
end

function calvuongts(logldiff::AbstractVector,dfdiff::Integer)
    n = length(logldiff)
    v = var(logldiff)
    vuong = sum(logldiff)
    # correction for the difference of number of parameters
    vuong -= log(n)*(dfdiff)/2
    vuong /= sqrt(n*v)
    vuong
end

function calvuongts!(kk::Integer,
    tseq::AbstractVector,
    logl::AbstractVector,    
    fwprob_logl::AbstractVector,    
    singlelogl::AbstractMatrix,
    reslogbw::AbstractVector, 
    bvdict::AbstractDict, 
    phaplo::AbstractMatrix, 
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector,
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    bvpair::NamedTuple; 
    phaseerror_midparent::Union{Nothing,Real},
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector, 
    overdispersionls::AbstractVector, 
    bindis::Real, 
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))
    # fwdict is modified inside of function
    dfdiff = -1
    isologl = singlelogl[:, snporder[tseq[kk]]]
    if kk==1
        snp  = snporder[tseq[2]]
        dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
        )
        startprobls=get_startprobls(priorprocess, bvdict)
        init_prob_logl = cal_initforward(startprobls, dataprobls)
        proplogl = calindlogl(init_prob_logl,reslogbw, tseq[2])
        proplogl .+= isologl
    elseif kk==length(tseq)        
        proplogl = fwprob_logl[2]        
        proplogl .+= isologl
    else
        snp  = snporder[tseq[kk+1]]
        dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
        )
        pri1=first(values(priorprocess))
        tnowdis=sum(pri1.markerdeltd[tseq[kk-1:kk]])
        tranprobls = cal_tranprobls(priorprocess,bvdict,tnowdis,bindis)
        prop_prob_logl = deepcopy(fwprob_logl) # fw_prob_logl refers to tseq[kk-1]
        cal_nextforward!(prop_prob_logl, dataprobls, tranprobls)
        proplogl = calindlogl(prop_prob_logl,reslogbw, tseq[kk+1])      
        proplogl .+= isologl 
        # updating  fwprob_logl for next iteration       
        snp  = snporder[tseq[kk]]
        dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
        )
        tranprobls = get_tranprobls(priorprocess,bvdict,tseq[kk-1])
        cal_nextforward!(fwprob_logl, dataprobls, tranprobls)        
    end
    logldiff = proplogl .- logl
    isempty(logldiff) ? 0 : calvuongts(logldiff,dfdiff)
end

function polymarkerdel!(bvpair::NamedTuple, 
    phaploindex::AbstractVector,
    phaploset::AbstractVector,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    delsiglevel::Real, 
    bindis::Real,
    snporder::Union{Nothing,AbstractVector}=nothing)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))
    phaplo=get_phaplo(phaploindex,phaploset)
    derivedgeno = get_derivedgeno(phaplo,chrmultiallele, priorspace)
    singlelogl = cal_singlelogl(phaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        priorspace, priorprocess, bvpair; phaseerror_midparent,doseerror, 
        baseerror, allelicbias,overdispersion
    )    
    reslogbw = cal_logbackward(phaplo, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        bvpair, priorspace, priorprocess; 
        phaseerror_midparent,doseerror, baseerror, allelicbias, overdispersion, snporder
    )
    # 
    nmarker = size(chroffgeno,1)
    doseerrorls = typeof(doseerror) <: Real ? doseerror*ones(nmarker) : doseerror
    baseerrorls = typeof(baseerror) <: Real ? baseerror*ones(nmarker) : baseerror
    allelicbiasls = typeof(allelicbias) <: Real ? allelicbias*ones(nmarker) : allelicbias       
    overdispersionls = typeof(overdispersion) <: Real ? overdispersion*ones(nmarker) : overdispersion
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)    
    snp = snporder[tseq[1]]
    dataprobls = cal_sitedataprob_singlephase(snp,phaplo, derivedgeno, bvpair,priorspace,chroffgeno; 
        phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
        allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
        site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
    )
    bvdict = get_bvdict(bvpair, priorspace)
    startprobls=get_startprobls(priorprocess, bvdict)
    fwprob_logl = cal_initforward(startprobls, dataprobls)
    logl = calindlogl(fwprob_logl,reslogbw, tseq[1])
    res=zeros(length(tseq))
    for kk=1:length(tseq)
        res[kk] = calvuongts!(kk,tseq,logl,fwprob_logl, singlelogl, reslogbw, bvdict, 
            phaplo, chroffgeno, chroffformat, chrmultiallele, derivedgeno, 
            priorspace,priorprocess,bvpair; 
            phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, bindis, snporder)
    end 
    threshold = abs(quantile(Normal(),delsiglevel))
    deltt = tseq[res .> threshold]    
    setpriorprocess!(priorprocess, deltt,bindis)
    deltt
end

function polymarkerdel!(bvpair::NamedTuple,    
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrparentgeno::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector;
    phaseerror_midparent::Union{Nothing,Real},
    doseerror::Union{Real,AbstractVector}, 
    baseerror::Union{Real,AbstractVector}, 
    allelicbias::Union{Real,AbstractVector},
    overdispersion::Union{Real,AbstractVector},
    delsiglevel::Real, 
    bindis::Real,   
    snporder::Union{Nothing,AbstractVector}=nothing)    
    phaploset=[[[i] for i=j] for j=eachcol(chrparentgeno)]
    phaploindex=[ones(Int, length(j)) for j=phaploset]
    polymarkerdel!(bvpair, phaploindex,phaploset, priorspace, priorprocess,
        chroffgeno, chroffformat, chrmultiallele; 
        phaseerror_midparent, doseerror, baseerror, allelicbias, overdispersion, delsiglevel, snporder,bindis
    )   
end
