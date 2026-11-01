
function stripchrend!(priorprocess;stripdis::Real=20)
    pri1=first(values(priorprocess))
    ndel = 0
    deltd = pri1.markerdeltd[1:end-1]
    nend = max(1,round(Int,0.05*length(deltd)))
    ii=findfirst(x-> !isnothing(x) && (x>stripdis/100),deltd[1:nend])
    if !isnothing(ii)
        for (strkey,pri) in priorprocess            
            ndel += sum(pri.markerincl[1:ii])
            pri.markerincl[1:ii] .= 0
            pri.markerdeltd[1:ii] .= nothing
        end
    end
    ii=findfirst(x-> !isnothing(x) && (x>stripdis/100),reverse(deltd[end-nend:end]))    
    if !isnothing(ii)
        for (strkey,pri) in priorprocess
            ndel += sum(pri.markerincl[end-ii+1:end])
            pri.markerincl[end-ii+1:end] .= 0
            pri.markerdeltd[end-ii:end] .= nothing
        end
    end
    ndel
end


function updatedistance!(snporder::AbstractVector,
    bvpair::NamedTuple,priorspace::AbstractDict,
    priorprocess::AbstractDict,
    chrparentgeno::AbstractMatrix, 
    chroffgeno::AbstractMatrix, 
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector, 
    derivedgeno::AbstractDict;
    phaseerror_midparent::Union{Nothing,Real},
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector, 
    overdispersionls::AbstractVector,
    priorlength::Real=100,     
    bindis::Real,
    temperature::Real=0.0)    
    # backward
    reslogbw = cal_logbackward(chrparentgeno, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        bvpair, priorspace, priorprocess; 
        phaseerror_midparent,doseerror=doseerrorls, baseerror=baseerrorls, 
        allelicbias=allelicbiasls, overdispersion=overdispersionls, snporder
    )
    # forward    
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)    
    snp = snporder[tseq[1]]
    dataprobls = cal_sitedataprob_singlephase(snp,chrparentgeno, derivedgeno, bvpair,priorspace,chroffgeno; 
        phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
        allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
        site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
    )
    bvdict = get_bvdict(bvpair, priorspace)
    startprobls=get_startprobls(priorprocess, bvdict)
    fwprob_logl = cal_initforward(startprobls, dataprobls)
    kkmax = length(tseq)-1
    for kk=1:kkmax
        snp = snporder[tseq[kk+1]]
        dataprobls = cal_sitedataprob_singlephase(snp,chrparentgeno, derivedgeno, bvpair,priorspace,chroffgeno; 
            phaseerror_midparent,doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp]
        )
        logbwprobls = reslogbw[tseq[kk+1]]
        tnowdis = calinterdis(tseq[kk],priorprocess,bvdict, dataprobls,fwprob_logl[1], logbwprobls,priorlength,bindis, temperature)
        setdistanceat!(priorprocess,tseq[kk],tnowdis,bindis)
        b = checktranprobseq(priorprocess,bindis)
        if !b
            println("kk=",kk, ",t=",tseq[kk], ",d=",tnowdis)
        end        
        tranprobls = get_tranprobls(priorprocess,bvdict,tseq[kk])
        cal_nextforward!(fwprob_logl, dataprobls, tranprobls)        
    end
    logl = sum(calindlogl(fwprob_logl,reslogbw, tseq[end]))    
    logl
end

function callogldis(d::Real,priorprocess::AbstractDict, bvdict::AbstractDict, dataprobls::AbstractVector,
    fwprobls::AbstractVector,bwprobls::AbstractVector,bindis::Real)    
    tranprobls = cal_tranprobls(priorprocess, bvdict,d,bindis)
    s = 0.0 
    for off=eachindex(fwprobls,bwprobls, dataprobls, tranprobls)
        tranls = tranprobls[off]
        if any(isnothing.(tranls))
            ls = fwprobls[off]
        else
            ls = kronvec(tranprobls[off]...,fwprobls[off])
        end
        s += log(sum(ls .* dataprobls[off] .* bwprobls[off])) 
    end
    s
end

function calinterdis(tnow::Integer,
    priorprocess::AbstractDict,
    bvdict::AbstractDict, 
    dataprobls::AbstractVector, # refers to tnow+1
    fwprobls::AbstractVector, # refers to tnow
    logbwprobls::AbstractVector, # refers to tnow+1    
    priorlength::Real,
    bindis::Real,
    temperature::Real)
    # Haldane map function, priorlength in unit of cM    
    markerdeltd = first(values(priorprocess)).markerdeltd
    accuracygoal, precisiongoal, maxiter = 5, 5, 20
    lowbound,upbound = log(10^(-7)), log(10.0)    
    bwprobls = [exp.(logbw .- maximum(logbw)) for logbw in logbwprobls]
    function loglfun(x::Real)
        d = exp(x)
        prior =  - d*100/priorlength        
        callogldis(d,priorprocess,bvdict, dataprobls,fwprobls,bwprobls,bindis) + prior
    end    
    xstart = log(max(exp(lowbound),markerdeltd[tnow]))
    if temperature ≈ 0
        # always using brent        
        res=brentMax(loglfun,lowbound,upbound, xstart; precisiongoal, accuracygoal,maxiter)        
        x = res[1]
    else
        constraint = x-> lowbound < x < upbound
        res2 = metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5)
        x = res2[end][1]        
    end
    d = round(exp(x),digits=6)        
    d < bindis/100 ? 0.0 : d
end
