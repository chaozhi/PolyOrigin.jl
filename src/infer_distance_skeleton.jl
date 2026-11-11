
function infer_skeletonmap(chrparentgeno::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict, 
    bvpair::NamedTuple,
    priorspace::AbstractDict,
    priorprocess::AbstractDict,
    polyped::PolyPed,
    targetls::AbstractVector,
    phaseerror_midparent::Union{Nothing,Real},
    isinfererror::Bool, 
    inputdoseerrorls::AbstractVector, 
    inputbaseerrorls::AbstractVector, 
    inputallelicbiasls::AbstractVector,         
    inputoverdispersionls::AbstractVector,    
    errorprior::NamedTuple, 
    errorsoftthresh::NamedTuple,    
    tukeyfence::Real,
    errorhardthresh::NamedTuple; 
    priorlength::Real=100, 
    bindis::Real,
    skeletonsize::Union{Nothing,Integer}, 
    snporder::Union{Nothing,AbstractVector}=nothing,        
    chrid::AbstractString, 
    startit::Integer=1, 
    io::Union{Nothing,IO},
    verbose::Bool=true)    
    errlsls = copy.([inputdoseerrorls, inputbaseerrorls, inputallelicbiasls, inputoverdispersionls])    
    skeleton, nbin = getskeleton_seg(errlsls,snporder,priorprocess,skeletonsize)
    skeletonprior = deepcopy(priorprocess)
    pri1=first(values(skeletonprior))
    if 0 in pri1.markerincl[skeleton]
        @error string("deleted markers are not allowed in skeleton")
    end
    deltt = trues(length(pri1.markerincl))
    deltt[skeleton] .= false
    setpriorprocess!(skeletonprior,deltt,bindis)
    pri1=first(values(skeletonprior))
    exclsnps = snporder[.!pri1.markerincl]
    errlsls = copy.([inputdoseerrorls, inputbaseerrorls, inputallelicbiasls, inputoverdispersionls])
    for errls in errlsls
        errls[exclsnps] .= NaN
    end
    doseerrorls, baseerrorls, allelicbiasls,overdispersionls = errlsls
    chrlenls = []
    for it= startit:startit+10
        startt = time()        
        # msg = string("chr=", chrid,", it=", it, ", #skeleton=",length(skeleton),", #bin=",nbin)
        msg = string("chr=", chrid,", it=", it, ", #skeleton=",length(skeleton))
        if isinfererror
            infer_errorls!(nothing, nothing, chrparentgeno, derivedgeno, chroffgeno, chroffformat, chrmultiallele, 
                bvpair, priorspace, skeletonprior, polyped, targetls;     
                phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, errorprior, 
                errorsoftthresh,tukeyfence, errorhardthresh, snporder, bindis, 
            )        
            msg *= get_msg_err(targetls, skeletonprior; snporder, 
                doseerror=doseerrorls, baseerror=baseerrorls, 
                allelicbias=allelicbiasls, overdispersion=overdispersionls)                
        end        
        updatedistance!(snporder,bvpair,priorspace,skeletonprior,
            chrparentgeno, chroffgeno, chroffformat, chrmultiallele, derivedgeno; 
            phaseerror_midparent,doseerrorls, baseerrorls, allelicbiasls, overdispersionls, priorlength, bindis
        )          
        chrlen = round(get_chrlen(skeletonprior),digits=1)                                
        msg *= string(", len=",chrlen, "cM",                        
            ", tused=", round(time()-startt, digits=1), "s"
        )
        printconsole(io,verbose,msg)
        if length(chrlenls) >= 3 
            diff_chrlen = max(abs.(chrlen .- chrlenls[end-2:end])...)
            if diff_chrlen<=1 || diff_chrlen/chrlen < 0.005
                break
            end
        end
        push!(chrlenls,chrlen)
    end
    skeletonprior
end


function getskeleton_seg(errlsls::AbstractVector,snporder::AbstractVector,
    priorprocess::AbstractDict,skeletonsize::Union{Nothing,Integer})
    doseerrorls, baseerrorls, allelicbiasls, allelicoverdispersionls, = errlsls
    tterrorls=[doseerrorls[i] + baseerrorls[i] + 0.01*abs(log10(allelicbiasls[i])) + 0.01*allelicoverdispersionls[i] for i in snporder]
    pri1=first(values(priorprocess))
    excl = findall(.!pri1.markerincl)
    deltd= copy(pri1.markerdeltd)
    deltd[isnothing.(deltd)] .= 0.0
    loc =accumulate(+,deltd[1:end-1])
    pushfirst!(loc,0)
    bins = splitindex(loc)
    binls = [setdiff(collect(i),excl) for i=bins]
    binls = binls[length.(binls) .> 0]
    nbin = length(binls)
    skeleton=[bin[argmin(tterrorls[bin])] for bin=binls]    
    if isnothing(skeletonsize) || length(skeleton)<=skeletonsize
        incl=findall(pri1.markerincl)
        incl[1] in skeleton || pushfirst!(skeleton,incl[1])
        incl[end] in skeleton || push!(skeleton,incl[end])
        return skeleton,nbin
    end    
    skeldeltd=deltd[[bin[end] for bin=binls]]
    skeldeltd[end]==0 || @error "unexpected last skeleton deltd"
    isapprox(sum(skeldeltd),sum(deltd);atol=1e-2) || @error string("unequal sum skeleton deltd, 
        sum(skeldeltd)=",sum(skeldeltd), ", sum(deltd)=",sum(deltd))
    skelloc=accumulate(+,circshift(skeldeltd,1))
    gridsize=skelloc[end]/(skeletonsize-2)
    grid=gridpartition(skelloc,gridsize)
    grid=grid[length.(grid) .> 0]
    grid=[skeleton[i] for i=grid]
    skeleton2=[bin[argmin(tterrorls[bin])] for bin=grid]
    incl=findall(pri1.markerincl)
    incl[1] in skeleton2 || pushfirst!(skeleton2,incl[1])
    incl[end] in skeleton2 || push!(skeleton2,incl[end])
    ndiff = skeletonsize-length(skeleton2)
    if ndiff>0
        ii=setdiff(skeleton,skeleton2)
        s=sortperm(tterrorls[ii])
        length(s) > ndiff && (s=s[1:ndiff])
        skeleton2=sort(vcat(skeleton2,ii[s]))
    end
    skeleton2,nbin
end

function rescale_map!(priorprocess::AbstractDict,skeletonprior::AbstractDict,bindis::Real)
    pri0=first(values(priorprocess))
    incl0 = findall(pri0.markerincl)
    deltd0=deepcopy(pri0.markerdeltd)
    pri1=first(values(skeletonprior))
    incl1 = findall(pri1.markerincl)
    deltd1=pri1.markerdeltd
    issubset(incl1,incl0) || @warn string("unexpected skeleton incl")
    if !in(incl0[1], incl1)
        dlshead = deltd0[1:incl0[1]-1]
        b = [isnothing(i) || i ≈ 0.0 for i in dlshead]
        all(b) || @warn string("first marker is not in skeleton")
    end
    if !in(incl0[end], incl1)
        dlstail = deltd0[incl0[end]:end]
        b = [isnothing(i) || i ≈ 0.0 for i in dlstail]
        all(b) || @warn string("last marker is not in skeleton")
    end
    for i=1:length(incl1)-1
        s=incl1[i]:incl1[i+1]-1        
        d0 = deltd0[s]
        d1 = deltd1[s]
        temp=sum(d0[.!isnothing.(d0)])
        scale = temp ≈ 0 ? 0 : sum(d1[.!isnothing.(d1)])/temp        
        deltd0[s] .= [isnothing(j) ? j : round(j*scale,digits=6) for j in deltd0[s]]
    end    
    for (_, pri) in priorprocess
        pri.markerdeltd = deltd0                
        pri.tranprobseq = [isnothing(i) ? i : getgametetran(i,pri1.nvalent,bindis) for i=deltd0]
    end
    priorprocess
end
