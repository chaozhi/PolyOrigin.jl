
function updateorder!(chrparentgeno::AbstractMatrix, 
    chroffgeno::AbstractMatrix, 
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    bvpair::NamedTuple,
    priorspace::AbstractDict, 
    priorprocess::AbstractDict;
    snporder::AbstractVector, 
    phaseerror_midparent::Union{Nothing,Real},
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector, 
    overdispersionls::AbstractVector,
    isreversechr::Bool=false,
    maxwinsize::Integer=50,
    temperature::Real=0)
    if isreversechr
        reverse!(snporder)        
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
            chrparentgeno, chroffgeno, chroffformat, chrmultiallele, 
            derivedgeno, bvpair, priorspace, priorprocess,
            phaseerror_midparent, snporder, doseerrorls, baseerrorls, allelicbiasls,overdispersionls)
        accept = mean(accepthis[1:end-winsize+1])
        deltlogl = loglhis[end]-loglhis[1]
        temp = vcat(round.([loglhis[end], deltlogl],digits=2),[round(accept,digits=4)])
        push!(logbook,vcat([string("#reverse_w",winsize)],temp))                
        if accept >= 0.01 && abs(deltlogl)>1 && winsize<min(maxwinsize,nsnp/2)
            winsize += 1
        else
            break
        end
    end
    if isreversechr
        reverse!(snporder)          
        for (key, val) in priorprocess
             reverseprior!(val)
        end
    end
    snporder,priorprocess,logbook
end


function reverseorder!(winsize::Integer, temperature::Real,    
    chrparentgeno::AbstractMatrix, 
    chroffgeno::AbstractMatrix, 
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict,
    bvpair::NamedTuple,
    priorspace::AbstractDict, 
    priorprocess::AbstractDict, 
    phaseerror_midparent::Union{Nothing,Real},
    snporder::AbstractVector,     
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector,
    overdispersionls::AbstractVector) 
    reslogbw = cal_logbackward(chrparentgeno, chroffgeno,chroffformat, chrmultiallele, derivedgeno, 
        bvpair, priorspace, priorprocess; 
        doseerror=doseerrorls, baseerror=baseerrorls, 
        allelicbias=allelicbiasls, overdispersion=overdispersionls,
        phaseerror_midparent, snporder
    )
    pri1 = first(values(priorprocess))
    tseq = findall(pri1.markerincl)
    snp = snporder[tseq[1]]
    dataprobls = cal_sitedataprob_singlephase(snp,chrparentgeno, derivedgeno, bvpair,priorspace,chroffgeno; 
        doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
        allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
        site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],phaseerror_midparent, 
    )
    bvdict = get_bvdict(bvpair, priorspace)        
    startprobls=get_startprobls(priorprocess, bvdict)
    fwprob_logl = cal_initforward(startprobls, dataprobls)        
    logl = sum(calindlogl(fwprob_logl,reslogbw, tseq[1]))    
    # println("winsize=", winsize, ", initial logl=",logl)
    loglhis = Vector{Float64}()
    push!(loglhis,logl)
    accepthis = Vector{Bool}()
    for kk=1:length(tseq)-1        
        kkmin=kk
        kkmax = min(kk+winsize-1,length(tseq))
        fwseg, fwseg_reverse = [cal_segforward(fwprob_logl, isreverse, tseq , kkmin, kkmax, bvdict, 
            chrparentgeno, chroffgeno, chroffformat, chrmultiallele, derivedgeno, bvpair, 
            priorspace, priorprocess, phaseerror_midparent,snporder, 
            doseerrorls, baseerrorls, allelicbiasls,overdispersionls) for isreverse in [false, true]]        
        logl = sum(calindlogl(last(fwseg),reslogbw, tseq[kkmax]))
        proplogl = sum(calindlogl(last(fwseg_reverse),reslogbw, tseq[kkmax]))
        tseg = tseq[kkmin:kkmax]
        if temperature <=0
            isaccept = proplogl >= logl
        else
            isaccept = rand() < min(1,exp((proplogl - logl)/temperature))
        end
        # println("[kkmin,kkmax]=",[kkmin,kkmax], ", logl=",logl, ", proplogl=",proplogl, ", isaccept=",isaccept)
        if isaccept            
            snporder[tseg] .= snporder[reverse(tseg)]
            set_reverse_seg!(priorprocess,tseg)
            logl = proplogl        
            fwprob_logl = first(fwseg_reverse)
        else
            fwprob_logl = first(fwseg)
        end        
        push!(accepthis,isaccept)
        push!(loglhis,logl)
    end
    snporder, priorprocess, loglhis, accepthis
end

function set_reverse_seg!(priorprocess::AbstractDict,tt::AbstractVector)
    tt2=tt[1:end-1]
    for (strkey, pri) in priorprocess
        pri.markerdeltd[tt2] = pri.markerdeltd[reverse(tt2)]
        pri.tranprobseq[tt2] = pri.tranprobseq[reverse(tt2)]
        pri.markerid[tt] = pri.markerid[reverse(tt)]
    end
    priorprocess
end



function cal_segforward(fwprob_logl::AbstractVector, 
    isreverse::Bool, 
    tseq::AbstractVector, 
    kkmin::Integer,
    kkmax::Integer,    
    bvdict::AbstractDict, 
    chrparentgeno::AbstractMatrix,
    chroffgeno::AbstractMatrix,
    chroffformat::AbstractVector, 
    chrmultiallele::AbstractVector,
    derivedgeno::AbstractDict, 
    bvpair::NamedTuple,
    priorspace::AbstractDict, 
    priorprocess::AbstractDict, 
    phaseerror_midparent::Union{Nothing,Real},
    snporder::AbstractVector, 
    doseerrorls::AbstractVector, 
    baseerrorls::AbstractVector, 
    allelicbiasls::AbstractVector,
    overdispersionls::AbstractVector) 
    # input fwprob_logl refers to kkmin-1 if kkmin > 1, and otherwise it refers to kkmn = 1    
    tls = tseq[kkmin:kkmax]
    tls_tran = tseq[kkmin:kkmax-1]
    if isreverse 
        reverse!(tls)
        reverse!(tls_tran)
    end
    if kkmin == 1
        t1 = popfirst!(tls)
        if isreverse                
            # recalculate fwprob_logl
            snp = snporder[t1]
            dataprobls = cal_sitedataprob_singlephase(snp,chrparentgeno, derivedgeno, bvpair,priorspace,chroffgeno; 
                doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
                allelicbias = allelicbiasls[snp], overdispersion =  overdispersionls[snp],
                site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],phaseerror_midparent, 
            )
            startprobls=get_startprobls(priorprocess, bvdict)
            fwprob_logl_last = cal_initforward(startprobls, dataprobls)                    
            fwprob_logl_next = copy_fwprob_logl(fwprob_logl_last) 
        else            
            fwprob_logl_last = copy_fwprob_logl(fwprob_logl)             
            fwprob_logl_next = copy_fwprob_logl(fwprob_logl) 
        end        
    else
        pushfirst!(tls_tran, tseq[kkmin-1])
        fwprob_logl_last = copy_fwprob_logl(fwprob_logl) 
        fwprob_logl_next = similar(fwprob_logl) 
    end    
    for i in eachindex(tls, tls_tran)        
        snp = snporder[tls[i]]
        dataprobls = cal_sitedataprob_singlephase(snp,chrparentgeno, derivedgeno, bvpair,priorspace,chroffgeno; 
            doseerror=doseerrorls[snp], baseerror = baseerrorls[snp], 
            allelicbias = allelicbiasls[snp], overdispersion = overdispersionls[snp],
            site_offformat = chroffformat[snp], site_nallele = chrmultiallele[snp],phaseerror_midparent, 
        )
        tranprobls = get_tranprobls(priorprocess,bvdict,tls_tran[i])                
        cal_nextforward!(fwprob_logl_last, dataprobls, tranprobls)    
        if i == 1 && kkmin > 1
            fwprob_logl_next .= copy_fwprob_logl(fwprob_logl_last) 
        end
    end    
    fwprob_logl_next, fwprob_logl_last
end


function copy_fwprob_logl(fwprob_logl::AbstractVector)
    [copy.(fwprob_logl[1]), copy(fwprob_logl[2])]     
end